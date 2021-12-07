/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2014 by Wenzel Jakob and others.

    Mitsuba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Mitsuba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <mitsuba/core/platform.h>

// Mitsuba's "Assert" macro conflicts with Xerces' XSerializeEngine::Assert(...).
// This becomes a problem when using a PCH which contains mitsuba/core/logger.h
#if defined(Assert)
# undef Assert
#endif
#include <xercesc/parsers/SAXParser.hpp>
#include <mitsuba/core/sched_remote.h>
#include <mitsuba/core/sstream.h>
#include <mitsuba/core/fresolver.h>
#include <mitsuba/core/fstream.h>
#include <mitsuba/core/appender.h>
#include <mitsuba/core/sshstream.h>
#include <mitsuba/core/shvector.h>
#include <mitsuba/core/statistics.h>
#include <mitsuba/render/renderjob.h>
#include <mitsuba/render/scenehandler.h>
#include <fstream>
#include <stdexcept>
#include <boost/algorithm/string.hpp>


#include <mitsuba/core/point.h>
#include <mitsuba/core/fwd.h>
#include <mitsuba/core/vector.h>

#include <mitsuba/core/fresolver.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/hw/basicshader.h>
#include "../src/bsdfs/ior.h"


#if defined(__WINDOWS__)
#include <mitsuba/core/getopt.h>
#include <winsock2.h>
#else
#include <signal.h>
#endif
using namespace std;
#include <iostream>
using std::cerr;
using std::endl;
#include <fstream>
using std::ofstream;
#include <cstdlib> // for exit function

//#include "../src/bsdfs/roughGGX.h"
using XERCES_CPP_NAMESPACE::SAXParser;

using namespace mitsuba;

// I have imported all the RoughGGX.h functions for "include" problems
static inline float generateRandomNumber()
{
	const float U = ((float)rand()) / (float)RAND_MAX;
	return U;
}

static inline bool IsFiniteNumber(float x)
{
	return (x <= std::numeric_limits<float>::max() && x >= -std::numeric_limits<float>::max());
}


static inline double  abgam(double x)
{
	double  gam[10],
		temp;

	gam[0] = 1. / 12.;
	gam[1] = 1. / 30.;
	gam[2] = 53. / 210.;
	gam[3] = 195. / 371.;
	gam[4] = 22999. / 22737.;
	gam[5] = 29944523. / 19733142.;
	gam[6] = 109535241009. / 48264275462.;
	temp = 0.5 * log(2 * M_PI) - x + (x - 0.5) * log(x)
		+ gam[0] / (x + gam[1] / (x + gam[2] / (x + gam[3] / (x + gam[4] /
			(x + gam[5] / (x + gam[6] / x))))));

	return temp;
}

inline double  gamma(double x)
{
	double  result;
	result = exp(abgam(x + 5)) / (x * (x + 1) * (x + 2) * (x + 3) * (x + 4));
	return result;
}

static inline double  beta(double m, double n)
{
	return (gamma(m) * gamma(n) / gamma(m + n));
}

#define vec3 Vector
#define vec2 Vector

struct RayInfo
{
	// direction
	vec3 w;
	float theta;
	float cosTheta;
	float sinTheta;
	float tanTheta;
	float alpha;
	float Lambda;

	void updateDirection(const vec3& w, const float alpha_x, const float alpha_y)
	{
		this->w = w;
		theta = acosf(w.z);
		cosTheta = w.z;
		sinTheta = sinf(theta);
		tanTheta = sinTheta / cosTheta;
		const float invSinTheta2 = 1.0f / (1.0f - w.z * w.z);
		const float cosPhi2 = w.x * w.x * invSinTheta2;
		const float sinPhi2 = w.y * w.y * invSinTheta2;
		alpha = sqrtf(cosPhi2 * alpha_x * alpha_x + sinPhi2 * alpha_y * alpha_y);
		// Lambda
		if (w.z > 0.9999f)
			Lambda = 0.0f;
		else if (w.z < -0.9999f)
			Lambda = -1.0f;
		else
		{
			const float a = 1.0f / tanTheta / alpha;
			Lambda = 0.5f * (-1.0f + ((a > 0) ? 1.0f : -1.0f) * sqrtf(1 + 1 / (a * a)));
		}
	}

	// height
	float h;
	float C1;
	float G1;

	void updateHeight(const float& h)
	{
		this->h = h;
		C1 = std::min(1.0f, std::max(0.0f, 0.5f * (h + 1.0f)));

		if (this->w.z > 0.9999f)
			G1 = 1.0f;
		else if (this->w.z <= 0.0f)
			G1 = 0.0f;
		else
			G1 = powf(this->C1, this->Lambda);
	}
};

inline float invC1(const float U)
{
	const float h = std::max(-1.0f, std::min(1.0f, 2.0f * U - 1.0f));
	return h;
}

inline float sampleHeight(const RayInfo& ray, const float U)
{
	if (ray.w.z > 0.9999f)
		return std::numeric_limits<Float>::max();
	if (ray.w.z < -0.9999f)
	{
		const float value = invC1(U * ray.C1);
		return value;
	}
	if (fabsf(ray.w.z) < 0.0001f)
		return ray.h;

	// probability of intersection
	if (U > 1.0f - ray.G1) // leave the microsurface
		return std::numeric_limits<Float>::max();

	const float h = invC1(
		ray.C1 / powf((1.0f - U), 1.0f / ray.Lambda)
	);
	return h;
}

float D_ggx(const vec3& wm, const float alpha_x, const float alpha_y) {
	if (wm.z <= 0.0f)
		return 0.0f;

	// slope of wm
	const float slope_x = -wm.x / wm.z;
	const float slope_y = -wm.y / wm.z;

	// P22
	const float tmp = 1.0f + slope_x * slope_x / (alpha_x * alpha_x) + slope_y * slope_y / (alpha_y * alpha_y);
	const float P22 = 1.0f / (M_PI * alpha_x * alpha_y) / (tmp * tmp);

	// value
	const float value = P22 / (wm.z * wm.z * wm.z * wm.z);
	return value;
}

vec2 sampleP22_11(const float theta_i, const float U, const float U_2, const float alpha_x, const float alpha_y)
{
	vec2 slope;

	if (theta_i < 0.0001f)
	{
		const float r = sqrtf(U / (1.0f - U));
		const float phi = 6.28318530718f * U_2;
		slope.x = r * cosf(phi);
		slope.y = r * sinf(phi);
		return slope;
	}

	// constant
	const float sin_theta_i = sinf(theta_i);
	const float cos_theta_i = cosf(theta_i);
	const float tan_theta_i = sin_theta_i / cos_theta_i;

	// slope associated to theta_i
	const float slope_i = cos_theta_i / sin_theta_i;

	// projected area
	const float projectedarea = 0.5f * (cos_theta_i + 1.0f);
	if (projectedarea < 0.0001f || projectedarea != projectedarea)
		return vec2(0, 0, 0);
	// normalization coefficient
	const float c = 1.0f / projectedarea;

	const float A = 2.0f * U / cos_theta_i / c - 1.0f;
	const float B = tan_theta_i;
	const float tmp = 1.0f / (A * A - 1.0f);

	const float D = sqrtf(std::max(0.0f, B * B * tmp * tmp - (A * A - B * B) * tmp));
	const float slope_x_1 = B * tmp - D;
	const float slope_x_2 = B * tmp + D;
	slope.x = (A < 0.0f || slope_x_2 > 1.0f / tan_theta_i) ? slope_x_1 : slope_x_2;

	float U2;
	float S;
	if (U_2 > 0.5f)
	{
		S = 1.0f;
		U2 = 2.0f * (U_2 - 0.5f);
	}
	else
	{
		S = -1.0f;
		U2 = 2.0f * (0.5f - U_2);
	}
	const float z = (U2 * (U2 * (U2 * 0.27385f - 0.73369f) + 0.46341f)) / (U2 * (U2 * (U2 * 0.093073f + 0.309420f) - 1.000000f) + 0.597999f);
	slope.y = S * z * sqrtf(1.0f + slope.x * slope.x);

	return slope;
}

vec3 sampleVNDF(const vec3& wi, const float alpha_x, const float alpha_y)
{
	const float U1 = generateRandomNumber();
	const float U2 = generateRandomNumber();

	// sample D_wi

	// stretch to match configuration with alpha=1.0	
	const vec3 wi_11 = normalize(vec3(alpha_x * wi.x, alpha_y * wi.y, wi.z));

	// sample visible slope with alpha=1.0
	vec2 slope_11 = sampleP22_11(acosf(wi_11.z), U1, U2, alpha_x, alpha_y);

	// align with view direction
	const float phi = atan2(wi_11.y, wi_11.x);
	vec2 slope(cosf(phi) * slope_11.x - sinf(phi) * slope_11.y, sinf(phi) * slope_11.x + cos(phi) * slope_11.y, 0);

	// stretch back
	slope.x *= alpha_x;
	slope.y *= alpha_y;

	// if numerical instability
	if ((slope.x != slope.x) || !IsFiniteNumber(slope.x))
	{
		if (wi.z > 0) return vec3(0.0f, 0.0f, 1.0f);
		else return normalize(vec3(wi.x, wi.y, 0.0f));
	}

	// compute normal
	const vec3 wm = normalize(vec3(-slope.x, -slope.y, 1.0f));

	return wm;
}

vec3 samplePhaseFunction_conductor(const vec3& wi, const float alpha_x, const float alpha_y, const Spectrum& m_eta, const Spectrum& m_k, Spectrum& weight)
{
	const float U1 = generateRandomNumber();
	const float U2 = generateRandomNumber();

	// sample D_wi
	// stretch to match configuration with alpha=1.0	
	const vec3 wi_11 = normalize(vec3(alpha_x * wi.x, alpha_y * wi.y, wi.z));

	// sample visible slope with alpha=1.0
	vec2 slope_11 = sampleP22_11(acosf(wi_11.z), U1, U2, alpha_x, alpha_y);

	// align with view direction
	const float phi = atan2(wi_11.y, wi_11.x);
	vec2 slope(cosf(phi) * slope_11.x - sinf(phi) * slope_11.y, sinf(phi) * slope_11.x + cos(phi) * slope_11.y, 0);

	// stretch back
	slope.x *= alpha_x;
	slope.y *= alpha_y;

	// compute normal
	vec3 wm;
	// if numerical instability
	if ((slope.x != slope.x) || !IsFiniteNumber(slope.x))
	{
		if (wi.z > 0) wm = vec3(0.0f, 0.0f, 1.0f);
		else wm = normalize(vec3(wi.x, wi.y, 0.0f));
	}
	else
		wm = normalize(vec3(-slope.x, -slope.y, 1.0f));

	// reflect
	const vec3 wo = -wi + 2.0f * wm * dot(wi, wm);
	weight = fresnelConductorExact(dot(wi, wm), m_eta, m_k);

	return wo;
}

Spectrum evalPhaseFunction_conductor(const RayInfo& ray, const vec3& wo, const float alpha_x, const float alpha_y, const Spectrum& m_eta, const Spectrum& m_k)
{
	if (ray.w.z > 0.9999f)
		return Spectrum(0.0f);

	// half vector 
	const vec3 wh = normalize(-ray.w + wo);
	if (wh.z < 0.0f)
		return Spectrum(0.0f);

	// projected area
	float projectedArea;
	if (ray.w.z < -0.9999f)
		projectedArea = 1.0f;
	else
		projectedArea = ray.Lambda * ray.w.z;

	// value
	const Spectrum value = fresnelConductorExact(dot(-ray.w, wh), m_eta, m_k) * std::max(0.0f, dot(-ray.w, wh)) * D_ggx(wh, alpha_x, alpha_y) / 4.0f / projectedArea / dot(-ray.w, wh);
	return value;
}

// MIS weights for bidirectional path tracing on the microsurface
float MISweight_conductor(const vec3& wi, const vec3& wo, const float alpha_x, const float alpha_y)
{
	if (wi.x == -wo.x && wi.y == -wo.y && wi.z == -wo.z)
		return 1.0f;
	const vec3 wh = normalize(wi + wo);
	const float value = D_ggx((wh.z > 0) ? wh : -wh, alpha_x, alpha_y);
	return value;
}

Spectrum eval_conductor(const vec3& wi, const vec3& wo, const float alpha_x, const float alpha_y, const Spectrum& m_eta, const Spectrum& m_k, const int scatteringOrderMax)
{
	if (wi.z <= 0 || wo.z <= 0)
		return Spectrum(0.0f);

	// init
	RayInfo ray;
	ray.updateDirection(-wi, alpha_x, alpha_y);
	ray.updateHeight(1.0f);
	Spectrum energy(1.0f);

	RayInfo ray_shadowing;
	ray_shadowing.updateDirection(wo, alpha_x, alpha_y);

	// eval single scattering	
	// half-vector
	const vec3 wh = normalize(wi + wo);
	const float D = D_ggx(wh, alpha_x, alpha_y);
	const float G2 = 1.0f / (1.0f + (-ray.Lambda - 1.0f) + ray_shadowing.Lambda);
	Spectrum singleScattering = fresnelConductorExact(dot(-ray.w, wh), m_eta, m_k) * D * G2 / (4.0f * wi.z);

	// MIS weight 
	float wi_MISweight;

	// multiple scattering
	Spectrum multipleScattering(0.0f);

	// random walk
	int current_scatteringOrder = 0;
	while (current_scatteringOrder < scatteringOrderMax)
	{
		// next height
		float U = generateRandomNumber();
		ray.updateHeight(sampleHeight(ray, U));

		// leave the microsurface?
		if (ray.h == std::numeric_limits<Float>::max())
			break;
		else
			current_scatteringOrder++;

		// next event estimation 
		if (current_scatteringOrder > 1) // single scattering is already computed
		{
			Spectrum phasefunction = evalPhaseFunction_conductor(ray, wo, alpha_x, alpha_y, m_eta, m_k);
			ray_shadowing.updateHeight(ray.h);
			float shadowing = ray_shadowing.G1;
			Spectrum I = energy * phasefunction * shadowing;

			// MIS
			const float MIS = wi_MISweight / (wi_MISweight + MISweight_conductor(-ray.w, wo, alpha_x, alpha_y));


			if (IsFiniteNumber(I[0]))
				multipleScattering += I * MIS;
		}

		// next direction
		Spectrum weight;
		ray.updateDirection(samplePhaseFunction_conductor(-ray.w, alpha_x, alpha_y, m_eta, m_k, weight), alpha_x, alpha_y);
		energy = energy * weight;
		ray.updateHeight(ray.h);

		if (current_scatteringOrder == 1)
			wi_MISweight = MISweight_conductor(wi, ray.w, alpha_x, alpha_y);

		// if NaN (should not happen, just in case)
		if ((ray.h != ray.h) || (ray.w.x != ray.w.x))
			return Spectrum(0.0f);
	}

	// 0.5f = MIS weight of singleScattering
	// multipleScattering already weighted by MIS
	return 0.5f * singleScattering + multipleScattering;
}

vec3 sample_conductor(const vec3& wi, const float alpha_x, const float alpha_y, const Spectrum& m_eta, const Spectrum& m_k, const int scatteringOrderMax, Spectrum& energy)
{
	energy = Spectrum(1.0f);

	// init
	RayInfo ray;
	ray.updateDirection(-wi, alpha_x, alpha_y);
	ray.updateHeight(1.0f);

	// random walk
	int current_scatteringOrder = 0;
	while (true)
	{
		// next height
		float U = generateRandomNumber();
		ray.updateHeight(sampleHeight(ray, U));

		// leave the microsurface?
		if (ray.h == std::numeric_limits<Float>::max())
			break;
		else
			current_scatteringOrder++;

		// next direction
		Spectrum weight;
		ray.updateDirection(samplePhaseFunction_conductor(-ray.w, alpha_x, alpha_y, m_eta, m_k, weight), alpha_x, alpha_y);
		energy = energy * weight;
		ray.updateHeight(ray.h);

		// if NaN (should not happen, just in case)
		if ((ray.h != ray.h) || (ray.w.x != ray.w.x))
		{
			energy = Spectrum(0.0f);
			return vec3(0, 0, 1);
		}

		if (current_scatteringOrder > scatteringOrderMax)
		{
			energy = Spectrum(0.0f);
			return vec3(0, 0, 1);
		}
	}

	return ray.w;
}


int main(int argc, char **argv) {
	cout << "Mitsuba Main\n";
	float isotropic_roughness = 0.7f;
	Spectrum m_eta = Spectrum(3.f);
	Spectrum m_k = Spectrum(0.5f);
	
	ofstream io_brdf; // io_brdf is like cin

	io_brdf.open("io_table_test.txt"); // opens the file
	if (!io_brdf) { // file couldn't be opened
		cerr << "Error: file could not be opened" << endl;
		cout << "NO";
		exit(1);
	}

	ofstream brdf_table; // io_brdf is like cin

	brdf_table.open("brdf_test.txt"); // opens the file
	if (!brdf_table) { // file couldn't be opened
		cerr << "Error: file could not be opened" << endl;
		cout << "NO";
		exit(1);
	}

	float ix;
	float iy;
	float iz;

	float ox;
	float oy;
	float oz;
	//varia per ogni componente di incoming e outgoing
	for ( ix = -1.f; ix <= 1.f; ix += 0.25) {
		for ( iy = -1.f; iy <= 1.f; iy += 0.25) {
			for ( iz = 0.f; iz <= 1.f; iz += 0.25) {
				for ( ox = -1.f; ox <= 1.f; ox += 0.25) {
					for ( oy = -1.f; oy <= 1.f; oy += 0.25) {
						for (oz = 0.f; oz <= 1.f; oz += 0.25) {
							
							io_brdf << "Incoming: (";
							io_brdf << ix;
							io_brdf << ",";
							io_brdf << iy;
							io_brdf << ",";
							io_brdf << iz;
							io_brdf << ") - ";
							io_brdf << "Outgoing: (";
							io_brdf << ox;
							io_brdf << ",";
							io_brdf << oy;
							io_brdf << ",";
							io_brdf << oz;
							io_brdf << ") ---> BRDF Value: ";
							
							vec3 incoming = vec3(ix, iy, iz);
							vec3 outgoing = vec3(ox, oy, oz);
							Spectrum brdf = eval_conductor(incoming, outgoing, isotropic_roughness, isotropic_roughness, m_eta, m_k, 5);
							
							io_brdf << brdf.toString();
							io_brdf << "\n";

							brdf_table << brdf[0];
							brdf_table << " ";
							brdf_table << brdf[1];
							brdf_table << " ";
							brdf_table << brdf[2];

							brdf_table << "\n";

						}
						oz = 0.f;
					}
					oy = 0.f;
				}
				oz = 0.f;
			}
			iz = 0.f;
		}
		iy = 0.f;
	}
	io_brdf.close();
			
    /*
	cout << "\n";
    cout << brdf2.toString();
	cout << "\n";  
	cout << "\n";
    cout << brdf6.toString();
	cout << "\n";
	
	
	float eps_xy = 0.002;
	float eps_z = 0.002;
	float count_x = 0;
	float count_y = 0;
	float count_z = 0;
	int count = 0;
	for (int i = 0; i < 10; i++)
		for(int j = 0; j < 10; j++ )
			for (int z = 0; z < 5; z++) {
				cout << "\n";
				cout << count_x;
				cout << ",";
				cout << count_y;
				cout << ",";
				cout << count_z;
				cout << "\n";
				count_x += eps_xy;
				count_y += eps_xy;
				count_z += eps_z;
				count++;
			}
	cout << "Numero di ripetizioni: ";
	cout << count;
	cout << "\n";
	/*
	Spectrum output2 = output + Spectrum(1.f);
    cout << output2.toString()+"\n";
    cout << "hello \n";
    cout << Spectrum((0.0f,0.5f,1.f)).toString()+"\n";
    cout << Spectrum(0.5f).toString()+"\n";
    cout << Spectrum(1.f).toString()+"\n";
	*/
	
}

