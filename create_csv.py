import numpy as np
from tqdm import tqdm

brdf3DegreesArray = np.zeros((32768,3),dtype = np.double)

idx = 0
count = 0
print("Extracting Brdf of 3 Freedom Degrees")
with open('Brdf_Table_3_Degrees.txt') as file:
    for line in tqdm(file):
        l = line.split(" ")
        for i in range(len(l)):
            l[i] = l[i].replace("\n","")
            brdf3DegreesArray[count][idx] = l[i]    
            idx+=1
        idx=0
        count+=1

#print("3 Degrees Array:\n",brdf3DegreesArray)

brdf4DegreesArray = np.zeros((1048576,3),dtype = np.double)

idx = 0
count = 0   
print("Extracting Brdf of 4 Freedom Degrees")     
with open('Brdf_Table_4_Degrees.txt') as file:
    for line in tqdm(file):
        l = line.split(" ")
        
        for i in range(len(l)):
            #print(l[i])
            l[i] = l[i].replace("\n","")
            brdf4DegreesArray[count][idx] = l[i]    
            idx+=1
        idx=0
        count+=1

#print("4 Degrees Array:\n",brdf4DegreesArray)

print("Saving into csv file...")

np.savetxt("brdf3.csv", np.asarray(brdf3DegreesArray), fmt='%a', delimiter=',')
np.savetxt("brdf4.csv", np.asarray(brdf4DegreesArray), fmt='%a', delimiter=',')

print("Beautifully done!")