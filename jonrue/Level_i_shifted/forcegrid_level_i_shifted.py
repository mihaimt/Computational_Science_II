# -*- coding: utf-8 -*-
"""
Created on Wed Mar 11 15:40:49 2015

@author: jonas
"""
import numpy as np
import matplotlib.pyplot as plt


a0=[]
b0=[]
a=np.zeros((256,128))
b=np.zeros((256,128))



data=open("level0.txt")
for line in data:
    a0.append(float(line))
data2=open("forcegrid_rcomponent.txt")
for line in data2:
    b0.append(float(line))
    

index=0   
for j in range(256):
    for i in range(128):
        a[j,i]=a0[index]-b0[index]
        b[j,i]=b0[index]
        index=index+1

   
subplot(2,2,1)
imshow(b,origin='lower')
colorbar()

   

subplot(2,2,2)      
imshow(a,origin='lower')
colorbar()

mean_theta=np.zeros(256)
for i in range(256):
    mean_theta[i]=mean(abs(a[i,:]))


mean_r=np.zeros(128)
for i in range(128):
    mean_r[i]=mean(abs(a[:,i]))

subplot(2,2,3)
plot(mean_theta)

subplot(2,2,4)
plot(mean_r)

index=0
mean_odd_odd=0
for i in range(0,256,2):
    for j in range(0,128,2):
        mean_odd_odd+=abs(a[i,j])
        index+=1
mean_odd_odd=mean_odd_odd/index
print "mean_odd_odd:"
print mean_odd_odd


index=0
mean_even_even=0
for i in range(1,256,2):
    for j in range(1,128,2):
        mean_even_even+=abs(a[i,j])
        index+=1
mean_even_even=mean_even_even/index
print "mean_even_even:"
print mean_even_even
        
index=0
mean_even_odd=0
for i in range(1,256,2):
    for j in range(0,128,2):
        mean_even_odd+=abs(a[i,j])
        index+=1
mean_even_odd=mean_even_odd/index
print "mean_even_odd:"
print mean_even_odd       

index=0
mean_odd_even=0
for i in range(0,256,2):
    for j in range(1,128,2):
        mean_odd_even+=abs(a[i,j])
        index+=1
mean_odd_even=mean_odd_even/index
print "mean_odd_even:"
print mean_odd_even   




    
#index=1
#for i in range(50):
#    print a[index,index]
#    index+=2



