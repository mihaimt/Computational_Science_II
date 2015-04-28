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



data=open("forcegrid_rcomponent_level0.txt")
for line in data:
    a0.append(float(line))
data2=open("forcegrid_rcomponent.txt")
for line in data2:
    b0.append(float(line))
    

index=0   
for j in range(256):
    for i in range(128):
        a[j,i]=a0[index]-b0[index]
        index=index+1
subplot(1,2,1)      
imshow(a,origin='lower')
colorbar()

subplot(1,2,2)
index2=50
hold(True)
for i in range(1):
    plot(a[index2,0:50])
    index2+=2
    
#index=1
#for i in range(50):
#    print a[index,index]
#    index+=2



