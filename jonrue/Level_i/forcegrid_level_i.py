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



data=open("../output/forcegrid_difference_r.txt")
for line in data:
    a0.append(float(line))
data2=open("../output/forcegrid_rcomponent.txt")
for line in data2:
    b0.append(float(line))
    

index=0   
for j in range(256):
    for i in range(128):
        a[j][i]=a0[index]
        index=index+1
subplot(1,2,1)      
imshow(a,origin='lower')
colorbar()

subplot(1,2,2)
index2=0
hold(True)
for i in range(5):
    plot(a[index2][:])
    index2+=1
    



