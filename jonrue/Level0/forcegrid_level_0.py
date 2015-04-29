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
c=np.zeros((256,128))
d=np.zeros((256,128))


data=open("../output/forcegrid_rcomponent.txt")
for line in data:
    a0.append(float(line))
data2=open("../output/forcegrid_thetacomponent.txt")
for line in data2:
    b0.append(float(line))
    

index=0   
for j in range(256):
    for i in range(128):
        a[j][i]=a0[index]
        index=index+1
imshow(a,origin='lower')
colorbar()
#
#index=0   
#for j in range(257):
#    for i in range(129):
#        c[j][i]=b0[index]
#        index=index+1
#imshow(c,origin='lower')
#colorbar()

#index=0   
#for j in range(257):
#    for i in range(129):
#        d[j][i]=sqrt(a0[index]*a0[index]+b0[index]*b0[index])
#        index=index+1
#imshow(d,origin='lower')
#colorbar()


#    
#data3=open("density_project.data")
#c0=[]
#for line in data3:
#    c0.append(float(line))
#
#b=np.zeros((256,128))
#index=0   
#for j in range(128):
#    for i in range(256):
#        b[i][j]=c0[index]
#        index=index+1
#subplot(1,2,2)
#imshow(b,origin='lower')
#colorbar()
    
