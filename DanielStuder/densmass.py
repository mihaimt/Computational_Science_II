# -*- coding: utf-8 -*-
"""
Created on Tue Apr 21 11:40:46 2015

@author: omer
"""
from pylab import*
density=loadtxt("../Data/density_project.data")
mass=loadtxt("Results/massl_project.data")
r=loadtxt("../Data/r_project.data")
theta=loadtxt("../Data/theta_project.data") 
dr=r[1]-r[0]
r=r-dr/2 
dtheta=theta[1]-theta[0]
theta=theta-dtheta/2 
d=zeros((256,128))
m=zeros((256,128))

for i in arange(0,128):
    for j in arange(0,256):
            d[j][i]=density[i*256+j]*dr*dtheta
            m[j][i]=mass[i*256+j]

print(mean(d-m))
pcolormesh(r,theta,d-m)
title('$\Delta F$ $Level0$')
xlabel('Radius $r$ [$r_0$]')
ylabel('Azimuth $\Theta$ [$rad$]')
xlim([r[0],r[127]])
ylim([theta[0],theta[255]])
colorbar()