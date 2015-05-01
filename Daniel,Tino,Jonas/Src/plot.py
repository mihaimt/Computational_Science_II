# -*- coding: utf-8 -*-
"""
Spyder Editor

This temporary script file is located here:
/home/omer/.spyder2/.temp.py
"""

from pylab import*
hold(False)
"""
Make plot of density
"""
density=loadtxt("../Data/density_project.data")
r=loadtxt("../Data/r_project.data")
theta=loadtxt("../Data/theta_project.data")
dim_r=len(r)
dim_theta=len(theta)
dens=zeros((dim_r,dim_theta))

for j in arange(0,dim_theta):
    for i in arange(0,dim_r):
        dens[i][j]=density[i*256+j]
        
#pcolormesh(theta,r,dens)
#title('$Density$')
#ylabel('Radius $r$ [$r_0$]')
#xlabel('Azimuth $\Theta$ [$rad$]')
#ylim([r[0],r[dim_r-1]])
#xlim([theta[0],theta[dim_theta-1]])
#colorbar()
#savefig("Res/Density.png")
        
"""
Make plots of force level 0 
"""
force0=loadtxt("Res/force0_project.data")

dr=r[1]-r[0]
r=r-dr/2
dtheta=theta[1]-theta[0]
theta=theta-dtheta/2 

ftheta0=zeros((dim_r,dim_theta))
fr0=zeros((dim_r,dim_theta))

for j in arange(0,dim_theta):
    for i in arange(0,dim_r):
        fr0[i][j]=force0[i*256+j,0]
        ftheta0[i][j]=force0[i*256+j,1]
        
F0=sqrt(ftheta0**2+fr0**2)

#pcolormesh(theta,r,F0)
#title('$Force$ $Level$ $0$')
#ylabel('$Radius$ $r$ [$r_0$]')
#xlabel('$Azimuth$ $\Theta$ [$rad$]')
#ylim([r[0],r[dim_r-1]])
#xlim([theta[0],theta[dim_theta-1]])
#colorbar()
#savefig("Res/Force_0.png")

#pcolormesh(theta,r,fr0)
#title('$Force_r$ $Level$ $0$')
#ylabel('$Radius$ $r$ [$r_0$]')
#xlabel('$Azimuth$ $\Theta$ [$rad$]')
#ylim([r[0],r[dim_r-1]])
#xlim([theta[0],theta[dim_theta-1]])
#colorbar()
#savefig("Res/Force_r_0.png")

#pcolormesh(theta,r,ftheta0)
#title('$Force_theta$ $Level$ $0$')
#ylabel('$Radius$ $r$ [$r_0$]')
#xlabel('$Azimuth$ $\Theta$ [$rad$]')
#ylim([r[0],r[dim_r-1]])
#xlim([theta[0],theta[dim_theta-1]])
#colorbar()
#savefig("Res/Force_theta_0.png")

"""
Make plots of force level l 
"""
forcel=loadtxt("Res/force_project.data")

fthetal=zeros((dim_r,dim_theta))
frl=zeros((dim_r,dim_theta))

for j in arange(0,dim_theta):
    for i in arange(0,dim_r):
        frl[i][j]=forcel[i*256+j,0]
        fthetal[i][j]=forcel[i*256+j,1]
        
Fl=sqrt(fthetal**2+frl**2)

#pcolormesh(theta,r,Fl)
#title('$Force$ $Level$ $l$')
#ylabel('$Radius$ $r$ [$r_0$]')
#xlabel('$Azimuth$ $\Theta$ [$rad$]')
#ylim([r[0],r[dim_r-1]])
#xlim([theta[0],theta[dim_theta-1]])
#colorbar()
#savefig("Res/Force_l.png")

#pcolormesh(theta,r,frl)
#title('$Force_r$ $Level$ $l$')
#ylabel('$Radius$ $r$ [$r_0$]')
#xlabel('$Azimuth$ $\Theta$ [$rad$]')
#ylim([r[0],r[dim_r-1]])
#xlim([theta[0],theta[dim_theta-1]])
#colorbar()
#savefig("Res/Force_r_l.png")

#pcolormesh(theta,r,fthetal)
#title('$Force_theta$ $Level$ $l$')
#ylabel('$Radius$ $r$ [$r_0$]')
#xlabel('$Azimuth$ $\Theta$ [$rad$]')
#ylim([r[0],r[dim_r-1]])
#xlim([theta[0],theta[dim_theta-1]])
#colorbar()
#savefig("Res/Force_theta_l.png")

"""
Make plots of force level l-0 differences 
"""
dfthetal=fthetal-ftheta0
dfrl=frl-fr0
dFl=Fl-F0


#pcolormesh(theta,r,dFl)
#title('$\Delta F$ $Level$ $l$')
#ylabel('$Radius$ $r$ [$r_0$]')
#xlabel('$Azimuth$ $\Theta$ [$rad$]')
#ylim([r[0],r[dim_r-1]])
#xlim([theta[0],theta[dim_theta-1]])
#colorbar()
#savefig("Res/dForce_l.png")

#pcolormesh(theta,r,dfrl)
#title('$\Delta F_r$ $Level$ $l$')
#ylabel('$Radius$ $r$ [$r_0$]')
#xlabel('$Azimuth$ $\Theta$ [$rad$]')
#ylim([r[0],r[dim_r-1]])
#xlim([theta[0],theta[dim_theta-1]])
#colorbar()
#savefig("Res/dForce_r_l.png")

pcolormesh(theta,r,fthetal)
title('$\Delta F_\Theta$ $Level$ $l$')
ylabel('$Radius$ $r$ [$r_0$]')
xlabel('$Azimuth$ $\Theta$ [$rad$]')
ylim([r[0],r[dim_r-1]])
xlim([theta[0],theta[dim_theta-1]])
colorbar()
savefig("Res/dForce_theta_l.png")



"""
Make frames of force level l-0 differences 
"""
## make dF(theta) for each r
#for i in range(128):
#    hold(False)
#    plot(theta,dFl[i,:])
#    xlim([theta[0],theta[255]])
#    ylim([-1.3,1.3])
#    title('$Level 1$ $Oscillations$')
#    xlabel('$\Theta [rad]$')
#    ylabel('$\Delta F(\Theta)$')
#    name=("$r =$%f" %r[i])
#    legend([name])
#    if i<10:
#        savefig("../Res/Anim/Level1/DF_theta/dfr00%i.png" %i)
#    elif i<100:
#        savefig("../Res/Anim/Level1/DF_theta/dfr0%i.png" %i)
#    else:
#        savefig("../Res/Anim/Level1/DF_theta/dfr%i.png" %i)


## make dF(r) for each theta
#for i in range(256):
#    hold(False)
#    plot(r,dFl[:,i])
#    xlim([r[0],r[127]])
#    ylim([-1.3,1.3])
#    title('$Level 1$ $Oscillations$')
#    xlabel('$r [r_0]$')
#    ylabel('$\Delta F(r)$')
#    name=("$\Theta =$%f" %theta[i])
#    legend([name])
#    if i<10:
#        savefig("../Res/Anim/Level1/DF_r/dfr00%i.png" %i)
#    elif i<100:
#        savefig("../Res/Anim/Level1/DF_r/dfr0%i.png" %i)
#    else:
#        savefig("../Res/Anim/Level1/DF_r/dfr%i.png" %i)


## make dFr(theta) for each r
#for i in range(128):
#    hold(False)
#    plot(theta,dfrl[i,:])
#    xlim([theta[0],theta[255]])
#    ylim([-1.3,1.3])
#    title('$Level 1$ $Oscillations$')
#    xlabel('$\Theta [rad]$')
#    ylabel('$\Delta F_{r}(\Theta)$')
#    name=("$r =$%f" %r[i])
#    legend([name])
#    if i<10:
#        savefig("../Res/Anim/Level1/DFr_theta/dfr00%i.png" %i)
#    elif i<100:
#        savefig("../Res/Anim/Level1/DFr_theta/dfr0%i.png" %i)
#    else:
#        savefig("../Res/Anim/Level1/DFr_theta/dfr%i.png" %i)

  
## make dFtheta(theta) for each r
#for i in range(128):
#    hold(False)
#    plot(theta,dfthetal[:,i])
#    xlim([theta[0],theta[255]])
#    ylim([-0.1,0.1])
#    title('$Level 1$ $Oscillations$')
#    xlabel('$\Theta [rad]$')
#    ylabel('$\Delta F_{\Theta}(\Theta)$')
#    name=("$r =$%f" %r[i])
#    legend([name])
#    if i<10:
#        savefig("../Res/Anim/Level1/DFtheta_theta/dftheta00%i.png" %i)
#    elif i<100:
#        savefig("../Res/Anim/Level1/DFtheta_theta/dftheta0%i.png" %i)
#    else:
#        savefig("../Res/Anim/Level1/DFtheta_theta/dftheta%i.png" %i)


## make dFr(r) for each theta
#for i in range(256):
#    hold(False)
#    plot(r,dfrl[i,:])
#    xlim([r[0],r[127]])
#    ylim([-1.3,1.3])
#    title('$Level 1$ $Oscillations$')
#    xlabel('$r [r_0]$')
#    ylabel('$\Delta F_{r}(r)$')
#    name=("$\Theta =$%f" %theta[i])
#    legend([name])
#    if i<10:
#        savefig("Animation/Level1/DFr_r/dfr00%i.png" %i)
#    elif i<100:
#        savefig("Animation/Level1/DFr_r/dfr0%i.png" %i)
#    else:
#        savefig("Animation/Level1/DFr_r/dfr%i.png" %i)

## make dFtheta(r) for each theta
#for i in range(256):
#    hold(False)
#    plot(r,dfthetal[i,:])
#    xlim([r[0],r[127]])
#    ylim([-0.1,0.1])
#    title('$Level 1$ $Oscillations$')
#    xlabel('$r [r_0]$')
#    ylabel('$\Delta F_{\Theta}(r)$')
#    name=("$\Theta =$%f" %theta[i])
#    legend([name])
#    if i<10:
#        savefig("Animation/Level1/DFtheta_r/dftheta00%i.png" %i)
#    elif i<100:
#        savefig("Animation/Level1/DFtheta_r/dftheta0%i.png" %i)
#    else:
#        savefig("Animation/Level1/DFtheta_r/dftheta%i.png" %i)






