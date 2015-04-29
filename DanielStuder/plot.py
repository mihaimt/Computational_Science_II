# -*- coding: utf-8 -*-
"""
Spyder Editor

This temporary script file is located here:
/home/omer/.spyder2/.temp.py
"""

from pylab import*

force0=loadtxt("Results/force0_project.data")
forcel=loadtxt("Results/forcel_project.data")
forcelosc=loadtxt("Results/forcelosc_project.data")
density=loadtxt("../Data/density_project.data")
f_r0=force0[:,0]
f_theta0=force0[:,1]

f_rl=forcel[:,0]
f_thetal=forcel[:,1]

f_rl=forcelosc[:,0]
f_thetal=forcelosc[:,1]
r=loadtxt("../Data/r_project.data")
theta=loadtxt("../Data/theta_project.data") 

dr=r[1]-r[0]
r=r-dr/2 

dtheta=theta[1]-theta[0]
theta=theta-dtheta/2 

fr0=zeros((256,128))
frl=zeros((256,128))
frlosc=zeros((256,128))
ftheta0=zeros((256,128))
fthetal=zeros((256,128))
fthetalosc=zeros((256,128))
dfrl=zeros((256,128))
dfthetal=zeros((256,128))
dfrlosc=zeros((256,128))
dfthetalosc=zeros((256,128))
dens=zeros((256,128))
F0=zeros((256,128))
Fl=zeros((256,128))
dFl=zeros((256,128))
for i in arange(0,128):
    for j in arange(0,256):
        ftheta0[j][i]=f_theta0[i*256+j]
        fr0[j][i]=f_r0[i*256+j]
        F0[j][i]=sqrt(ftheta0[j][i]**2+fr0[j][i]**2)
        dens[j][i]=density[i*256+j]

for i in arange(0,128):
    for j in arange(0,256):
        fthetal[j][i]=f_thetal[i*256+j]
        frl[j][i]=f_rl[i*256+j]
        Fl[j][i]=sqrt(fthetal[j][i]**2+frl[j][i]**2)


        

dfrl=frl-fr0
dfthetal=fthetal-ftheta0

dFl=Fl-F0

#plot(r,dfl[256/2,:])
#plot(r,dfl[0,:])
#plot(r,dfl[255,:])
#plot(r,dfl[255/4,:])
#plot(r,dfl[3*255/4,:])
#xlim([r[0],r[127]])
#title('Level 1 Oscillations')
#xlabel('$r$')
#ylabel('$\Delta F_{\Theta}(r)$')
#legend(["$\Theta = \pi$","$\Theta = 0$","$\Theta = 2\pi$","$\Theta = \pi/2$","$\Theta = 3\pi/2$"])

## make dF(theta) for each r
#for i in range(128):
#    hold(False)
#    plot(theta,dFl[:,i])
#    xlim([theta[0],theta[255]])
#    ylim([-1.3,1.3])
#    title('$Level 1$ $Oscillations$')
#    xlabel('$\Theta [rad]$')
#    ylabel('$\Delta F(\Theta)$')
#    name=("$r =$%f" %r[i])
#    legend([name])
#    if i<10:
#        savefig("Animation/Level1/DF_theta/dfr00%i.png" %i)
#    elif i<100:
#        savefig("Animation/Level1/DF_theta/dfr0%i.png" %i)
#    else:
#        savefig("Animation/Level1/DF_theta/dfr%i.png" %i)
#
## make dF(r) for each theta
#for i in range(256):
#    hold(False)
#    plot(r,dFl[i,:])
#    xlim([r[0],r[127]])
#    ylim([-1.3,1.3])
#    title('$Level 1$ $Oscillations$')
#    xlabel('$r [r_0]$')
#    ylabel('$\Delta F(r)$')
#    name=("$\Theta =$%f" %theta[i])
#    legend([name])
#    if i<10:
#        savefig("Animation/Level1/DF_r/dfr00%i.png" %i)
#    elif i<100:
#        savefig("Animation/Level1/DF_r/dfr0%i.png" %i)
#    else:
#        savefig("Animation/Level1/DF_r/dfr%i.png" %i)
#
## make dFr(theta) for each r
#for i in range(128):
#    hold(False)
#    plot(theta,dfrl[:,i])
#    xlim([theta[0],theta[255]])
#    ylim([-1.3,1.3])
#    title('$Level 1$ $Oscillations$')
#    xlabel('$\Theta [rad]$')
#    ylabel('$\Delta F_{r}(\Theta)$')
#    name=("$r =$%f" %r[i])
#    legend([name])
#    if i<10:
#        savefig("Animation/Level1/DFr_theta/dfr00%i.png" %i)
#    elif i<100:
#        savefig("Animation/Level1/DFr_theta/dfr0%i.png" %i)
#    else:
#        savefig("Animation/Level1/DFr_theta/dfr%i.png" %i)
#        
#    
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
#        savefig("Animation/Level1/DFtheta_theta/dftheta00%i.png" %i)
#    elif i<100:
#        savefig("Animation/Level1/DFtheta_theta/dftheta0%i.png" %i)
#    else:
#        savefig("Animation/Level1/DFtheta_theta/dftheta%i.png" %i)
#
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

#fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
#ax.contourf(theta, r, dens.T)
#ylim([0,2])


pcolormesh(r,theta,frl)
title('$\Delta F$ $Level0$')
xlabel('Radius $r$ [$r_0$]')
ylabel('Azimuth $\Theta$ [$rad$]')
xlim([r[0],r[127]])
ylim([theta[0],theta[255]])
colorbar()

