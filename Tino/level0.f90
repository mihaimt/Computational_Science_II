program disc_force
implicit none

real, dimension(128) :: lenR                 !Distance between R Gridpoints         - From File
real, dimension(129) :: lenRShift            !Shifted Points of above
real, dimension(256) :: lenT                 !Distance between Theta Gridpoints     - From File
real, dimension(128,256) :: dens             !Density Grid                          - From File
real, dimension(129,257) :: forceR           !Force in R Direction                  - To be calculated and writen to file    
real, dimension(129,257) :: forceT           !Force in Theta Direction              - To be calculated and writen to file 
real, dimension(257,256) :: cosL             !List of possible CosinusValues
real, dimension(257,256) :: sinL             !List of possible SinusValues
real, dimension(129,128) :: rL               !New dimensionless integrationValue r/rP 
real, dimension(129,128) :: drL              !Integrant of above (-r/rP**2)
real :: difR,difT,area                       !Differece between grid Points
real :: scal,denom                           !Scalar in front of ForceVektor(consisting of denom (=denominator) and values defined below)
real :: rrT,cosT,sinT                        !Temporary Values - used to not always to go back to the big Arrays defined above
real :: temp1,temp2                          !Temporary Values - used to split the Lists for Cos and Sin - and to split the dimension of ForceVektor
real :: time1, time2                         !To Measure Time
integer::r,t,rP,tP                           !LoopIndex 
integer::nR,nT                               !Length of R and Theta


call CPU_Time(time1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Reading Files
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

nR = 128
nT = 256

open (unit = 1, file = "r_project.data") 
do r=1,nR
	read (1, *) lenR(r)
end do
close(1)


open (unit = 2, file = "theta_project.data")
do t=1,nT
	read (2,*) lenT(t)
end do
close(2)



open (unit = 3, file = "radl.data")
do r=1,nR
  	do t=1,nT
    	read(3,*) dens(r,t)
    end do
end do
close(3)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Calculating
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


difR = (lenR(2)-lenR(1))                        !Difference between Gridpoints and Area
difT = (lenT(2)-lenT(1))
area = difR*difT

do t = 1,nT+1
    do tP = 1,nT
        temp1 = difT*(+t-tP-0.5)                !Listing possible sinus and cosinus Values
        cosL(t,tP) = cos(temp1)
        sinL(t,tP) = sin(temp1)
    end do
end do

difR = difR/2
difT = difT/2

do r = 1,128
    lenRShift(r) = lenR(r)-difR                 !Shifting the R Gridpoints to the edge 
end do
lenRShift(nR+1) = lenR(nR)+difR 


do rP=1,nR
    do r=1,nR+1
        rL(r,rP)  = (r-0.5)/rP                  !Calculating the possible Values for the Integrant
        drL(r,rP) = (r-0.5)/(rP*rP)
    end do
end do



do t=1,nT                                             
  do r=1,nR+1
    temp1 = 0
    temp2 = 0
    do tP=1,nT
        cosT = cosL(t,tP)                                   !Buffer Cosinus and Sinus Values                
        sinT = sinL(t,tP)
        do rP=1,nR
            rrT      = rL(r,rP)                             !Buffer Integrant Value
            denom    = (rrT*(rrT-2*cosT)+1)                 !Calculating Denominator
            denom    = rrT*denom*sqrt(denom)
            scal     = dens(rP,tP)*drL(r,rP)/denom          !Calculating Scalar
            temp1    = temp1 + scal*(rrT-cosT)              !Splitting RValues
            temp2    = temp2 + scal*sinT                    !Splitting ThetaValued
      end do
    end do
    forceR(r,t) = area*temp1                                !Multiplying Area 
    forceT(r,t) = area*temp2                                !!! HERE SHOULD ALSO BE THE GRAVITATION VALUE !!!
  end do
end do

call CPU_Time(time2)
write (*,*) "We are done with calculating my dear friend - duration:", time2-time1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Writing Files
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

open (unit=4, file="force_t.txt")
open (unit=5, file="force_r.txt")
do t=1,nT
    do r=1,nR+1
        write(4,*)forceR(r,t)
        write(5,*)forceT(r,t)
    end do
end do
close(4)
close(5)

      

end program disc_force
