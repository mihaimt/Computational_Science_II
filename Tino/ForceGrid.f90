program disc_force
implicit none

integer, parameter :: nR=128                 !Length of R                           - to be known 
integer, parameter :: nT=256                 !Length of T                           - to be known

real, dimension(nR) :: lenR                 !Distance between R Gridpoints         - From File
real, dimension(nT) :: lenT                 !Distance between Theta Gridpoints     - From File
real, dimension(nR,nT) :: dens             !Density Grid                          - From File
real, dimension(nR+1,nT) :: forceR           !Force in R Direction                  - To be calculated and writen to file    
real, dimension(nR+1,nT) :: forceT           !Force in Theta Direction              - To be calculated and writen to file 
real, dimension(nT,nT) :: cosL               !List of possible CosinusValues
real, dimension(nT,nT) :: sinL               !List of possible SinusValues
real, dimension(nR+1,nR) :: rL               !New dimensionless integrationValue r/rP 
real, dimension(nR+1,nR) :: drL              !Integrant of above (-r/rP**2)
real, dimension(:), allocatable :: lenRP     !Distance between RP Gridpoints        - since size level depending - dynamic memory allocation
real, dimension(:), allocatable :: lenTP     !Distance between TP Gridpoints        - since size level depending - dynamic memory allocation
real, dimension(:,:), allocatable::densP     !Prime density

real :: difR,difT,area                       !Differece between grid Points
real :: scal,denom                           !Scalar in front of ForceVektor(consisting of denom (=denominator) and values defined below)
real :: rrT,cosT,sinT                        !Temporary Values - used to not always to go back to the big Arrays defined above
real :: temp1,temp2                          !Temporary Values - used to split the Lists for Cos and Sin - and to split the dimension of ForceVektor
real :: time1, time2                         !To Measure Time
real :: eps                                  !Avoiding Singularity
integer::r,t,rP,tP                           !LoopIndex 
integer::nRP,nTP                             !Length of RPrime and ThetaPrime       
integer:: lvl, jump, start, ind              !needed Variables


call CPU_Time(time1)
lvl = 1                                      !Choosing level


jump = 2**lvl                                !Distance between two Points (Prime)                               
nRP = nR/jump                               
nTP = nT/jump
eps = 0.000001
start = 1+2**(lvl-1)                         !Where to start - PrimeView


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Reading Files
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    

allocate(lenRP(nRP))
allocate(lenTP(nTP))
allocate(densP(nRP,nTP))



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


open (unit = 3, file = "density_project.data")
do r=1,nR
  	do t=1,nT
    	read(3,*) dens(r,t)
    end do
end do
close(3)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Calculating
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


if (lvl==0) then
    
    difR = (lenR(2)-lenR(1))                        !Difference between Gridpoints and Area
    difT = (lenT(2)-lenT(1))
    area = difR*difT
    
    do t = 1,nT
        do tP = 1,nTP
            temp1 = difT*(+t-tP-0.5)                !Listing possible sinus and cosinus Values
            cosL(t,tP) = cos(temp1)
            sinL(t,tP) = sin(temp1)
        end do
    end do

    difR = difR/2
    difT = difT/2

    do rP=1,nRP
        do r=1,nR+1
            rL(r,rP)  = (r-0.5)/rP                  !Calculating the possible Values for the Integrant
            drL(r,rP) = (r-0.5)/(rP*rP)
        end do
    end do



else !if Level!=0 

 
    
    
    ind = 0                                         !Index Value
    do rP=start,nR-1,jump                           !Listing Possible RP Values
        ind = ind + 1
        lenRP(ind)=lenR(rP)
    end do

    
    ind = 0                                         
    do tP=start,nT-1,jump
        ind = ind + 1                               !Listing Possible TP Values
        lenTP(ind)=lenT(tP)
    end do
    
    difR = (lenRP(2)-lenR(1))
    difT = (lenTP(2)-lenR(1))                       !Calculating Area
    area = difR*difT

    do t=1,nT
        do tP = 1,nTP
            temp1 = lenT(t)-lenTP(tP)
            cosL(t,tP) = cos(temp1)                 !Listing possible Cos/Sin Values
            sinL(t,tP) = sin(temp1)
        end do
    end do


    do rP=1,nRP
        do r=1,nR
            temp1     = lenR(r)/lenRP(rP)
            rL(r,rP)  = temp1                    !Calculating the possible Values for the Integrant
            drL(r,rP) = temp1/lenRP(rP)
        end do
    end do
 
  
    
    do rP=1,nRP
        do tP=1,nTP
            temp1 = 0
            do r=0,jump-1
                do t=0,jump-1
                    temp1 = temp1+ dens(jump*rP+r,jump*tP+t)!Summing density into one point
                end do
            end do
            densP(rP,tP) = temp1
        end do
    end do

   
end if ! level0 or level else

do t=1,nT                                             
  do r=1,nR+1
    temp1 = 0
    temp2 = 0
    do tP=1,nTP
        cosT = cosL(t,tP)                                   !Buffer Cosinus and Sinus Values                
        sinT = sinL(t,tP)
        do rP=1,nRP
            rrT      = rL(r,rP)                             !Buffer Integrant Value
            denom    = (rrT*(rrT-2*cosT)+1)                 !Calculating Denominator
            denom    = rrT*denom*sqrt(denom)+eps
            scal     = densP(rP,tP)*drL(r,rP)/denom          !Calculating Scalar
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

!deallocate(densP)
!deallocate(lenRP(nRP))
!deallocate(lenTP(nTP))

open (unit=4, file="force_t1.txt")
open (unit=5, file="force_r1.txt")
do t=1,nT
    do r=1,nR+1
        write(4,*)forceR(r,t)
        write(5,*)forceT(r,t)
    end do
end do
close(4)
close(5)

end program disc_force

