program disc_force
implicit none

integer, parameter :: dimR = 128, dimT = 256
real, parameter :: eps = 1e-6

real, dimension(dimR) :: rGrid
real, dimension(dimT) :: tGrid
real, dimension(dimR,dimT) :: dns0
real, dimension(dimR,dimT) :: dns1
real, dimension(dimR,dimT) :: dnsTemp
real, dimension(dimR,dimR) :: rList
real, dimension(dimR,dimR) :: drList
real, dimension(dimT,dimT) :: cosList
real, dimension(dimT,dimT) :: sinList

real, dimension(dimR,dimT) :: forceR
real, dimension(dimR,dimT) :: forceT

real :: temp1, temp2, cosTemp, sinTemp, rTemp, denom, scal
real :: deltaR,deltaT,volumina

integer :: r,t,rP,tP

integer :: jump, tStart, rStart




open (unit = 1, file = "density/r_project.data") 
do rP=1,dimR
	read (1, *) rGrid(rP)
end do
close(1)

open (unit = 2, file = "density/theta_project.data")
do tP=1,dimT
	read (2,*) tGrid(tP)
end do
close(2)

deltaR = rGrid(2)-rGrid(1)
deltaT = tGrid(2)-tGrid(1)
volumina = deltaT*deltaR

do r=1,dimR
    do rP=1,dimR
        rList(r,rP) = (r-0.5)/(rP-0.5)
        drList(r,rP) = -(r-0.5)/((rP-0.5)*(rP-0.5)*deltaR)
    end do
end do

do t = 1,dimT
    do tP = 1,dimT
        temp1 = deltaT*(t-tP-1.0)            
        cosList(t,tP) = cos(temp1)
        sinList(t,tP) = sin(temp1)
    end do
end do
        

!!!!!!!!!!!!!!!!!! MAIN LOOP  !!!!!!!!!!!!!!!!!!!!!


open (unit = 3, file = "density/density_project.data")
do rP=1,dimR
  	do tP=1,dimT
    	read(3,*) dns0(rP,tP)
        dns0(rP,tP) = dns0(rP,tP) * volumina
    end do
end do
close(3)


do rP=1,dimR-1
  	do tP=1,dimT-1
        dnsTemp(rP,tP) = dns0(rP,tP)+dns0(rP,tP+1)
    end do
end do

write(*,*) "Here"

do rP=1,dimR-1
  	do tP=1,dimT-1
        dns1(rP,tP) = dnsTemp(rP,tP)+dnsTemp(rP+1,tP)
    end do
end do

do rP=1,dimR
    dnsTemp(rP,dimT) = dns0(rP,dimT)+dns0(rP,1)
end do

do rP=1,dimR
    dns1(rP,dimT) = dnsTemp(rP+1,dimT)+dnsTemp(rP,dimT)
end do

do tP=1,dimT-1
    dns1(dimR,tP)=dns0(dimR,tP)+dns0(dimR,tP+1)
end do

jump = 2 !2**(lvl-1)
write(*,*) "Here"


do t = 1,dimT
!write(*,*) t
    do r = 1,dimR
        temp1  = 0
        temp2  = 0
        tStart = 2-mod(t,jump)
	rStart  = 2-mod(r,jump)
        !write(*,*) t, "   ", jump, "    ", tStart
     
     
        do tP = tStart,dimT,jump
            cosTemp = cosList(t,tP)
            sinTemp = sinList(t,tP)
            
            
            do rP = rStart,dimR,jump
                rTemp = rList(r,rP)
                denom = (rTemp*(rTemp-2*cosTemp)+1)
                denom = denom*sqrt(denom)*rTemp+eps
                scal  = dns1(rP,tP)*drList(r,rP)/denom
                temp1 = temp1 + scal*(rTemp-cosTemp)
                temp2 = temp2 + scal*sinTemp
                
            end do
        end do
        forceR(r,t) = temp1
        forceT(r,t) = temp2
        
    end do
end do
  
  
!!!!!!!!!!Write!!!!!!

open (unit=4, file="output/force_t.txt")
open (unit=5, file="output/force_r.txt")
do t=1,dimT
    do r=1,dimR
        write(4,*)forceT(r,t)
        write(5,*)forceR(r,t)
    end do
end do
close(4)
close(5)





write(*,*) "done Bravo"
end program disc_force

