program test
implicit none


integer :: dim_r
integer :: dim_theta

integer :: i
integer :: j
integer :: index_theta
integer :: index_thetaPrime
integer :: index_r
integer :: index_rPrime

real :: t1
real :: t2
real :: t3
real :: delta_r
real :: delta_theta
real :: volumina
real :: r_component
real :: theta_component
real :: R_temp, cos_temp, sin_temp
real :: numerator
real :: denominator
real :: fraction_value

real, dimension(128) :: r
real, dimension(128) :: r_prime
real, dimension(256) :: theta
real, dimension(256) :: theta_prime
real, dimension(128,256) :: density_grid
real, dimension(128,256) :: force_grid_r
real, dimension(128,256) :: force_grid_t
real, dimension(256,256) :: sinus_List
real, dimension(256,256) :: cosinus_List
real, dimension(128,128) :: R_List
real, dimension(128,128) :: dR_List


dim_r=128
dim_theta=256


call CPU_Time(t1)

!Begin: Read in Data.........................................................
open(unit=1,file="../density/r_project.data")
do i=1,dim_r
	read (1,'(e20.10)') r_prime(i)
end do
close(1)
open(unit=2,file="../density/theta_project.data")
do i=1,dim_theta
	read(2,'(e20.10)') theta_prime(i)
end do
close(2)
open(unit=3,file="../density/density_project.data")
do i=1,dim_r
	do j=1,dim_theta
	read(3,'(e20.10)') density_grid(i,j)
	end do
end do
close(3)
!End: Read in Data........................................................

write(*,*) theta_prime(1)

delta_r=r_prime(2)-r_prime(1) 
delta_theta=theta_prime(2)-theta_prime(1)
volumina=delta_r*delta_theta

!Begin:Shit the grid.....................................................
do i=1,dim_r
	r(i)=r_prime(i)-delta_r*0.5
end do
do i=1,dim_theta
	theta(i)=theta_prime(i)-delta_theta*0.5
end do
!End: Shift the grid.....................................................

!Begin: Pre-Calculation of Values needed inside the main Loop.........................
do i=1,dim_theta
	do j=1,dim_theta
	cosinus_List(i,j)=cos((-0.5+i-j)*delta_theta)
	sinus_List(i,j)=sin((-0.5+i-j)*delta_theta)
	end do
end do
do i=1,dim_r
	do j=1,dim_r
		R_List(i,j)=(i-0.5)/j
		dR_List(i,j)=-(i-0.5)/(j*j*delta_r)
	end do
end do
!End: Pre-Calculation of Values needed inside the main Loop............................

call cpu_time(t2)


!Begin: Main LOOP......................................
do index_theta=1,dim_theta
	do index_r=1,dim_r
		r_component=0
		theta_component=0
		do index_thetaPrime=1,dim_theta
			cos_temp=cosinus_List(index_theta,index_thetaPrime)
			sin_temp=sinus_List(index_theta,index_thetaPrime)
			do index_rPrime=1,dim_r
				R_temp=R_List(index_r,index_rPrime)
				numerator=density_grid(index_rPrime,index_thetaPrime)*dR_List(index_r,index_rPrime)
				denominator=(1+R_temp*(R_temp-2*cos_temp))
				denominator=denominator*sqrt(denominator)*R_temp
				fraction_value=numerator/denominator
				r_component=(R_temp-cos_temp)*fraction_value+r_component
				theta_component=sin_temp*fraction_value+theta_component

			end do
		end do
		force_grid_r(index_r,index_theta)=volumina*r_component
		force_grid_t(index_r,index_theta)=volumina*theta_component

	end do
end do
!End: Main LOOP...........................................
			
call cpu_time(t3)

write(*,*) t3-t1
write(*,*) t3-t2

open(unit=4,file="../output/forcegrid.txt")
open(unit=5,file="../output/forcegrid2.txt")
do i=1,256
	do j=1,128
		write(4,'(e20.10)') force_grid_r(j,i)
		write(5,'(e20.10)') force_grid_t(j,i)
	end do
end do
close(4)
close(5)


end program
