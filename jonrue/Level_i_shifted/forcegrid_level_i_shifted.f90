Include 'PreCalculations.f90'
Include 'force_odd_odd.f90'
Include 'CreateGrids.f90'
Include 'force_even_even.f90'
Include 'force_even_odd.f90'
Include 'force_odd_even.f90'

program test

implicit none

integer, parameter :: dim_r=128
integer, parameter :: dim_theta=256
integer :: dim_r_primeLevel1
integer :: dim_theta_primeLevel1
integer :: LevelNr

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


real, dimension(dim_r) :: r
real, dimension(dim_r) :: r_prime
real, dimension(dim_theta) :: theta
real, dimension(dim_theta) :: theta_prime
real, dimension(dim_r,dim_theta) :: density_grid
real, dimension(dim_r,dim_theta) :: force_grid_r
real, dimension(dim_r,dim_theta) :: force_grid_t
real, dimension(dim_r,dim_theta) :: force_grid_rLevel1
real, dimension(dim_r,dim_theta) :: force_grid_tLevel1
real, dimension(dim_r,dim_theta) :: difference_FGR
real, dimension(dim_r,dim_theta) :: difference_FGT
real, dimension(dim_theta,dim_theta) :: sinus_List
real, dimension(dim_theta,dim_theta) :: cosinus_List
real, dimension(dim_r,dim_r) :: R_List
real, dimension(dim_r,dim_r) :: dR_List

real, dimension(:), allocatable :: r_prime_Level1
real, dimension(:), allocatable :: theta_prime_Level1
real, dimension(:,:), allocatable :: sinus_List_Level1
real, dimension(:,:), allocatable :: cosinus_List_Level1
real, dimension(:,:), allocatable :: R_List_Level1
real, dimension(:,:), allocatable :: dR_List_Level1
real, dimension(:,:), allocatable :: density_grid_Level1


LevelNr=1


call CPU_Time(t1)

!Read in Data.........................................................
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
!........................................................

call get_dimension(LevelNr, dim_r, dim_r_primeLevel1)
call get_dimension(LevelNr, dim_theta, dim_theta_primeLevel1)



allocate(r_prime_Level1(dim_r_primeLevel1))
allocate(theta_prime_Level1(dim_theta_primeLevel1))
call create_empty_DensityGrid(dim_r,LevelNr, r_prime, dim_r_primeLevel1, r_prime_Level1)
call create_empty_DensityGrid(dim_theta, LevelNr, theta_prime, dim_theta_primeLevel1, theta_prime_Level1)



delta_r=r_prime(2)-r_prime(1) 
delta_theta=theta_prime(2)-theta_prime(1)
volumina=delta_r*delta_theta




!Shift the grid.....................................................
do i=1,dim_r
	r(i)=r_prime(i)-delta_r*0.5
end do
do i=1,dim_theta
	theta(i)=theta_prime(i)-delta_theta*0.5
end do


call cpu_time(t2)



!Main LOOP......................................
open(unit=3,file="../density/density_project.data")
do i=1,dim_r
	do j=1,dim_theta
		read(3,'(e20.10)') density_grid(i,j)
		density_grid(i,j)=density_grid(i,j)*volumina
	end do
end do
close(3)




allocate(density_grid_Level1(dim_r_primeLevel1,dim_theta_primeLevel1))
call create_real_DensityGrid(density_grid,dim_r,dim_theta,density_grid_Level1,dim_r_primeLevel1,dim_theta_primeLevel1)



do i=2,256,2
	do j=2,128,2
		call force_even_even(j,i,dim_r,dim_theta, r,theta,dim_theta_primeLevel1, &
			dim_r_primeLevel1, density_grid_Level1, force_grid_r, force_grid_t, &
			r_prime_Level1,theta_prime_Level1)
	end do
end do


do i=2,256,2
	do j=1,128,2
		call force_even_odd(j,i,dim_r,dim_theta, r,theta,dim_theta_primeLevel1, &
			dim_r_primeLevel1, density_grid_Level1, force_grid_r, force_grid_t, &
			r_prime_Level1,theta_prime_Level1)
	end do
end do

do i=1,256,2
	do j=2,128,2
		call force_odd_even(j,i,dim_r,dim_theta, r,theta,dim_theta_primeLevel1, &
			dim_r_primeLevel1, density_grid_Level1, force_grid_r, force_grid_t, &
			r_prime_Level1,theta_prime_Level1)
	end do
end do

do i=1,256,2
	do j=1,128,2
		call force_odd_odd(j,i,dim_r,dim_theta, r,theta,dim_theta_primeLevel1, &
			dim_r_primeLevel1, density_grid_Level1, force_grid_r, force_grid_t, &
			r_prime_Level1,theta_prime_Level1)
	end do
end do


!.................................................................................
			
call cpu_time(t3)


write(*,*) t3-t1
write(*,*) t3-t2

open(unit=4,file="forcegrid_rcomponent.txt")
open(unit=5,file="forcegrid_thetacomponent.txt")

do i=1,256
	do j=1,128
		write(4,'(e20.10)') force_grid_r(j,i)
		write(5,'(e20.10)') force_grid_t(j,i)
	
	end do
end do
close(4)
close(5)


end program
