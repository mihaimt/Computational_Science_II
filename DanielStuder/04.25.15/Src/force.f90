program forcelevel
implicit none
!cd Documents/Uni/FS15/Computational_Science_II/DanielStuder/04.25.15/
!gfortran -o force Src/force.f90
! allocate memory
integer::dim_r=128
integer::dim_theta=256

integer::level
integer::levels=8

integer::level_size
integer::level_halfsize
integer::pos_theta
integer::pos_r

integer::i,j,k,l
real::t_start,t_end,t
real::var0,var1,var2,var3,var4

real::G=1
real::dr,dtheta
real::dforce_r,dforce_theta

real, dimension(:), allocatable::r
real, dimension(:), allocatable::r0
real, dimension(:), allocatable::theta
real, dimension(:), allocatable::theta0

real, dimension(:,:), allocatable::density

real, dimension(:,:), allocatable::r_sub
real, dimension(:,:), allocatable::cos_dtheta
real, dimension(:,:), allocatable::sin_dtheta

real, dimension(:,:,:), allocatable::mass

real, dimension(:,:,:), allocatable::force


allocate(r(dim_r))
allocate(r0(dim_r))

allocate(theta(dim_theta))
allocate(theta0(dim_theta))

allocate(density(dim_r,dim_theta))

allocate(r_sub(dim_r,dim_theta))
allocate(cos_dtheta(dim_theta,dim_theta))
allocate(sin_dtheta(dim_theta,dim_theta))

allocate(mass(levels,dim_r,dim_theta))

allocate(force(2,dim_r,dim_theta))


! load data 
open(unit=1,file="../../Data/r_project.data")
do k=1,dim_r
	read(1,'(e20.10)') r0(k)
end do
close(1)
open(unit=2,file="../../Data/theta_project.data")
do l=1,dim_theta
	read(2,'(e20.10)') theta0(l)
end do
close(2)
open(unit=3,file="../../Data/density_project.data")
do k=1,dim_r
	do l=1,dim_theta
		read(3,'(e20.10)') density(k,l)
	end do
end do
close(3)


! build theta coordinate
dtheta=theta0(2)-theta0(1)
do j=1,dim_theta
	theta(j)=theta0(j)-dtheta*0.5
end do


! build r coordinate
dr=r0(2)-r0(1)
do i=1,dim_r
	r(i)=r0(i)-dr*0.5
end do


! build cos_dtheta and sin_dtheta grid
do j=1,dim_theta
	do l=1,dim_theta
		var0=theta(j)-theta(l)
		cos_dtheta(j,l)=cos(var0)
		sin_dtheta(j,l)=sin(var0)
	end do
end do


! build r_subl grid
do i=1,dim_r
	do k=1,dim_r
		r_sub(i,k)=r(i)/r(k)
	end do
end do


! build mass grid levels 
level=0
write(*,*) "calculate force level:"
READ(*,"(1I1)") level

call CPU_Time(t_start)	

level_size=2**level
level_halfsize=level_size/2
do i=1,dim_r
	do j=1,dim_theta
		mass(level,i,j)=0
	end do
end do	
do k=1+level_halfsize,dim_r,level_size
	do l=1+level_halfsize,dim_theta,level_size
	!write(*,*) k,l

	do i=k-level_halfsize,k+level_halfsize-1
		do j=l-level_halfsize,l+level_halfsize-1
			mass(level,k,l)=mass(level,k,l)-G*density(i,j)*dr*dtheta
		end do
	end do

	!write(*,*) mass(level,k,l)
	end do
end do
call CPU_Time(t_end)				
t=t_end-t_start
write(*,*) "required time generating mass term leves"
write(*,*) t,"seconds"


!level=2
!level_size=2**level
!do i=1,dim_r
!	write(*,*) 1+mod(mod(i,level_size)+2**(level-1)-1,level_size)
!end do
!level=1
!level_size=2**level



! calculate force
call CPU_Time(t_start)	
do j=1,dim_theta
	pos_theta=mod(mod(j,level_size)+level_size-1,level_size)
	write(*,*) pos_theta
	do i=1,dim_r
		dforce_r=0.0
		dforce_theta=0.0
		pos_r=mod(mod(i,level_size)+level_size-1,level_size)
		if(pos_r==0.and.pos_theta==0) then
			do l=1+level_halfsize,dim_theta,level_size
				do k=1+level_halfsize,dim_r,level_size
					var0=cos_dtheta(j,l)
					var1=sin_dtheta(j,l)
					var2=r_sub(i,k)
					var3=1+var2*var2-2*var2*var0
					var3=var3*sqrt(var3)
					var3=mass(level,k,l)/(var3)/r(k)
					dforce_r=dforce_r+(var2-var0)*var3
					dforce_theta=dforce_theta+var1*var3
				end do
			end do
		end if
		if(pos_r==1.and.pos_theta==0.and.k+level_size<=dim_r) then
			do l=1+level_halfsize,dim_theta,level_size
				do k=1+level_halfsize,dim_r,level_size
					var0=cos_dtheta(j,l)
					var1=sin_dtheta(j,l)
					var2=r_sub(i,k+1)
					var3=1+var2*var2-2*var2*var0
					var3=var3*sqrt(var3)
					var3=(mass(level,k,l)+mass(level,k+level_size,l))/2/(var3)/r(k+1)
					dforce_r=dforce_r+(var2-var0)*var3
					dforce_theta=dforce_theta+var1*var3
				end do
			end do
		end if
		if(pos_r==0.and.pos_theta==1.and.l+level_size<=dim_theta) then
			do l=1+level_halfsize,dim_theta,level_size
				do k=1+level_halfsize,dim_r,level_size
					var0=cos_dtheta(j,l+1)
					var1=sin_dtheta(j,l+1)
					var2=r_sub(i,k)
					var3=1+var2*var2-2*var2*var0
					var3=var3*sqrt(var3)
					var3=(mass(level,k,l)+mass(level,k,l+level_size))/2/(var3)/r(k)
					dforce_r=dforce_r+(var2-var0)*var3
					dforce_theta=dforce_theta+var1*var3
				end do
			end do
		end if
		if(pos_r==1.and.pos_theta==1.and.k+level_size<=dim_r.and.l+level_size<=dim_theta) then
			do l=1+level_halfsize,dim_theta,level_size
				do k=1+level_halfsize,dim_r,level_size
					var0=cos_dtheta(j,l+1)
					var1=sin_dtheta(j,l+1)
					var2=r_sub(i,k+1)
					var3=1+var2*var2-2*var2*var0
					var3=var3*sqrt(var3)
					var4=mass(level,k,l)
					var4=var4+mass(level,k+level_size,l)
					var4=var4+mass(level,k+level_size,l+level_size)
					var4=var4+mass(level,k,l+level_size)
					var4=var4/4/(var3)/r(k+1)
					dforce_r=dforce_r+(var2-var0)*var4
					dforce_theta=dforce_theta+var1*var4
				end do
			end do
		end if
					
		force(1,i,j)=dforce_r
		force(2,i,j)=dforce_theta
		!write(*,*) force(1,i,j),force(2,i,j)
	end do 
end do
call CPU_Time(t_end)				
t=t_end-t_start
write(*,*) "required time for level",level
write(*,*) t,"seconds"


! save force_project.data______________________________________________
open(unit=4,file="Res/force_project.data")
do i=1,dim_r
	do j=1,dim_theta
		write(4,'(2e20.10)') force(1,i,j) , force(2,i,j)
	end do
end do
close(4)


! deallocate memory
deallocate(r)
deallocate(r0)

deallocate(theta)
deallocate(theta0)

deallocate(density)

deallocate(r_sub)
deallocate(cos_dtheta)
deallocate(sin_dtheta)

deallocate(mass)


deallocate(force)

end program
