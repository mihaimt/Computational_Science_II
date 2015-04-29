program forcefield
implicit none

! allocate memory
integer::dim_r=128
integer::dim_theta=256

integer::level
integer::levels=8

integer::levelsize
integer::halfsize
integer::start_theta
integer::start_r

integer::i,j,k,l,s,m
real::t_start,t_end,t
real::var0,var1,var2,var3,var4,var5,var6

real::G
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

real, dimension(:,:,:), allocatable::mass_term

real, dimension(:,:,:), allocatable::force

allocate(r(dim_r))
allocate(r0(dim_r))

allocate(theta(dim_theta))
allocate(theta0(dim_theta))

allocate(density(dim_r,dim_theta))

allocate(r_sub(dim_r,dim_theta))
allocate(cos_dtheta(dim_theta,dim_theta))
allocate(sin_dtheta(dim_theta,dim_theta))

allocate(mass_term(levels,dim_r,dim_theta))

allocate(force(2,dim_r,dim_theta))


! load data 
open(unit=1,file="../Data/r_project.data")
do k=1,dim_r
	read(1,'(e20.10)') r0(k)
end do
close(1)
open(unit=2,file="../Data/theta_project.data")
do l=1,dim_theta
	read(2,'(e20.10)') theta0(l)
end do
close(2)
open(unit=3,file="../Data/density_project.data")
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

! build mass_term grid
call CPU_Time(t_start)	
do level=1,1
	halfsize=2**(level-1)
	do k=1,dim_r
		do l=1,dim_theta
			mass_term(level,k,l)=0
		end do
	end do			 

	do k=1,dim_r
		do l=1,dim_theta
			var1=-G*density(k,l)*dtheta*dr
			do i=k-halfsize+1,k+halfsize
				do j=l-halfsize+1,l+halfsize
					if(i>0.and.i<=dim_r) then
						if(j>0.and.j<=dim_theta) then
							mass_term(level,i,j)=mass_term(level,i,j)+var1!/r(i)
						end if
					end if
				end do
			end do
		end do
	end do
end do
call CPU_Time(t_end)				
t=t_end-t_start
write(*,*) "required time generating mass term level 1-8"
write(*,*) t,"seconds"



! calculate force
call CPU_Time(t_start)	
level=1
levelsize=2**level
do j=1,dim_theta
	start_theta=1+mod(mod(j,levelsize)+2**(level-1)-1,levelsize)
	do i=1,dim_r
		start_r=1+mod(mod(i,levelsize)+2**(level-1)-1,levelsize)
		dforce_r=0.0
		dforce_theta=0.0
		do l=2,dim_theta,2!levelsize
			var2=cos_dtheta(j,l)
			var3=sin_dtheta(j,l)
			do k=2,dim_r,2!levelsize
				var4=r_sub(i,k)
				var5=1+var4*var4-2*var4*var2
				var5=var5*sqrt(var5)
				var6=mass_term(level,k,l)/(var5+0.000001)
				dforce_r=dforce_r+(var4-var2)*var6
				dforce_theta=dforce_theta+var3*var6
			end do
		end do
		force(1,i,j)=dforce_r
		force(2,i,j)=dforce_theta
	end do
end do
call CPU_Time(t_end)				
t=t_end-t_start
write(*,*) "required time for level",level
write(*,*) t,"seconds"

! save force_project.data
open(unit=4,file="Results/forcel_project.data")
do i=1,dim_r
	do j=1,dim_theta
		write(4,'(2e20.10)') force(1,i,j) , force(2,i,j)
	end do
end do
close(4)

! save mass_project.data
open(unit=4,file="Results/massl_project.data")
do i=1,dim_r
	do j=1,dim_theta
		write(4,'(2e20.10)') mass_term(1,i,j)
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

deallocate(mass_term)

deallocate(force)

end program

subroutine make_r_sub(a,b,c)

real, dimension(:) :: a,b,c

end subroutine

