program forcelevel
implicit none
!cd Documents/Uni/FS15/Computational_Science_II/DanielStuder/04.25.15/
!gfortran -o force Src/force.f90
! allocate memory
integer, parameter::dim_r=128
integer, parameter::dim_theta=256

integer::level
integer::levels=8

integer::l_size
integer::l_hsize


integer::i,j,k,l,u,v
real::t_start,t_end,t
real::var0,var1,var2,var3,var4
real::a,b,c,d,e

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


! build mass grid levels 
level=0
write(*,*) "calculate force level:"
READ(*,"(1I1)") level

call CPU_Time(t_start)	

l_size=2**level
l_hsize=l_size/2

do i=1,dim_r
	do j=1,dim_theta
		mass(level,i,j)=0
	end do
end do
	
do k=1+l_hsize,dim_r,l_size
	do l=1+l_hsize,dim_theta,l_size

        do i=k-l_hsize,k+l_hsize-1
            do j=l-l_hsize,l+l_hsize-1
                mass(level,k,l)=mass(level,k,l)-G*density(i,j)*dr*dtheta
            end do
        end do
    
	end do
end do

call CPU_Time(t_end)				
t=t_end-t_start
write(*,*) "required time generating mass term leves"
write(*,*) t,"seconds"


! calculate force
call CPU_Time(t_start)	
do j=1,dim_theta
	v=mod(j-1,l_size)
	do i=1,dim_theta
		dforce_r=0.0
		dforce_theta=0.0
		u=mod(i-1,l_size)

		a=real((l_size-v)*(l_size-u))
		b=real(u*(l_size-v))
		c=real(v*(l_size-u))
		d=real(v*u)
		e=1.0/l_size**2
	
		var4=0.0
		do l=1+l_hsize,dim_theta,l_size
			do k=1+l_hsize,dim_r,l_size
				if(k+l_size<=dim_r.and.l+l_size<=dim_theta) then
					var0=cos_dtheta(j,l+v)
					var1=sin_dtheta(j,l+v)
					if(a>0) then
						var4=var4+a*mass(level,k,l)
					end if
					if(b>0) then
						var4=var4+b*mass(level,k+l_size,l)
					end if
					if(c>0) then
						var4=var4+c*mass(level,k,l+l_size)
					end if
					if(d>0) then
						var4=var4+d*mass(level,k+l_size,l+l_size)
					end if		
					var4=var4*e
					var2=r_sub(i,k+u)
					var3=1+var2*var2-2*var2*var0
					var3=var3*sqrt(var3)
					var3=var4/var3/r(k+u)
					dforce_r=dforce_r+(var2-var0)*var3
					dforce_theta=dforce_theta+var1*var3
				end if
				if(k+l_size<=dim_r.and.l+l_size>dim_theta) then

					if(l+v<=dim_theta) then
						var0=cos_dtheta(j,l+v)
						var1=sin_dtheta(j,l+v)
					else
						var0=cos_dtheta(j,l+v-dim_theta)
						var1=sin_dtheta(j,l+v-dim_theta)
					end if

					if(a>0) then
						var4=var4+a*mass(level,k,l)
					end if
					if(b>0) then
						var4=var4+b*mass(level,k+l_size,l)
					end if
					if(c>0) then
						var4=var4+c*mass(level,k,l+l_size-dim_theta)
					end if
					if(d>0) then
						var4=var4+d*mass(level,k+l_size,l+l_size-dim_theta)
					end if		
					var4=var4*e
					var2=r_sub(i,k+u)
					var3=1+var2*var2-2*var2*var0
					var3=var3*sqrt(var3)
					var3=var4/var3/r(k+u)
					dforce_r=dforce_r+(var2-var0)*var3
					dforce_theta=dforce_theta+var1*var3	
				end if		
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


! save force_project.data______________________________________________
open(unit=4,file="Res/force_project.data")
do i=1,dim_r
	do j=1,dim_theta
		write(4,'(2e20.10)') force(1,i,j) , force(2,i,j)
	end do
end do
close(4)

call system("python Src/plot.py")

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
