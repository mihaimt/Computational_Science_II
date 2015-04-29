program forcefield
implicit none

!=======================================================================
! initialization 
! 
!=======================================================================


! initialize variables and grids________________________________________
integer::dim_r,dim_theta
integer::dim_r0,dim_theta0

integer::i,j,k,l
real::t_start,t_end,t
real::var0,var1,var2,var3,var4,var5,var6
real::var0s,var1s,var2s,var3s,var4s,var5s,var6s
real::G
real::dr0,dtheta0

real::dforce_r,dforce_theta

real,dimension(128,128)::r_sub0
real,dimension(256,256)::cos_dtheta0
real,dimension(256,256)::sin_dtheta0 

real,dimension(128)::r
real,dimension(128)::r0

real,dimension(256)::theta
real,dimension(256)::theta0

real,dimension(128,256)::density
real,dimension(128,256)::mass_term0

real,dimension(128,256,2)::force0

integer::level
integer::step
integer::start_r
G=1
dim_theta0=256
dim_r0=128
dim_theta=dim_theta0
dim_r=dim_r0


!=======================================================================
!=======================================================================



!=======================================================================
! build force grid (r,theta)
! loading data level 0 grids (density,r0,theta0)
!=======================================================================
		
	
! load data from Computational_Science_II/Data/_________________________
open(unit=1,file="../../Data/r_project.data")
do i=1,dim_r0
	read(1,'(e20.10)') r0(i)
end do
close(1)
open(unit=2,file="../../Data/theta_project.data")
do i=1,dim_theta0
	read(2,'(e20.10)') theta0(i)
end do
close(2)
open(unit=3,file="../../Data/density_project.data")
do i=1,dim_r0
	do k=1,dim_theta0
		read(3,'(e20.10)') density(i,k)
	end do
end do
close(3)


! build theta grid______________________________________________________
dtheta0=theta0(2)-theta0(1)
do j=1,dim_theta0
	theta(j)=theta0(j)-dtheta0*0.5
end do


! build r grid__________________________________________________________
dr0=r0(2)-r0(1)
do i=1,dim_r0
	r(i)=r0(i)-dr0*0.5
end do


!=======================================================================
!=======================================================================


!=======================================================================
! build level 0 grids (cos_dtheta0,sin_dtheta0,r_sub0,mass_term0)
! calculate force0
!=======================================================================


! build cos_dtheta0 and sin_dtheta0 grid________________________________
do j=1,dim_theta
	do l=1,dim_theta0
		var0=theta(j)-theta0(l)
		cos_dtheta0(j,l)=cos(var0)
		sin_dtheta0(j,l)=sin(var0)
	end do
end do


! build r_sub0 grid_____________________________________________________
do i=1,dim_r
	do k=1,dim_r0
		r_sub0(i,k)=r(i)/r0(k)
	end do
end do


! build mass_term0 grid_________________________________________________
do k=1,dim_r0
	do l=1,dim_theta0
		mass_term0(k,l)=-G*density(k,l)*dtheta0*dr0/r0(k)
	end do
end do


! calculate time________________________________________________________
call CPU_Time(t_start)	


! calculate force0______________________________________________________
do j=1,dim_theta
	do i=1,dim_r
		dforce_r=0.0
		dforce_theta=0.0
		do l=1,dim_theta0
			var2=cos_dtheta0(j,l)
			var3=sin_dtheta0(j,l)
			do k=1,dim_r0
				var4=r_sub0(i,k)
				var5=1+var4*var4-2*var4*var2
				var5=var5*sqrt(var5)
				var6=mass_term0(k,l)/(var5+0.0000001)
				dforce_r=dforce_r+(var4-var2)*var6
				dforce_theta=dforce_theta+var3*var6
			end do
		end do
		force0(i,j,1)=dforce_r
		force0(i,j,2)=dforce_theta
	end do
end do


! calculate time________________________________________________________
call CPU_Time(t_end)				
t=t_end-t_start
write(*,*) "required time for level 0:"
write(*,*) t


! save force0_project.data______________________________________________
open(unit=4,file="Res/force0_project.data")
do i=1,dim_r
	do j=1,dim_theta
		write(4,'(2e20.10)') force0(i,j,1) , force0(i,j,2)
	end do
end do
close(4)	



!=======================================================================
!=======================================================================



end program

