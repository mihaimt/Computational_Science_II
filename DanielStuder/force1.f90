program forcefield
implicit none

!=======================================================================
! initialization 
! 
!=======================================================================


! initialize variables and grids________________________________________
integer::dim_r,dim_theta
integer::dim_r0,dim_theta0
integer::dim_rl,dim_thetal

integer::level

integer::i,j,k,l
real::t_start,t_end,t
real::var0,var1,var2,var3,var4,var5,var6
real::G
real::dr0,dtheta0
real::drl,dthetal

real::dforce_r,dforce_theta

real,dimension(128,128)::r_subl
real,dimension(256,256)::cos_dthetal
real,dimension(256,256)::sin_dthetal

real,dimension(128)::r
real,dimension(128)::r0
real,dimension(128)::rl

real,dimension(256)::theta
real,dimension(256)::theta0
real,dimension(256)::thetal

real,dimension(128,256)::density
real,dimension(128,256)::mass_terml

real,dimension(128,256,2)::forcel

G=1
dim_theta0=256
dim_r0=128
dim_theta=dim_theta0
dim_r=dim_r0

level=2
dim_rl=dim_r0/(2**level)
dim_thetal=dim_theta0/(2**level)


!=======================================================================
!=======================================================================


!=======================================================================
! build force grid (r,theta)
! loading data level 0 grids (density,r0,theta0)
!=======================================================================
		
	
! load data from Computational_Science_II/Data/_________________________
open(unit=1,file="Data/r_project.data")
do i=1,dim_r0
	read(1,'(e20.10)') r0(i)
end do
close(1)
open(unit=2,file="Data/theta_project.data")
do i=1,dim_theta0
	read(2,'(e20.10)') theta0(i)
end do
close(2)
open(unit=3,file="Data/density_project.data")
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
! build level l grids (thetal,rl,cos_dthetal,sin_dthetal,r_subl,mass_terml)
! calculate forcel
!=======================================================================

! build thetal grid_____________________________________________________
dthetal=2**level*dtheta0
do j=1,dim_thetal
	thetal(j)=theta0(1)+(2**level-1)*dtheta0*0.5+dthetal*(j-1)
end do


! build rl grid_________________________________________________________
drl=2**level*dr0
do i=1,dim_rl
	rl(i)=r0(1)+(2**level-1)*dr0*0.5+drl*(i-1)
end do


! build cos_dthetal and sin_dthetal grid________________________________
do j=1,dim_theta
	do l=1,dim_thetal
		var0=theta(j)-thetal(l)
		cos_dthetal(j,l)=cos(var0)
		sin_dthetal(j,l)=sin(var0)
	end do
end do


! build r_subl grid_____________________________________________________
do i=1,dim_r
	do k=1,dim_rl
		r_subl(i,k)=r(i)/rl(k)
	end do
end do


! build mass_terml grid_________________________________________________
do i=1,dim_rl
	do k=1,dim_thetal
		mass_terml(i,k)=0
		do j=1,(2**level)
			do l=1,(2**level)
				mass_terml(i,k)=mass_terml(i,k)+G*density((i-1)*(2**level)+j,(k-1)*(2**level)+l)*dtheta0*dr0/rl(i)
			end do
		end do
	end do
end do


! calculate time________________________________________________________
call CPU_Time(t_start)	


! calculate forcel______________________________________________________
do j=1,dim_theta
	do i=1,dim_r
		dforce_r=0.0
		dforce_theta=0.0
		do l=1,dim_thetal
			var2=cos_dthetal(j,l)
			var3=sin_dthetal(j,l)
			do k=1,dim_rl
				var4=r_subl(i,k)
				var5=1+var4*var4-2*var4*var2
				var5=var5*sqrt(var5)
				var6=mass_terml(k,l)/(var5+0.0000001)
				dforce_r=dforce_r+(var4-var2)*var6
				dforce_theta=dforce_theta+var3*var6
			end do
		end do
		forcel(i,j,1)=dforce_r
		forcel(i,j,2)=dforce_theta
	end do
end do


! calculate time________________________________________________________
call CPU_Time(t_end)				
t=t_end-t_start
write(*,*) "required time for level l:"
write(*,*) t


! save force0_project.data______________________________________________
open(unit=4,file="Results/forcel_project.data")
do i=1,dim_r
	do j=1,dim_theta
		write(4,'(2e20.10)') forcel(i,j,1) , forcel(i,j,2)
	end do
end do
close(4)	


!=======================================================================
!=======================================================================


end program
