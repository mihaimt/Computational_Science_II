program forcefield
implicit none

!=======================================================================
! initialization 
! 
!=======================================================================

! initialize variables and grids________________________________________
integer::dim_r,dim_theta
integer::dim_r0=128
integer::dim_theta0=256
integer::dim_rl,dim_thetal

integer::level

integer::i,j
integer::k,l
real::t_start,t_end,t
real::var0,var1,var2,var3,var4,var5,var6

real::G
real::dr0,dtheta0
real::drl,dthetal

real::dforce_r,dforce_theta
integer::levelsize
integer::step
integer::startl
integer::start_r
integer::start_theta
real,dimension(128,128)::r_sub0
real,dimension(128,128)::r_subl
real,dimension(256,256)::cos_dtheta0
real,dimension(256,256)::sin_dtheta0
real,dimension(256,256)::cos_dthetal
real,dimension(256,256)::sin_dthetal

real,dimension(128)::r
real,dimension(128)::r0

real,dimension(256)::theta
real,dimension(256)::theta0


real,dimension(128,256)::density
real,dimension(128,256)::mass_term0
real,dimension(8,128,256)::mass_terml

real,dimension(128,256,2)::forcel
integer::shift
G=1
!dim_theta0=256
!dim_r0=128
dim_theta=dim_theta0
dim_r=dim_r0




dim_rl=dim_r0/(2**level)
dim_thetal=dim_theta0/(2**level)


!=======================================================================
!=======================================================================

!level = 1
!write(*,*) "calculate force level:"
!READ(*,"(1I1)") level
!=======================================================================
! build force grid (r,theta)
! loading data level 0 grids (density,r0,theta0)
!=======================================================================
		
	
! load data from Computational_Science_II/Data/_________________________
open(unit=1,file="../Data/r_project.data")
do i=1,dim_r0
	read(1,'(e20.10)') r0(i)
end do
close(1)
open(unit=2,file="../Data/theta_project.data")
do i=1,dim_theta0
	read(2,'(e20.10)') theta0(i)
end do
close(2)
open(unit=3,file="../Data/density_project.data")
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
!
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

!=======================================================================
!=======================================================================


!=======================================================================
! build level l grids (cos_dthetal,sin_dthetal,r_subl,mass_terml)
!
!=======================================================================

! build cos_dthetal and sin_dthetal grid________________________________
do j=1,dim_theta
	do l=1,dim_theta
		var0=theta(j)-theta(l)
		cos_dthetal(j,l)=cos(var0)
		sin_dthetal(j,l)=sin(var0)
	end do
end do

! build r_subl grid_____________________________________________________
do i=1,dim_r
	do k=1,dim_r
		r_subl(i,k)=r(i)/r(k)
	end do
end do

! calculate time________________________________________________________
call CPU_Time(t_start)

! build mass_terml grid_________________________________________________
do level=1,8
	do j=1,dim_r
		do l=1,dim_theta
			mass_terml(level,j,l)=0
		end do
	end do

	levelsize= 2*level/2
	do j=1,dim_r
		do l=1,dim_theta
	
			mass_terml(level,j,l)=mass_terml(level,j,l)-G*density(j,l)*dtheta0*dr0/r(j)

			do step=1,levelsize

				if(l+step<=dim_theta) then
					mass_terml(level,j,l+step)=mass_terml(level,j,l+step)-G*density(j,l)*dtheta0*dr0/r(j)
				end if

				if(j+step<=dim_r) then
					mass_terml(level,j+step,l)=mass_terml(level,j+step,l)-G*density(j,l)*dtheta0*dr0/r(j+step)
				end if

				if(j+step<=dim_r.and.l+step<=dim_theta) then
					mass_terml(level,j+step,l+step)=mass_terml(level,j+step,l+step)-G*density(j,l)*dtheta0*dr0/r(j+step)
				end if

			end do
	
		end do
	end do
end do



! calculate time________________________________________________________
call CPU_Time(t_end)				
t=t_end-t_start
write(*,*) "required time generating mass term level 1-8"
write(*,*) t,"seconds"
call CPU_Time(t_start)	


! calculate forcel______________________________________________________

level=1
step=2**level
do j=1,dim_theta
	start_theta=1+mod(mod(j,step)+2**(level-1)-1,step)
	write(*,*) start_theta
	do i=1,dim_r
		start_r=1+mod(mod(i,step)+2**(level-1)-1,step)
		write(*,*) start_r
		dforce_r=0.0
		dforce_theta=0.0
		! level0
		do l=1,dim_theta0
			var2=cos_dtheta0(j,l)
			var0=sin_dtheta0(j,l)
			do k=1,dim_r0
				var0=r_sub0(j,l)
			end do
		end do
	end do
end do
do j=1,dim_theta
	!start_theta=1+mod(mod(j,step)+2**(level-1)-1,step)
	do i=1,dim_r
		!start_r=1+mod(mod(i,step)+2**(level-1)-1,step)
		dforce_r=0.0
		dforce_theta=0.0

		! level0
		!do l=1,dim_theta0
		!	var2=cos_dtheta0(j,l)
		!	var3=sin_dtheta0(j,l)
		!	do k=1,dim_r0
		!		var4=r_sub0(i,k)
		!		var5=1+var4*var4-2*var4*var2
		!		var5=var5*sqrt(var5)
		!		var6=mass_term0(k,l)/(var5+0.0000001)
		!		dforce_r=dforce_r+(var4-var2)*var6
		!		dforce_theta=dforce_theta+var3*var6
		!	end do
		!end do

		! level lev with oscillations
		do l=3,dim_theta,2
			var2=cos_dthetal(j,l)
			var3=sin_dthetal(j,l)
			do k=3,dim_r,2
				var4=r_subl(i,k)
				var5=1+var4*var4-2*var4*var2
				var5=var5*sqrt(var5)
				var6=mass_terml(1,k,l)/(var5+0.0000001)
				dforce_r=dforce_r+(var4-var2)*var6
				dforce_theta=dforce_theta+var3*var6
			end do
		end do
	
		! level lev without oscillations
		!do l=start_theta,dim_theta,step
		!	var2=cos_dthetal(j,l)
		!	var3=sin_dthetal(j,l)
		!	do k=start_r,dim_r,step
		!		var4=r_subl(i,k)
		!		var5=1+var4*var4-2*var4*var2
		!		var5=var5*sqrt(var5)
		!		var6=mass_terml(lev,k,l)/(var5)
		!		dforce_r=dforce_r+(var4-var2)*var6
		!		dforce_theta=dforce_theta+var3*var6
		!	end do
		!end do

	

		forcel(i,j,1)=dforce_r
		forcel(i,j,2)=dforce_theta
	end do
end do


! calculate time________________________________________________________
call CPU_Time(t_end)				
t=t_end-t_start
write(*,*) "required time for level",level
write(*,*) t,"seconds"


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


