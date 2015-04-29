program forcefield
implicit none

!=======================================================================
! initialization 
! 
!=======================================================================


! initialize variables and grids________________________________________
integer::dim_r,dim_theta
integer::dim_r0,dim_theta0
integer::dim_r1,dim_theta1

integer::level

integer::i,j,k,l
real::t_start,t_end,t
real::var0,var1,var2,var3,var4,var5,var6
real::var0s,var1s,var2s,var3s,var4s,var5s,var6s
real::G
real::dr0,dtheta0
real::dr1,dtheta1

real::dforce_r,dforce_theta

real,dimension(128,128)::r_sub1a
real,dimension(128,128)::r_sub1b
real,dimension(128,128)::r_sub1c
real,dimension(128,128)::r_sub1d

real,dimension(256,256)::cos_dtheta1a
real,dimension(256,256)::sin_dtheta1a
real,dimension(256,256)::cos_dtheta1b
real,dimension(256,256)::sin_dtheta1b
real,dimension(256,256)::cos_dtheta1c
real,dimension(256,256)::sin_dtheta1c
real,dimension(256,256)::cos_dtheta1d
real,dimension(256,256)::sin_dtheta1d

real,dimension(128)::r
real,dimension(128)::r0
real,dimension(128)::r1a
real,dimension(128)::r1b
real,dimension(128)::r1c
real,dimension(128)::r1d

real,dimension(256)::theta
real,dimension(256)::theta0
real,dimension(256)::theta1a
real,dimension(256)::theta1b
real,dimension(256)::theta1c
real,dimension(256)::theta1d

real,dimension(128+1,256+1)::density
real,dimension(128,256)::mass_term1a
real,dimension(128,256)::mass_term1b
real,dimension(128,256)::mass_term1c
real,dimension(128,256)::mass_term1d

real,dimension(128,256,2)::force1

G=1
dim_theta0=256
dim_r0=128
dim_theta=dim_theta0
dim_r=dim_r0

level=1

dim_r1=dim_r0/(2**level)
dim_theta1=dim_theta0/(2**level)


!=======================================================================
!=======================================================================


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
do i=1,dim_r0
	density(i,dim_theta0+1)=density(i,1)
end do
do k=1,dim_theta0+1
	density(dim_r0+1,k)=0
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
! build level 1 grids (theta1,r1,cos_dtheta1,sin_dtheta1,r_sub1,mass_term1)
! calculate forcel
!=======================================================================

! build theta1 grid_____________________________________________________
dtheta1=2**level*dtheta0
do j=1,dim_theta1
	theta1a(j)=theta0(1)+(2**level-1)*dtheta0*0.5+dtheta1*(j-1)
	theta1b(j)=theta1a(j)+dtheta0
	theta1c(j)=theta1a(j)
	theta1d(j)=theta1a(j)+dtheta0
end do


! build r1 grid_________________________________________________________
dr1=2**level*dr0
do i=1,dim_r1
	r1a(i)=r0(1)+(2**level-1)*dr0*0.5+dr1*(i-1)
	r1b(i)=r1a(i)+dr0
	r1c(i)=r1a(i)
	r1d(i)=r1c(i)+dr0
end do


! build cos_dtheta1 and sin_dtheta1 grid________________________________
do j=1,dim_theta
	do l=1,dim_theta1
		var0=theta(j)-theta1a(l)
		cos_dtheta1a(j,l)=cos(var0)
		sin_dtheta1a(j,l)=sin(var0)
		var0=theta(j)-theta1b(l)
		cos_dtheta1b(j,l)=cos(var0)
		sin_dtheta1b(j,l)=sin(var0)
		var0=theta(j)-theta1c(l)
		cos_dtheta1c(j,l)=cos(var0)
		sin_dtheta1c(j,l)=sin(var0)
		var0=theta(j)-theta1d(l)
		cos_dtheta1d(j,l)=cos(var0)
		sin_dtheta1d(j,l)=sin(var0)
	end do
end do


! build r_sub1 grid_____________________________________________________
do i=1,dim_r
	do k=1,dim_r1
		r_sub1a(i,k)=r(i)/r1a(k)
		r_sub1b(i,k)=r(i)/r1b(k)
		r_sub1c(i,k)=r(i)/r1c(k)
		r_sub1d(i,k)=r(i)/r1d(k)
	end do
end do


! build mass_term1 grid_________________________________________________
do i=1,dim_r1
	do k=1,dim_theta1
		mass_term1a(i,k)=0
		mass_term1b(i,k)=0
		mass_term1c(i,k)=0
		mass_term1d(i,k)=0
		do j=1,(2**level)
			do l=1,(2**level)
				mass_term1a(i,k)=mass_term1a(i,k)-G*density((i-1)*(2**level)+j,(k-1)*(2**level)+l)*dtheta0*dr0/r1a(i)
				mass_term1b(i,k)=mass_term1b(i,k)-G*density((i-1)*(2**level)+j+1,(k-1)*(2**level)+l)*dtheta0*dr0/r1b(i)
				mass_term1c(i,k)=mass_term1c(i,k)-G*density((i-1)*(2**level)+j,(k-1)*(2**level)+l+1)*dtheta0*dr0/r1c(i)
				mass_term1d(i,k)=mass_term1d(i,k)-G*density((i-1)*(2**level)+j+1,(k-1)*(2**level)+l+1)*dtheta0*dr0/r1d(i)
			end do
		end do
	end do
end do


! calculate time________________________________________________________+(var4-var2s)*var6s
call CPU_Time(t_start)	


! calculate forcel______________________________________________________+var3s*var6s
do j=1,dim_theta
	do i=1,dim_r
		dforce_r=0.0
		dforce_theta=0.0
		! take level grid 1a
		if (Mod(i,2) == 1 .and. Mod(j,2)==1) then
			do l=1,dim_theta1
				var2=cos_dtheta1a(j,l)
				var3=sin_dtheta1a(j,l)
				do k=1,dim_r1
					var4=r_sub1a(i,k)
					var5=1+var4*var4-2*var4*var2
					var5=var5*sqrt(var5)
					var6=mass_term1a(k,l)/(var5)
					dforce_r=dforce_r+(var4-var2)*var6
					dforce_theta=dforce_theta+var3*var6
				end do
			end do
			force1(i,j,1)=dforce_r
			force1(i,j,2)=dforce_theta
		end if

		! take level grid 1b
		if (Mod(i,2) == 0 .and. Mod(j,2) ==1) then
			do l=1,dim_theta1
				var2=cos_dtheta1b(j,l)
				var3=sin_dtheta1b(j,l)
				do k=1,dim_r1
					var4=r_sub1b(i,k)
					var5=1+var4*var4-2*var4*var2
					var5=var5*sqrt(var5)
					var6=mass_term1b(k,l)/(var5)
					dforce_r=dforce_r+(var4-var2)*var6
					dforce_theta=dforce_theta+var3*var6
				end do
			end do
			force1(i,j,1)=dforce_r
			force1(i,j,2)=dforce_theta
		end if

		! take level grid 1c
		if (Mod(i,2) == 1 .and. Mod(j,2) ==0) then
			do l=1,dim_theta1
				var2=cos_dtheta1c(j,l)
				var3=sin_dtheta1c(j,l)
				do k=1,dim_r1
					var4=r_sub1c(i,k)
					var5=1+var4*var4-2*var4*var2
					var5=var5*sqrt(var5)
					var6=mass_term1c(k,l)/(var5)
					dforce_r=dforce_r+(var4-var2)*var6
					dforce_theta=dforce_theta+var3*var6
				end do
			end do
			force1(i,j,1)=dforce_r
			force1(i,j,2)=dforce_theta
		end if

		! take level grid 1d
		if (Mod(i,2) == 0 .and. Mod(j,2) ==0) then
			do l=1,dim_theta1
				var2=cos_dtheta1d(j,l)
				var3=sin_dtheta1d(j,l)
				do k=1,dim_r1
					var4=r_sub1d(i,k)
					var5=1+var4*var4-2*var4*var2
					var5=var5*sqrt(var5)
					var6=mass_term1d(k,l)/(var5)
					dforce_r=dforce_r+(var4-var2)*var6
					dforce_theta=dforce_theta+var3*var6
				end do
			end do
			force1(i,j,1)=dforce_r
			force1(i,j,2)=dforce_theta
		end if

	end do
end do


! calculate time________________________________________________________
call CPU_Time(t_end)				
t=t_end-t_start
write(*,*) "required time for level",level
write(*,*) t,"seconds"


! save force0_project.data______________________________________________
open(unit=4,file="Results/force1_project.data")
do i=1,dim_r
	do j=1,dim_theta
		write(4,'(2e20.10)') force1(i,j,1) , force1(i,j,2)
	end do
end do
close(4)	


!=======================================================================
!=======================================================================


end program
