program forcefield
implicit none

! initialize variables and grids
integer::dim_r,dim_theta
integer::dim_r1,dim_theta1
integer::i,j,k,l
real::t_start,t_end,t
real::var0,var1,var2,var3,var4,var5,var6
real::var2s,var3s,var5s,var6s
real::G
real::dr,dtheta
real::dforce_r,dforce_theta
real,dimension(256+1,256)::cos_dtheta 
real,dimension(256+1,256)::sin_dtheta 
real,dimension(128+1,128)::r_sub
real,dimension(128+1,128)::dr_sub 
real,dimension(128)::r
real,dimension(256)::theta
real,dimension(128+1)::r1
real,dimension(256+1)::theta1
real,dimension(128,256)::density
real,dimension(128+1,256+1,2)::force
G=1
dim_theta=256
dim_r=128
dim_theta1=dim_theta+1
dim_r1=dim_r+1
!_______________________________________________________________________
			
! load data from Computational_Science_II/Data/_________________________
open(unit=1,file="Data/r_project.data")
do i=1,dim_r
	read(1,'(e20.10)') r(i)
end do
close(1)
open(unit=2,file="Data/theta_project.data")
do i=1,dim_theta
	read(2,'(e20.10)') theta(i)
end do
close(2)
open(unit=3,file="Data/density_project.data")
do i=1,dim_r
	do k=1,dim_theta
		read(3,'(e20.10)') density(i,k)
	end do
end do
close(3)
!_______________________________________________________________________

! build theta1, cos_dtheta and sin_dtheta grid__________________________
dtheta=theta(2)-theta(1)
do j=1,dim_theta
	theta1(j)=theta(j)-dtheta*0.5
end do
theta1(dim_theta1)=theta(dim_theta)+dtheta*0.5
do j=1,dim_theta1
	do l=1,dim_theta
		var0=theta1(j)-theta(l)
		cos_dtheta(j,l)=cos(var0)
		sin_dtheta(j,l)=sin(var0)
	end do
end do
!_______________________________________________________________________

! build r1, r_sub and dr_sub grid_______________________________________
dr=r(2)-r(1)
do i=1,dim_r
	r1(i)=r(i)-dr*0.5
end do
r1(dim_r1)=r(dim_r)+dr*0.5
do i=1,dim_r1
	do k=1,dim_r
		var1=r(k)
		r_sub(i,k)=r1(i)/var1
		dr_sub(i,k)=-r_sub(i,k)/var1*dr
	end do
end do
!_______________________________________________________________________

call CPU_Time(t_start)	

! calculate force_______________________________________________________
do j=1,dim_theta1
	do i=1,dim_r1
		dforce_r=0.0
		dforce_theta=0.0
		do l=1,dim_theta
			
			
			do k=1,dim_r
				
				var5=r1(i)*r1(i)+r(k)*r(k)-2*r(k)*r1(i)*cos(theta1(j)-theta(l))
				
				var5=var5*sqrt(var5)
			
				var6=density(k,l)*r(k)
				
				var6=var6/var5
				
				dforce_r=dforce_r+(r1(i)-r(k)*cos(theta1(j)-theta(l)))*var6
				dforce_theta=dforce_theta+r(k)*sin(theta1(j)-theta(l))*var6
			end do
		end do
		force(i,j,1)=dforce_r*dtheta*dr*G
		force(i,j,2)=dforce_theta*dtheta*dr*G
	end do
end do
!_______________________________________________________________________

call CPU_Time(t_end)				
t=t_end-t_start
write(*,*) "required time:"
write(*,*) t

! save force_project.data_______________________________________________
open(unit=4,file="Results/force_project.data")
do i=1,dim_r
	do j=1,dim_theta
		write(4,'(2e20.10)') force(i,j,1) , force(i,j,2)
	end do
end do
close(4)	
!_______________________________________________________________________

end program
