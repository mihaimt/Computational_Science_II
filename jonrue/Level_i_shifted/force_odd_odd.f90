subroutine force_odd_odd(r,theta,dim_r, dim_theta, r_achse, theta_achse,dim_theta_prime,&
			dim_r_prime, density_grid, f1, f2, r_prime, theta_prime)

implicit none
integer, intent(in) :: r
integer, intent(in) :: theta
integer, intent(in) :: dim_theta
integer, intent(in) :: dim_r
integer, intent(in) :: dim_theta_prime
integer, intent(in) :: dim_r_prime
real,dimension(dim_r_prime,dim_theta_prime), intent(in) :: density_grid
real,dimension(dim_r_prime), intent(in) :: r_prime
real,dimension(dim_theta_prime), intent(in) :: theta_prime
real,dimension(dim_r), intent(in) :: r_achse
real,dimension(dim_theta), intent(in) :: theta_achse
real, dimension(dim_r,dim_theta), intent(out) :: f1
real, dimension(dim_r,dim_theta), intent(out) :: f2


integer :: index_thetaPrime
integer :: index_rPrime

real :: r_component
real :: theta_component
real :: R_temp, cos_temp, sin_temp
real :: numerator
real :: denominator
real :: fraction_value

real :: t_p
real :: r_p
real :: dR
real :: dens


real :: surpress_divergent
surpress_divergent=1e-6

r_component=0
theta_component=0
do index_thetaPrime=1,dim_theta_prime,1
	t_p=theta_prime(index_thetaPrime)
	cos_temp=cos(theta_achse(theta)-t_p)
	sin_temp=sin(theta_achse(theta)-t_p)
	do index_rPrime=1,dim_r_prime,1
		r_p=r_prime(index_rPrime)
		R_temp=r_achse(r)/r_p
		dR=-r_achse(r)/(r_p*r_p)
		dens=density_grid(index_rPrime,index_thetaPrime)
		dens=dens*1
		numerator=dens*dR
		denominator=(1+R_temp*(R_temp-2*cos_temp))
		denominator=denominator*sqrt(denominator)*R_temp
		fraction_value=numerator/denominator
		r_component=(R_temp-cos_temp)*fraction_value+r_component
		theta_component=sin_temp*fraction_value+theta_component

	end do
end do
f1(r,theta)=r_component
f2(r,theta)=theta_component



end subroutine force_odd_odd
