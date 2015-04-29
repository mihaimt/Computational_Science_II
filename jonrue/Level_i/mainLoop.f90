subroutine mainLoop(dim_theta,dim_r,dim_theta_prime,dim_r_prime, density_grid, cosinus_List, sinus_List, R_List, dR_List, f1, f2)
implicit none
integer, intent(in) :: dim_theta
integer, intent(in) :: dim_r
integer, intent(in) :: dim_theta_prime
integer, intent(in) :: dim_r_prime
real, dimension(dim_r_prime,dim_theta_prime), intent(in) :: density_grid
real, dimension(dim_theta,dim_theta_prime), intent(in) :: cosinus_List
real, dimension(dim_theta,dim_theta_prime), intent(in) :: sinus_List
real, dimension(dim_r,dim_r_prime), intent(in) :: R_List
real, dimension(dim_r,dim_r_prime), intent(in) :: dR_List
real, dimension(dim_r,dim_theta), intent(out) :: f1
real, dimension(dim_r,dim_theta), intent(out) :: f2

integer :: index_theta
integer :: index_thetaPrime
integer :: index_r
integer :: index_rPrime

real :: r_component
real :: theta_component
real :: R_temp, cos_temp, sin_temp
real :: numerator
real :: denominator
real :: fraction_value

real :: surpress_divergent
surpress_divergent=1e-6

do index_theta=1,dim_theta
	do index_r=1,dim_r
		r_component=0
		theta_component=0
		do index_thetaPrime=1,dim_theta_prime
			cos_temp=cosinus_List(index_theta,index_thetaPrime)
			sin_temp=sinus_List(index_theta,index_thetaPrime)
			do index_rPrime=1,dim_r_prime
				R_temp=R_List(index_r,index_rPrime)
				numerator=density_grid(index_rPrime,index_thetaPrime)*dR_List(index_r,index_rPrime)
				denominator=(surpress_divergent+1+R_temp*(R_temp-2*cos_temp))
				denominator=denominator*sqrt(denominator)*R_temp
				fraction_value=numerator/denominator
				r_component=(R_temp-cos_temp)*fraction_value+r_component
				theta_component=sin_temp*fraction_value+theta_component

			end do
		end do
		f1(index_r,index_theta)=r_component
		f2(index_r,index_theta)=theta_component

	end do
end do

end subroutine mainLoop
