subroutine sinus_cosinusPreCalculations(dim_theta, dim_theta_zoom, theta, theta_zoom, cosinusList, sinusList)
implicit none
integer, intent(in) :: dim_theta
integer, intent(in) :: dim_theta_zoom
real, dimension(dim_theta), intent(in) :: theta
real, dimension(dim_theta_zoom), intent(in) :: theta_zoom
real, dimension(dim_theta,dim_theta_zoom), intent(out) :: cosinusList
real, dimension(dim_theta,dim_theta_zoom), intent(out) :: sinusList
integer :: i
integer :: j

do i=1,dim_theta
	do j=1,dim_theta_zoom
		cosinusList(i,j)=cos(theta(i)-theta_zoom(j))
		sinusList(i,j)=sin(theta(i)-theta_zoom(j))
	end do
end do

end subroutine sinus_cosinusPreCalculations

subroutine R_dR_PreCalculations(dim_r, dim_r_zoom, r, r_zoom, R_List, dR_List)
implicit none 
integer, intent(in) :: dim_r
integer, intent(in) :: dim_r_zoom
real, dimension(dim_r), intent(in) :: r
real, dimension(dim_r_zoom), intent(in) :: r_zoom
real, dimension(dim_r,dim_r_zoom), intent(out) :: R_List
real, dimension(dim_r,dim_r_zoom), intent(out) :: dR_List
integer :: i
integer :: j
do i=1,dim_r
	do j=1,dim_r_zoom
		R_List(i,j)=r(i)/r_zoom(j)
		dR_List(i,j)=-r(i)/(r_zoom(j)*r_zoom(j))
	end do
end do
end subroutine R_dR_PreCalculations
