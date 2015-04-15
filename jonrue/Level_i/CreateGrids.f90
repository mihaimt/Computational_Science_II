subroutine create_empty_DensityGrid(dim_r, LevelNr, r, dim_r_out, r_out)
implicit none
integer, intent(in) :: dim_r
integer, intent(in) :: LevelNr
integer, intent(in) :: dim_r_out
real, dimension(dim_r), intent(in) :: r
real, dimension(dim_r_out), intent(out) :: r_out
integer :: iteration_variable
integer :: step_size
step_size=dim_r/dim_r_out
do iteration_variable=1,dim_r_out
	r_out(iteration_variable)=(r(1+(iteration_variable-1)*step_size)+r(step_size+(iteration_variable-1)*step_size))/2
end do
end subroutine create_empty_DensityGrid

subroutine get_dimension(LevelNr,dim_r,dim_r_out)
integer, intent(in) :: LevelNr
integer, intent(in) :: dim_r
integer, intent(out) :: dim_r_out
integer :: temp_LevelNr
dim_r_out=dim_r
temp_LevelNr=LevelNr
do while(temp_LevelNr/=0)
	dim_r_out=dim_r_out/2
	temp_LevelNr=temp_LevelNr-1
end do
end subroutine get_dimension

subroutine create_real_DensityGrid(inputGrid, dim_IGR, dim_IGT, outputGrid,dim_OGR, dim_OGT)
implicit none
integer, intent(in) :: dim_IGR
integer, intent(in) :: dim_IGT
real, dimension(dim_IGR,dim_IGT) :: inputGrid
integer, intent(in) :: dim_OGR
integer, intent(in) :: dim_OGT
real, dimension(dim_OGR,dim_OGT), intent(out) :: outputGrid

integer :: i
integer :: j
integer :: k
integer :: l
real :: temp
integer :: stepsize
stepsize=dim_IGR/dim_OGR
		
do i=1,dim_OGR
	do j=1,dim_OGT
		temp=0
		do k=(1+(i-1)*stepsize),(i*stepsize)
			do l=(1+(j-1)*stepsize),(j*stepsize)
				
				temp=temp+inputGrid(k,l)
		
			end do
		end do
		outputGrid(i,j)=temp
	end do
end do

end subroutine create_real_DensityGrid
	

