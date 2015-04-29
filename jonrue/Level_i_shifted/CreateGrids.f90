subroutine create_empty_DensityGrid(dim_r, LevelNr, r, dim_r_out, r_out)
implicit none
integer, intent(in) :: dim_r
integer, intent(in) :: LevelNr
integer, intent(in) :: dim_r_out
real, dimension(dim_r), intent(in) :: r
real, dimension(dim_r_out), intent(out) :: r_out
integer :: ind
integer :: zaehler
zaehler=1

do ind=1,dim_r-1,2
	r_out(zaehler)=(r(ind)+r(ind+1))*0.5
	zaehler=zaehler+1
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

integer :: ii
integer :: jj
integer :: ind_r
integer :: ind_t

ind_r=1
ind_t=1
do ii=1,dim_IGT-1,2
	ind_r=1
	do jj=1,dim_IGR-1,2
		outputGrid(ind_r,ind_t)=inputGrid(jj,ii)+inputGrid(jj+1,ii)+inputGrid(jj,ii+1)+inputGrid(jj+1,ii+1)
		ind_r=ind_r+1
	end do
	ind_t=ind_t+1
end do

end subroutine create_real_DensityGrid
	

