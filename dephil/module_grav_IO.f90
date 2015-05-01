!------------------------------------------------------------------------
! MODULE module_grav_IO
!
!     subroutine read_input                    reads input files r_project, theta_project and density_project and stores them in arrays of given dimension
!     subroutine grid_param_lvlx               creates the higher level parameters to r0, theta0 and density0 and stores them in arrays of given dimension to a given level N_lvl
!
!------------------------------------------------------------------------ 
MODULE grav_IO

    CONTAINS

        SUBROUTINE read_input(N_r0, N_theta0, r0, theta0, sigma0)
            IMPLICIT NONE
            ! DECLARATIONS
            INTEGER, intent(in) :: N_r0, N_theta0
            REAL(8), DIMENSION(N_r0), intent(out) :: r0
            REAL(8), DIMENSION(N_theta0), intent(out) :: theta0
            REAL(8), DIMENSION(N_r0, N_theta0), intent(out) :: sigma0
            ! Local Variables
            INTEGER :: i, j     ! iterators
            REAL(8) :: temp     ! temporary variable

            ! OPEN INPUT FILES
            open(unit=8, file='./data/r_project.data', status='old', action='read')
            open(unit=9, file='./data/theta_project.data', status='old', action='read')
            open(unit=10, file='./data/density_project.data', status='old', action='read')

            ! READ FILES
            ! read <radii> into 1-D array for LEVEL0   - indicates center of the cells
            do i = 1, N_r0
                read(8, '(e20.10)') temp
                r0(i) = temp
            end do
            close(unit=8)
            ! read <angles> into 1-D array for LEVEL0   - indicates center of the cells
            do j = 1, N_theta0
                read(9, '(e20.10)') temp
                theta0(j) = temp
            end do
            close(unit=9)
            ! read <density> into 2-D array for LEVEL0
            do i = 1, N_r0
                do j = 1, N_theta0
                    read(10, '(e20.10)') temp
                    sigma0(i,j) = temp
                end do
            end do
            close(unit=10)

        END SUBROUTINE read_input


        SUBROUTINE grid_param_lvlx(N_lvl, N_r0, N_theta0, r0, theta0, sigma0, rlvl, thetalvl, sigmalvl)    ! allocatable: rlvl, thetalvl, sigmalvl
            IMPLICIT NONE
            ! DECLARATIONS
            INTEGER, intent(in) :: N_lvl
            INTEGER, intent(in) :: N_r0, N_theta0
            REAL(8), DIMENSION(N_r0), intent(in) :: r0
            REAL(8), DIMENSION(N_theta0), intent(in) :: theta0
            REAL(8), DIMENSION(N_r0, N_theta0), intent(in) :: sigma0

            REAL(8), DIMENSION(:), ALLOCATABLE, intent(out) :: rlvl         ! allocated according to level
            REAL(8), DIMENSION(:), ALLOCATABLE, intent(out) :: thetalvl     ! allocated according to level
            REAL(8), DIMENSION(:,:), ALLOCATABLE, intent(out) :: sigmalvl   ! allocated according to level
            ! Local Variables
            INTEGER :: i, j, permute_i, permute_j     ! iterators
            INTEGER :: factor_lvl                     ! = 2^N_lvl
            INTEGER :: N_rlvl, N_thetalvl             ! dimensions of parameters in Level N_lvl
            REAL(8) :: temp=0

            factor_lvl=2**N_lvl
            N_rlvl=N_r0/factor_lvl
            N_thetalvl=N_theta0/factor_lvl

            allocate(rlvl(N_rlvl))
            allocate(thetalvl(N_thetalvl))
            allocate(sigmalvl(N_rlvl, N_thetalvl))

            ! READ IN LEVEL N_lvl ------------------
            ! read <radii> into 1-D array for LEVEL X   - indicates center of the cells
            do i = factor_lvl, N_r0, factor_lvl
                rlvl(i/factor_lvl) = (r0(i)+r0(i-(factor_lvl-1)))*.5
            end do
            ! read <angles> into 1-D array for LEVEL X   - indicates center of the cells
            do j = factor_lvl, N_theta0, factor_lvl
                thetalvl(j/factor_lvl) = (theta0(j)+theta0(j-(factor_lvl-1)))*.5
            end do
            ! read <density> into 2-D array for LEVEL X
            do i = factor_lvl, N_r0, factor_lvl
                do j = factor_lvl, N_theta0, factor_lvl
                    do permute_i = 0, factor_lvl-1
                        do permute_j = 0, factor_lvl-1
                            temp = temp + sigma0(i-permute_i, j-permute_j)
                        end do
                    end do
                sigmalvl(i/factor_lvl, j/factor_lvl) = temp
                temp = 0
                end do
            end do

        END SUBROUTINE grid_param_lvlx

END MODULE grav_IO
