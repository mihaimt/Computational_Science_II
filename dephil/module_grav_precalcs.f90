!----------------------------------------------------------------------
! MODULE module_grav_precalcs
!
!     subroutine precalcs_lvl0         precalculations to r0, theta0 and density0
!     subroutine precalcs_lvlx         precalculations of a given level N_lvl
!     subroutine precalcs_denomx       not used in main program; precalculation of the sqrt in the denominator of the force calculation to a given level N_lvl
!
!----------------------------------------------------------------------
MODULE grav_precalcs

    CONTAINS

        SUBROUTINE precalcs_lvl0(N_r0, N_theta0, r0, theta0, dr0, dtheta0, r0_squared, r0_prime_squared,&
                         & r0_ratio, r0_ratio_squared, cos_table0, sin_table0)
            IMPLICIT NONE
            ! DECLARATIONS
            INTEGER, intent(in) :: N_r0, N_theta0
            REAL(8), DIMENSION(:), ALLOCATABLE, intent(in) :: r0
            REAL(8), DIMENSION(:), ALLOCATABLE, intent(in) :: theta0
            REAL(8), DIMENSION(:), ALLOCATABLE, intent(out) :: dr0
            REAL(8), DIMENSION(0:N_r0+1, 0:N_r0+1), intent(out) :: r0_ratio, r0_ratio_squared
            REAL(8), DIMENSION(:), ALLOCATABLE, intent(out) :: dtheta0
            REAL(8), DIMENSION(1-N_theta0:N_theta0-1), intent(out) :: cos_table0, sin_table0
            REAL(8), DIMENSION(:), ALLOCATABLE, intent(out) :: r0_squared, r0_prime_squared
            ! Local variables
            INTEGER :: i, j, iprime, jprime     ! iterators
            REAL(8) :: r0_i, r0_iprime, theta0_j, diff_theta

            ! different <dr/dtheta> for LEVEL 0   - distances from centers of cells == cell widths
            allocate(dr0(-1:N_r0+2))
            do i = 1, N_r0-1
                dr0(i) = r0(i+1)-r0(i)
            end do
            dr0(N_r0) = dr0(N_r0-1)
            ! add ghost cells
            dr0(0) = dr0(1)
            dr0(-1) = dr0(0)
            dr0(N_r0+1) = dr0(N_r0)
            dr0(N_r0+2) = dr0(N_r0+1)
            
            allocate(dtheta0(-1:N_theta0+2))
            do j = 1, N_theta0-1
                dtheta0(j) = theta0(j+1)-theta0(j)
            end do
            dtheta0(N_theta0) = dtheta0(N_theta0-1)
            ! add ghost cells
            dtheta0(0) = dtheta0(N_theta0)
            dtheta0(-1) = dtheta0(N_theta0-1)
            dtheta0(N_theta0+1) = dtheta0(1)
            dtheta0(N_theta0+2) = dtheta0(2)

            ! <r²> in the corners and <r'²> in the centres for LEVEL0   - r² indicates corners and r'² indicates centers of the cells
            allocate(r0_squared(N_r0), r0_prime_squared(N_r0))
            do i = 1, N_r0
                r0_i = r0(i)-.5*dr0(i) ! shift r to the corners
                r0_squared(i) = r0_i*r0_i
                r0_prime_squared(i) = r0(i)*r0(i)
            end do

            ! <r ratios> and <r ratios squared> for LEVEL0   - ratios of r to r' from corners to centers of the cells
            do i = 1, N_r0
                r0_i = r0(i)-.5*dr0(i) ! shift r to the corners
                do iprime = 1, N_r0
                    r0_iprime = r0_i/r0(iprime) ! r_iprime only used as temp
                    r0_ratio(i, iprime) = r0_iprime
                    r0_ratio_squared(i, iprime) = r0_iprime*r0_iprime
                end do
            end do

            ! the <sin/cos> tables for LEVEL 0   - of differences of thetas from corner to centers of the cells
            do j = 1, N_theta0
                theta0_j = theta0(j)-.5*dtheta0(j)
                do jprime = 1, N_theta0
                    diff_theta = theta0_j-theta0(jprime)
                    cos_table0(j-jprime) = cos(diff_theta)
                    sin_table0(j-jprime) = sin(diff_theta)
                end do
            end do
        END SUBROUTINE precalcs_lvl0


        SUBROUTINE precalcs_lvlx(N_lvl, N_r0, N_theta0, r0, theta0, dr0, dtheta0, rlvl, thetalvl,&
                         & drlvl, dthetalvl, rlvl_squared, rlvl_prime_squared,&
                         & rlvl_ratio, rlvl_ratio_squared,&
                         & cos_tablelvl, sin_tablelvl)
                         ! allocatable: rlvl, thetalvl, drlvl, rlvl_squared, rlvl_prime_squared, dthetalvl, rlvl_ratio, rlvl_ratio_squared, cos_tablelvl, sin_tablelvl
            IMPLICIT NONE
            ! DECLARATIONS
            INTEGER, intent(in) :: N_lvl
            INTEGER, intent(in) :: N_r0, N_theta0
            REAL(8), DIMENSION(:), ALLOCATABLE, intent(in) :: r0, dr0
            REAL(8), DIMENSION(:), ALLOCATABLE, intent(in) :: theta0, dtheta0
            REAL(8), DIMENSION(:), ALLOCATABLE, intent(in) :: rlvl         ! at this point already allocated
            REAL(8), DIMENSION(:), ALLOCATABLE, intent(in) :: thetalvl     ! at this point already allocated

            REAL(8), DIMENSION(:), ALLOCATABLE, intent(out) :: drlvl, rlvl_squared, rlvl_prime_squared    ! allocated according to level
            REAL(8), DIMENSION(:), ALLOCATABLE, intent(out) :: dthetalvl                                  ! allocated according to level
            REAL(8), DIMENSION(:, :), ALLOCATABLE, intent(out) :: rlvl_ratio, rlvl_ratio_squared          ! allocated according to level
            REAL(8), DIMENSION(:, :), ALLOCATABLE, intent(out) :: cos_tablelvl, sin_tablelvl              ! allocated according to level

            INTEGER :: i, j, iprime, jprime           ! iterators
            INTEGER :: factor_lvl                     ! = 2^N_lvl
            INTEGER :: N_rlvl, N_thetalvl             ! dimensions of parameters in Level N_lvl
            REAL(8) :: r0_i, rlvl_i, rlvl_iprime, theta0_j, diff_theta

            factor_lvl=2**N_lvl
            N_rlvl=N_r0/factor_lvl
            N_thetalvl=N_theta0/factor_lvl

            allocate(drlvl(-1:N_rlvl+2), rlvl_squared(0:N_rlvl+1), rlvl_prime_squared(0:N_rlvl+1))
            allocate(dthetalvl(-1:N_thetalvl+2))
            allocate(rlvl_ratio(N_r0, 0:N_rlvl+1), rlvl_ratio_squared(N_r0, 0:N_rlvl+1))
            allocate(cos_tablelvl(N_theta0, 0:N_thetalvl+1), sin_tablelvl(N_theta0, 0:N_thetalvl+1))

            ! different <dr/dtheta> for LEVEL X   - distances from centers of cells == cell widths
            do i = 0, N_rlvl
                drlvl(i) = rlvl(i+1)-rlvl(i)
            end do
            ! add ghost cells
            drlvl(-1) = drlvl(0)
            drlvl(N_rlvl+1) = drlvl(N_rlvl)
            drlvl(N_rlvl+2) = drlvl(N_rlvl+1)

            do j = 1, N_thetalvl-1
                dthetalvl(j) = thetalvl(j+1)-thetalvl(j)
            end do
            dthetalvl(N_thetalvl) = dthetalvl(N_thetalvl-1)
            ! add ghost cells
            dthetalvl(0) = dthetalvl(N_thetalvl)
            dthetalvl(-1) = dthetalvl(N_thetalvl-1)
            dthetalvl(N_thetalvl+1) = dthetalvl(1)
            dthetalvl(N_thetalvl+2) = dthetalvl(2)

            ! <r²> in the corners and <r'²> in the centres for LEVEL X   - r² indicates corners and r'² indicates centers of the cells
            do i = 0, N_rlvl+1
                rlvl_i = rlvl(i)-.5*drlvl(i) ! shift r to the corners
                rlvl_squared(i) = rlvl_i*rlvl_i
                rlvl_prime_squared(i) = rlvl(i)*rlvl(i)
            end do

            ! <r ratios> and <r ratios squared> for LEVEL X   - ratios of r to r' from corners to centers of the cells
            do i = 1, N_r0
                r0_i = r0(i)-.5*dr0(i) ! shift r to the corners
                do iprime = 0, N_rlvl+1
                    rlvl_iprime = r0_i/rlvl(iprime) ! r_iprime only used as temp
                    rlvl_ratio(i, iprime) = rlvl_iprime
                    rlvl_ratio_squared(i, iprime) = rlvl_iprime*rlvl_iprime
                end do
            end do

            ! fill the <sin/cos table> for LEVEL X   - of differences of thetas from corner to centers of the cells
            do j = 1, N_theta0
                theta0_j = theta0(j)-.5*dtheta0(j)
                do jprime = 0, N_thetalvl+1
                    diff_theta = theta0_j-thetalvl(jprime)
                    cos_tablelvl(j, jprime) = cos(diff_theta)
                    sin_tablelvl(j, jprime) = sin(diff_theta)
                end do
            end do

        END SUBROUTINE precalcs_lvlx


        SUBROUTINE precalc_denomx(N_lvl, N_r0, N_theta0, rlvl_ratio,&
                          & rlvl_ratio_squared, cos_tablelvl, denomlvl) ! allocatable: rlvl_ratio, rlvl_ration_squared, cos_tablelvl, denomlvl
            IMPLICIT NONE
            ! DECLARATIONS
            INTEGER, intent(in) :: N_lvl, N_r0, N_theta0
            REAL(8), DIMENSION(:, :), ALLOCATABLE, intent(in) :: rlvl_ratio, rlvl_ratio_squared
            REAL(8), DIMENSION(:, :), ALLOCATABLE, intent(in) :: cos_tablelvl
            REAL(8), DIMENSION(:, :, :, :), ALLOCATABLE, intent(out) :: denomlvl

            INTEGER :: i, j, iprime, jprime           ! iterators
            INTEGER :: factor_lvl                     ! = 2^N_lvl
            INTEGER :: N_rlvl, N_thetalvl             ! dimensions of parameters in Level N_lvl

            factor_lvl=2**N_lvl
            N_rlvl=N_r0/factor_lvl
            N_thetalvl=N_theta0/factor_lvl

            allocate(denomlvl(N_r0, N_theta0, N_rlvl, N_thetalvl))

            do i = 1, N_r0
                do j = 1, N_theta0
                    do iprime = 1, N_rlvl
                        do jprime = 1, N_thetalvl
                            denomlvl(i, j, iprime, jprime) = sqrt(1+rlvl_ratio_squared(i, iprime)&
                                                                 &-2*rlvl_ratio(i, iprime)*cos_tablelvl(j, jprime))
                        end do
                    end do
                end do
            end do
        END SUBROUTINE precalc_denomx

END MODULE grav_precalcs
