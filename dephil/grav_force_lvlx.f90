!--------------------------------------------------------------------------------------
! PROGRAM GRAV_FORCE_LVLX
!	compile: with make (Makefile included), e.g. make grav_force_lvlx
!       
!	use:	 ./grav_force_lvlx [LEVEL], e.g.  ./grav_force_lvlx 1
!--------------------------------------------------------------------------------------
PROGRAM grav_force_lvlx
    IMPLICIT NONE

! DECLARATIONS
    INTEGER :: i, iprime, j, jprime
    INTEGER, PARAMETER :: N_r0=128, N_theta0=256
    INTEGER :: N_rx, N_thetax
    INTEGER :: LEVEL
    REAL(8)  :: epsilon2 ! this value seems to be the closest to the LEVEL0 solution
                         ! the value for LEVEL 2 seems to be close to 0.000336

    REAL(8) :: force_point, denom_point,&
             & f_rx=0., f_thetax=0.

    REAL(8), DIMENSION(N_r0) :: r0, dr0, r0_squared, r0_prime_squared
    REAL(8), DIMENSION(:), ALLOCATABLE :: rx, drx, rx_squared, rx_prime_squared
    REAL(8), DIMENSION(N_theta0) :: theta0, dtheta0
    REAL(8), DIMENSION(:), ALLOCATABLE :: thetax, dthetax
    REAL(8), DIMENSION(N_r0, N_r0) :: r0_ratio, r0_ratio_squared
    REAL(8), DIMENSION(:, :), ALLOCATABLE :: rx_ratio, rx_ratio_squared
    REAL(8), DIMENSION(N_r0, N_theta0) :: sigma0, mass0
    REAL(8), DIMENSION(:, :), ALLOCATABLE :: sigmax, massx
    REAL(8), DIMENSION(1-N_theta0:N_theta0-1) :: cos_table0, sin_table0
    REAL(8), DIMENSION(:, :), ALLOCATABLE :: cos_tablex, sin_tablex

    REAL :: start, finish

    CHARACTER(len=1) :: arg
!--------------------------------------------------------------------------------------
! READ THE LEVEL (FROM THE ARGUMENT LINE)

    call getarg(1, arg)
    read (arg, '(i5)') LEVEL

    N_rx = N_r0/2**LEVEL
    N_thetax = N_theta0/2**LEVEL

    if (LEVEL==0) then
        epsilon2 = 0
    else if (LEVEL==1) then
        epsilon2 = .000001 !.000118
    else if (LEVEL==2) then
        epsilon2 = .000001 !.000346
    else if (LEVEL==3) then
        epsilon2 = .000001 !.000902      ! try different values
    else
        epsilon2 = .000001       ! try different values
    end if

!--------------------------------------------------------------------------------------
! OPEN FILES
    open(unit=11, file=TRIM('./data/f_radial_lvl'//arg//'.data'), action='write')
    open(unit=12, file=TRIM('./data/f_angular_lvl'//arg//'.data'), action='write')

!--------------------------------------------------------------------------------------
! READ INPUT FILES
    call read_input(N_r0, N_theta0, r0, theta0, sigma0)

! READ IN LEVEL ONE ------------------
    call grid_param_lvlx(LEVEL, N_r0, N_theta0, r0, theta0, sigma0, rx, thetax, sigmax)

!--------------------------------------------------------------------------------------
! PRECALCULATIONS
    call precalcs_lvl0(N_r0, N_theta0, r0, theta0, dr0, dtheta0, r0_squared, r0_prime_squared,&
                       & r0_ratio, r0_ratio_squared, cos_table0, sin_table0)

    call precalcs_lvlx(LEVEL, N_r0, N_theta0, r0, theta0, rx, thetax,&
                                 & drx, dthetax, rx_squared, rx_prime_squared,&
                                 & rx_ratio, rx_ratio_squared,&
                                 & cos_tablex, sin_tablex)

!--------------------------------------------------------------------------------------
! PROPAGATION starts here
    call CPU_TIME(start)

! fill the <mass table> for LEVEL 0
    call calc_masslvl0(N_r0, N_theta0, sigma0, r0, dr0, dtheta0, mass0)

! fill the <mass table> for LEVEL 1
    call calc_masslvlx(LEVEL, N_r0, N_theta0, mass0, massx)

! write force components for every corner in grid
    do i = 1, N_r0
        do j = 1, N_theta0!, 2**LEVEL
            ! sum up the forces on the point (i, j)
            do iprime = 1, N_rx
                do jprime = 1, N_thetax
                    ! formula for the gravitational force split into four parts for faster calculation
                            ! expressed in terms of ratios of r
                            ! F_grav = Sum ( Mass / (r_iprime*sqrt(1+r_ratio_squared-2*r_ratio*cos))³/² ) * r_iprime | (r_ratio - cos), (sin)

                    !denom_point = InvSqrt(1+r1_ratio_squared(i, iprime)-2*r1_ratio(i, iprime)*cos_table1(j, jprime)+epsilon2)
                    denom_point = (1+rx_ratio_squared(i, iprime)-2*rx_ratio(i, iprime)*cos_tablex(j, jprime)+epsilon2)
                    !denom_point = c_invsqrt64(1+r1_ratio_squared(i, iprime)-2*r1_ratio(i, iprime)*cos_table1(j, jprime)+epsilon2)   ! returns somehow Nan's

                    force_point = massx(iprime,jprime)/(sqrt(denom_point)*denom_point*rx_prime_squared(iprime))

                    f_rx = f_rx + force_point&
                            &*(rx_ratio(i, iprime)-cos_tablex(j, jprime))
                    f_thetax = f_thetax + force_point&
                            &*sin_tablex(j, jprime)
                end do
            end do
            !do jprime = 1, 2**LEVEL
                write(11, '(e20.10)') f_rx
                write(12, '(e20.10)') f_thetax
            !end do
            f_rx = 0.
            f_thetax = 0.
        end do
    end do

! here the calculation has ended
    call CPU_TIME(finish)
    write(*,*) 'Calculation time =', finish-start, 'seconds'

    close(unit=11)
    close(unit=12)

! difference of the forces
!    open(unit=13, file='./data/f_radial_diffto'//arg//'.data', action='write')
!    open(unit=14, file='./data/f_radial.data', action='read')
!    open(unit=15, file='./data/f_radial_lvl'//arg//'.data', action='read')

!    do i = 1, N_r0*N_theta0
!        read(14, '(e20.10)') f_r
!        read(15, '(e20.10)') f_r1
!        write(13, '(e20.10)') f_r-f_r1
!    end do

!    close(unit=13)
!    close(unit=14)
!    close(unit=15)

    CONTAINS

        SUBROUTINE read_input(N_r0, N_theta0, r0, theta0, sigma0)
            IMPLICIT NONE
            ! DECLARATIONS
            INTEGER, intent(in) :: N_r0, N_theta0
            REAL(8), DIMENSION(N_r0), intent(out) :: r0
            REAL(8), DIMENSION(N_theta0), intent(out) :: theta0
            REAL(8), DIMENSION(N_r0, N_theta0), intent(out) :: sigma0
            INTEGER :: i, j     ! iterators
            REAL(8) :: temp     ! temporary variable

            ! OPEN INPUT FILES
            open(unit=8, file='./data/r_project.data', status='old', action='read')
            open(unit=9, file='./data/theta_project.data', status='old', action='read')
            open(unit=10, file='./data/density_project.data', status='old', action='read')

            ! READ FILES
            ! read <radii> into 1-D array for LEVEL0
            do i = 1, N_r0
                read(8, '(e20.10)') temp
                r0(i) = temp
            end do
            close(unit=8)
            ! read <angles> into 1-D array for LEVEL0
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

            REAL(8), DIMENSION(:), ALLOCATABLE, intent(out) :: rlvl         ! allocate according to level
            REAL(8), DIMENSION(:), ALLOCATABLE, intent(out) :: thetalvl     ! allocate according to level
            REAL(8), DIMENSION(:,:), ALLOCATABLE, intent(out) :: sigmalvl   ! allocate according to level

            INTEGER :: i, j, permute_i, permute_j     ! iterators
            INTEGER :: lvl_factor                     ! = 2^N_lvl
            INTEGER :: N_rlvl, N_thetalvl             ! dimensions of parameters in Level N_lvl
            REAL(8) :: temp=0

            lvl_factor=2**N_lvl
            N_rlvl=N_r0/lvl_factor
            N_thetalvl=N_theta0/lvl_factor

            allocate(rlvl(N_rlvl))
            allocate(thetalvl(N_thetalvl))
            allocate(sigmalvl(N_rlvl, N_thetalvl))

            ! READ IN LEVEL N_lvl ------------------
            ! read <radii> into 1-D array for LEVEL X
            do i = lvl_factor, N_r0, lvl_factor
                rlvl(i/lvl_factor) = (r0(i)+r0(i-(lvl_factor-1)))*.5
            end do
            ! read <angles> into 1-D array for LEVEL X
            do j = lvl_factor, N_theta0, lvl_factor
                thetalvl(j/lvl_factor) = (theta0(j)+theta0(j-(lvl_factor-1)))*.5
            end do
            ! read <density> into 2-D array for LEVEL X
            do i = lvl_factor, N_r0, lvl_factor
                do j = lvl_factor, N_theta0, lvl_factor
                    do permute_i = 0, lvl_factor-1
                        do permute_j = 0, lvl_factor-1
                            temp = temp + sigma0(i-permute_i, j-permute_j)
                        end do
                    end do
                sigmalvl(i/lvl_factor, j/lvl_factor) = temp
                temp = 0
                end do
            end do

        END SUBROUTINE grid_param_lvlx

        SUBROUTINE precalcs_lvl0(N_r0, N_theta0, r0, theta0, dr0, dtheta0, r0_squared, r0_prime_squared,&
                                 & r0_ratio, r0_ratio_squared, cos_table0, sin_table0)
            IMPLICIT NONE
            ! DECLARATIONS
            INTEGER, intent(in) :: N_r0, N_theta0
            REAL(8), DIMENSION(N_r0), intent(in) :: r0
            REAL(8), DIMENSION(N_theta0), intent(in) :: theta0
            REAL(8), DIMENSION(N_r0), intent(out) :: dr0, r0_squared, r0_prime_squared
            REAL(8), DIMENSION(N_r0, N_r0), intent(out) :: r0_ratio, r0_ratio_squared
            REAL(8), DIMENSION(N_theta0), intent(out) :: dtheta0
            REAL(8), DIMENSION(1-N_theta0:N_theta0-1), intent(out) :: cos_table0, sin_table0
            INTEGER :: i, j, iprime, jprime     ! iterators
            REAL(8) :: r0_i, r0_iprime, theta0_j, diff_theta

            ! different <dr/dtheta> for LEVEL 0 ------------------
            do i = 1, N_r0-1
                dr0(i) = r0(i+1)-r0(i)
            end do
            dr0(N_r0) = dr0(N_r0-1)
            do j = 1, N_theta0-1
                dtheta0(j) = theta0(j+1)-theta0(j)
            end do
            dtheta0(N_theta0) = dtheta0(N_theta0-1)
            ! <r²> in the corners and <r'²> in the centres for LEVEL0 ---------
            do i = 1, N_r0
                r0_i = r0(i)-.5*dr0(i) ! shift r to the corners
                r0_squared(i) = r0_i*r0_i
                r0_prime_squared(i) = r0(i)*r0(i)
            end do
            ! <r ratios> and <r ratios squared> for LEVEL0 ------------------
            do i = 1, N_r0
                r0_i = r0(i)-.5*dr0(i) ! shift r to the corners
                do iprime = 1, N_r0
                    r0_iprime = r0_i/r0(iprime) ! r_iprime only used as temp
                    r0_ratio(i, iprime) = r0_iprime
                    r0_ratio_squared(i, iprime) = r0_iprime*r0_iprime
                end do
            end do
            ! fill the <sin/cos table> for LEVEL 0 ------------------
            do j = 1, N_theta0
                theta0_j = theta0(j)-.5*dtheta0(j)
                do jprime = 1, N_theta0
                    diff_theta = theta0_j-theta0(jprime)
                    cos_table0(j-jprime) = cos(diff_theta)
                    sin_table0(j-jprime) = sin(diff_theta)
                end do
            end do

        END SUBROUTINE precalcs_lvl0

        SUBROUTINE precalcs_lvlx(N_lvl, N_r0, N_theta0, r0, theta0, rlvl, thetalvl,&
                                 & drlvl, dthetalvl, rlvl_squared, rlvl_prime_squared,&
                                 & rlvl_ratio, rlvl_ratio_squared,&
                                 & cos_tablelvl, sin_tablelvl)
                                 ! allocatable: rlvl, thetalvl, drlvl, rlvl_squared, rlvl_prime_squared, dthetalvl, rlvl_ratio, rlvl_ratio_squared, cos_tablelvl, sin_tablelvl
            IMPLICIT NONE
            ! DECLARATIONS
            INTEGER, intent(in) :: N_lvl
            INTEGER, intent(in) :: N_r0, N_theta0
            REAL(8), DIMENSION(N_r0), intent(in) :: r0
            REAL(8), DIMENSION(N_theta0), intent(in) :: theta0
            REAL(8), DIMENSION(:), ALLOCATABLE, intent(in) :: rlvl         ! at this point already allocated
            REAL(8), DIMENSION(:), ALLOCATABLE, intent(in) :: thetalvl     ! at this point already allocated

            REAL(8), DIMENSION(:), ALLOCATABLE, intent(out) :: drlvl, rlvl_squared, rlvl_prime_squared
            REAL(8), DIMENSION(:), ALLOCATABLE, intent(out) :: dthetalvl
            REAL(8), DIMENSION(:, :), ALLOCATABLE, intent(out) :: rlvl_ratio, rlvl_ratio_squared
            REAL(8), DIMENSION(:, :), ALLOCATABLE, intent(out) :: cos_tablelvl, sin_tablelvl

            INTEGER :: i, j                           ! iterators
            INTEGER :: lvl_factor                     ! = 2^N_lvl
            INTEGER :: N_rlvl, N_thetalvl             ! dimensions of parameters in Level N_lvl
            REAL(8) :: r0_i, rlvl_i, rlvl_iprime, theta0_j, diff_theta

            lvl_factor=2**N_lvl
            N_rlvl=N_r0/lvl_factor
            N_thetalvl=N_theta0/lvl_factor

            allocate(drlvl(N_rlvl), rlvl_squared(N_rlvl), rlvl_prime_squared(N_rlvl))
            allocate(dthetalvl(N_thetalvl))
            allocate(rlvl_ratio(N_r0, N_rlvl), rlvl_ratio_squared(N_r0, N_rlvl))
            allocate(cos_tablelvl(N_theta0, N_thetalvl), sin_tablelvl(N_theta0, N_thetalvl))

            ! different <dr/dtheta> for LEVEL X ------------------
            do i = 1, N_rlvl-1
                drlvl(i) = rlvl(i+1)-rlvl(i)
            end do
            drlvl(N_rlvl) = drlvl(N_rlvl-1)
            do j = 1, N_thetalvl-1
                dthetalvl(j) = thetalvl(j+1)-thetalvl(j)
            end do
            dthetalvl(N_thetalvl) = dthetalvl(N_thetalvl-1)

            ! <r²> in the corners and <r'²> in the centres for LEVEL X ---------
            do i = 1, N_rlvl
                rlvl_i = rlvl(i)-.5*drlvl(i) ! shift r to the corners
                rlvl_squared(i) = rlvl_i*rlvl_i
                rlvl_prime_squared(i) = rlvl(i)*rlvl(i)
            end do

            ! <r ratios> and <r ratios squared> for LEVEL X ------------------
            do i = 1, N_r0
                r0_i = r0(i)-.5*dr0(i) ! shift r to the corners
                do iprime = 1, N_rlvl
                    rlvl_iprime = r0_i/rlvl(iprime) ! r_iprime only used as temp
                    rlvl_ratio(i, iprime) = rlvl_iprime
                    rlvl_ratio_squared(i, iprime) = rlvl_iprime*rlvl_iprime
                end do
            end do

            ! fill the <sin/cos table> for LEVEL X ------------------
            do j = 1, N_theta0
                theta0_j = theta0(j)-.5*dtheta0(j)
                do jprime = 1, N_thetalvl
                    diff_theta = theta0_j-thetalvl(jprime)
                    cos_tablelvl(j, jprime) = cos(diff_theta)
                    sin_tablelvl(j, jprime) = sin(diff_theta)
                end do
            end do

        END SUBROUTINE precalcs_lvlx

        SUBROUTINE calc_masslvl0(N_r0, N_theta0, sigma0, r0, dr0, dtheta0, mass0)
            IMPLICIT NONE
            ! DECLARATION
            INTEGER, intent(in) :: N_r0, N_theta0
            REAL(8), DIMENSION(N_r0), intent(in) :: r0, dr0
            REAL(8), DIMENSION(N_theta0), intent(in) :: dtheta0
            REAL(8), DIMENSION(N_r0, N_theta0), intent(in) :: sigma0
            REAL(8), DIMENSION(N_r0, N_theta0), intent(out) :: mass0
            INTEGER :: iprime, jprime
            REAL(8) :: r0dr0_iprime

            ! fill the <mass table> for LEVEL 0
            do iprime = 1, N_r0
                r0dr0_iprime = r0(iprime)*dr0(iprime) ! multiply dr already here to save some operations
                do jprime = 1, N_theta0
                    mass0(iprime,jprime) = -sigma0(iprime,jprime)*r0dr0_iprime*dtheta0(jprime)
                end do
            end do

        END SUBROUTINE calc_masslvl0

        SUBROUTINE calc_masslvlx(N_lvl, N_r0, N_theta0, mass0, masslvl)
            IMPLICIT NONE
            ! DECLARATION
            INTEGER, intent(in) :: N_lvl
            INTEGER, intent(in) :: N_r0, N_theta0
            REAL(8), DIMENSION(N_r0, N_theta0), intent(in) :: mass0
            REAL(8), DIMENSION(:, :), ALLOCATABLE, intent(out) :: masslvl

            INTEGER :: i, j, permute_i, permute_j     ! iterators
            INTEGER :: lvl_factor                     ! = 2^N_lvl
            INTEGER :: N_rlvl, N_thetalvl             ! dimensions of parameters in Level N_lvl
            REAL(8) :: temp=0

            lvl_factor=2**N_lvl
            N_rlvl=N_r0/lvl_factor
            N_thetalvl=N_theta0/lvl_factor

            allocate(masslvl(N_rlvl, N_thetalvl))

            ! fill the <mass table> for LEVEL 1
            do i = lvl_factor, N_r0, lvl_factor
                do j = lvl_factor, N_theta0, lvl_factor
                    do permute_i = 0, lvl_factor-1
                        do permute_j = 0, lvl_factor-1
                            temp = temp + mass0(i-permute_i, j-permute_j)
                        end do
                    end do
                masslvl(i/lvl_factor, j/lvl_factor) = temp
                temp = 0
                end do
            end do

        END SUBROUTINE

END PROGRAM grav_force_lvlx