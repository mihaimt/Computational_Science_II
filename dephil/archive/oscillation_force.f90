!--------------------------------------------------------------------------------------
! PROGRAM OSCILLATION_FORCE
!	compile: with make (Makefile included), e.g. make oscillation_force
!       
!	use:	 ./oscillation_force [LEVEL], e.g.  ./oscillation_force 2
!--------------------------------------------------------------------------------------
PROGRAM oscillation_force
    USE grav_parameters
    USE grav_IO
    USE grav_precalcs
    USE grav_mass
    IMPLICIT NONE

! DECLARATIONS
    INTEGER :: i, iprime, j, jprime
    INTEGER :: N_rx, N_thetax
    INTEGER :: LEVEL, factor_x
    REAL(8) :: force_point, denom_point
    REAL :: start, finish

    CHARACTER(len=1) :: arg
!--------------------------------------------------------------------------------------
! READ THE LEVEL (FROM THE ARGUMENT LINE)

    call getarg(1, arg)
    read (arg, '(i5)') LEVEL

    factor_x = 2**LEVEL
    N_rx = N_r0/factor_x
    N_thetax = N_theta0/factor_x

    allocate(force_r(0:N_r0+factor_x, 0:N_theta0+factor_x), force_theta(0:N_r0+factor_x, 0:N_theta0+factor_x))     ! allocate with ghost cells

!--------------------------------------------------------------------------------------
! OPEN FILES
    open(unit=11, file=TRIM('./data/radial_osc_force_lvl'//arg//'.data'), action='write')
    open(unit=12, file=TRIM('./data/angular_osc_force_lvl'//arg//'.data'), action='write')

!--------------------------------------------------------------------------------------
! READ INPUT FILES
    call read_input(N_r0, N_theta0, r0, theta0, sigma0)

! READ IN LEVEL ONE ------------------
    call grid_param_lvlx(LEVEL, N_r0, N_theta0, r0, theta0, sigma0, rx, thetax, sigmax)

!--------------------------------------------------------------------------------------
! PRECALCULATIONS
    call precalcs_lvl0(N_r0, N_theta0, r0, theta0, dr0, dtheta0, r0_squared, r0_prime_squared,&
                       & r0_ratio, r0_ratio_squared, cos_table0, sin_table0)
    call precalcs_lvlx(LEVEL, N_r0, N_theta0, r0, theta0, dr0, dtheta0, rx, thetax,&
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

    write(*,*)"Calculating level "//arg//" without oscillation..."
! write force components for every corner in grid
    do i = 1, N_r0, factor_x
        do j = 1, N_theta0, factor_x
            ! sum up the forces on the point (i, j)
            do iprime = 1, N_rx
                do jprime = 1, N_thetax
                    ! formula for the gravitational force split into four parts for faster calculation expressed in terms of ratios of r
                    ! F_grav = Sum ( Mass / (r_iprime*sqrt(1+r_ratio_squared-2*r_ratio*cos))³/² ) * r_iprime | (r_ratio - cos), (sin)
                    denom_point = (1+rx_ratio_squared(i, iprime)-2*rx_ratio(i, iprime)*cos_tablex(j, jprime))

                    force_point = massx(iprime,jprime)/(sqrt(denom_point)*denom_point*rx_prime_squared(iprime))

                    force_r(i, j) = force_r(i, j) + force_point&
                            &*(rx_ratio(i, iprime)-cos_tablex(j, jprime))
                    force_theta(i, j) = force_theta(i, j) + force_point&
                            &*sin_tablex(j, jprime)
                end do
            end do
        end do
    end do

! fill the gaps in the force
    call fill_gaps(LEVEL, N_r0, N_theta0, force_r, force_theta)
    
! here the calculation has ended
    call CPU_TIME(finish)
    write(*,*) 'Calculation time =', finish-start, 'seconds'

    do i = 1, N_r0
        do j = 1, N_theta0
            write(11, '(e20.10)') force_r(i, j)
            write(12, '(e20.10)') force_theta(i, j)
        end do
    end do
    close(unit=11)
    close(unit=12)


    CONTAINS

        SUBROUTINE fill_radial_gaps(N_lvl, N_r0, N_theta0, f_r, f_theta)
            IMPLICIT NONE
            ! DECLARATIONS
            INTEGER, intent(in) :: N_lvl
            INTEGER, intent(in) :: N_r0, N_theta0
            REAL(8), DIMENSION(0:N_r0+2**N_lvl, 0:N_theta0+2**N_lvl), intent(inout) :: f_r, f_theta
            INTEGER :: i, j, ishift

            ! assign boundary conditions / ghost cells in force
            f_r(0, :) = 0                                               ! ghost cell
            f_theta(0, :) = 0                                           ! ghost cell
            f_r(:, 0) = f_r(:, N_theta0)                                ! periodicity
            f_theta(:, 0) = f_theta(:, N_theta0)                        ! periodicity
            f_r(N_r0+2**N_lvl, :) = 0                                   ! ghost cell
            f_theta(N_r0+2**N_lvl, :) = 0                               ! ghost cell
            f_r(:, N_theta0+2**N_lvl) = f_r(:, 2**N_lvl)                ! periodicity
            f_theta(:, N_theta0+2**N_lvl) = f_theta(:, 2**N_lvl)        ! periodicity
            
            do i = 1, N_r0, 2**N_lvl
                do j = 1, N_theta0, 2**N_lvl
                    do ishift = 1, 2**N_lvl-1
                        f_r(i+ishift, j) = f_r(i, j)*(1-ishift*1./2**N_lvl) + f_r(i+2**N_lvl, j)*ishift*1./2**N_lvl
                        f_theta(i+ishift, j) = f_theta(i, j)*(1-ishift*1./2**N_lvl) + f_theta(i+2**N_lvl, j)*ishift*1./2**N_lvl
                    end do
                end do
            end do
        END SUBROUTINE fill_radial_gaps

        SUBROUTINE fill_angular_gaps(N_lvl, N_r0, N_theta0, f_r, f_theta)
            IMPLICIT NONE
            ! DECLARATIONS
            INTEGER, intent(in) :: N_lvl
            INTEGER, intent(in) :: N_r0, N_theta0
            REAL(8), DIMENSION(0:N_r0+2**N_lvl, 0:N_theta0+2**N_lvl), intent(inout) :: f_r, f_theta
            INTEGER :: i, j, jshift
            
            do i = 1, N_r0, 2**N_lvl
                do j = 1, N_theta0, 2**N_lvl
                    do jshift = 1, 2**N_lvl-1
                        f_r(i, j+jshift) = f_r(i, j)*(1-jshift*1./2**N_lvl) + f_r(i, j+2**N_lvl)*jshift*1./2**N_lvl
                        f_theta(i, j+jshift) = f_theta(i, j)*(1-jshift*1./2**N_lvl) + f_theta(i, j+2**N_lvl)*jshift*1./2**N_lvl
                    end do
                end do
            end do
        END SUBROUTINE fill_angular_gaps

        SUBROUTINE fill_gaps(N_lvl, N_r0, N_theta0, f_r, f_theta)
            IMPLICIT NONE
            ! DECLARATIONS
            INTEGER, intent(in) :: N_lvl
            INTEGER, intent(in) :: N_r0, N_theta0
            REAL(8), DIMENSION(0:N_r0+2**N_lvl, 0:N_theta0+2**N_lvl), intent(inout) :: f_r, f_theta
            INTEGER :: i, j, ishift, jshift, factor_lvl
            REAL(8) :: inv_factor_2
            
            factor_lvl = 2**N_lvl
            inv_factor_2 = 1./(factor_lvl*factor_lvl)

            ! assign boundary conditions / ghost cells in force
            f_r(0, :) = 0                                             ! zero overdensity
            f_theta(0, :) = 0                                         ! zero overdensity
            f_r(:, 0) = f_r(:, N_theta0)                              ! periodicity
            f_theta(:, 0) = f_theta(:, N_theta0)                      ! periodicity
            f_r(N_r0+factor_lvl, :) = 0                               ! zero overdensity
            f_theta(N_r0+factor_lvl, :) = 0                           ! zero overdensity
            f_r(:, N_theta0+factor_lvl) = f_r(:, factor_lvl)          ! periodicity
            f_theta(:, N_theta0+factor_lvl) = f_theta(:, factor_lvl)  ! periodicity

            do i = 1, N_r0, factor_lvl
                do j = 1, N_theta0, factor_lvl
                    do ishift = 0, factor_lvl-1
                        do jshift = 0, factor_lvl-1
                            f_r(i+ishift, j+jshift) = ( (factor_lvl-ishift) * (factor_lvl-jshift) * f_r(i, j) + &
                                                       & ishift  * (factor_lvl-jshift) * f_r(i+factor_lvl, j) + &
                                                       &(factor_lvl-ishift) * jshift * f_r(i, j+factor_lvl)   + &
                                                       &       ishift * jshift * f_r(i+factor_lvl, j+factor_lvl)&
                                                       &)*inv_factor_2
                            f_theta(i+ishift, j+jshift) = ( (factor_lvl-ishift) * (factor_lvl-jshift) * f_theta(i, j) + &
                                                           & ishift  * (factor_lvl-jshift) * f_theta(i+factor_lvl, j) + &
                                                           &(factor_lvl-ishift) * jshift * f_theta(i, j+factor_lvl)   + &
                                                           &       ishift * jshift * f_theta(i+factor_lvl, j+factor_lvl)&
                                                           &)*inv_factor_2
                        end do
                    end do
                end do
            end do
        END SUBROUTINE fill_gaps

END PROGRAM oscillation_force
