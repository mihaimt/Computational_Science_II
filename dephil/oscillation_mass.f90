!--------------------------------------------------------------------------------------
! PROGRAM OSCILLATION_MASS
!	compile: with make (Makefile included), e.g. make oscillation_mass
!       
!	use:	 ./oscillation_mass [LEVEL], e.g.  ./oscillation_mass 2
!--------------------------------------------------------------------------------------
PROGRAM oscillation_mass
    USE grav_parameters
    USE grav_IO
    USE grav_precalcs
    USE grav_mass
    IMPLICIT NONE

! DECLARATIONS
    INTEGER :: i, iprime, ishift, j, jprime, jshift, shifted_i, shifted_j
    INTEGER :: N_rx, N_thetax
    INTEGER :: LEVEL, factor_x

    REAL(8) :: inv_factor2, vmass, dr_shift, dtheta_shift, force_point, denom_point,&
               &ratio, ratio_squared, cosine, sine, rprime, delta_theta

    REAL :: start, finish

    CHARACTER(len=1) :: arg
!--------------------------------------------------------------------------------------
! READ THE LEVEL (FROM THE ARGUMENT LINE)

    call getarg(1, arg)
    read (arg, '(i5)') LEVEL

    factor_x = 2**LEVEL
    inv_factor2 = 1./(factor_x*factor_x)
    N_rx = N_r0/factor_x
    N_thetax = N_theta0/factor_x

    allocate(force_r(0:N_r0+factor_x, 0:N_theta0+factor_x), force_theta(0:N_r0+factor_x, 0:N_theta0+factor_x))     ! allocate with ghost cells

!--------------------------------------------------------------------------------------
! OPEN FILES
    open(unit=11, file=TRIM('./data/radial_osc_mass_lvl'//arg//'.data'), action='write')
    open(unit=12, file=TRIM('./data/angular_osc_mass_lvl'//arg//'.data'), action='write')

!--------------------------------------------------------------------------------------
! READ INPUT FILES
    call read_input(N_r0, N_theta0, r0, theta0, sigma0)   ! r0 and theta0 indicate the centers of the grid cells

! READ IN LEVEL ONE ------------------
    call grid_param_lvlx(LEVEL, N_r0, N_theta0, r0, theta0, sigma0, rx, thetax, sigmax)   ! rx and thetax indicate the centers of the grid cells in LEVEL X

!--------------------------------------------------------------------------------------
! PRECALCULATIONS
    call precalcs_lvl0(N_r0, N_theta0, r0, theta0, dr0, dtheta0, r0_squared, r0_prime_squared,&   ! r0_squared indicates the corners, r0_prime_squared indicates the centers
                       & r0_ratio, r0_ratio_squared, cos_table0, sin_table0)                      ! cos and sin with differences of corners to centers
    call precalcs_lvlx(LEVEL, N_r0, N_theta0, r0, theta0, dr0, dtheta0, rx, thetax,&
                                 & drx, dthetax, rx_squared, rx_prime_squared,&
                                 & rx_ratio, rx_ratio_squared,&
                                 & cos_tablex, sin_tablex)                                        ! same as above for LEVEL X, ratios indicate the corners to the centers

!--------------------------------------------------------------------------------------
! PROPAGATION starts here
    call CPU_TIME(start)

! fill the <mass table> for LEVEL 0
    call calc_masslvl0(N_r0, N_theta0, sigma0, r0, dr0, dtheta0, mass0)

! fill the <mass table> for LEVEL 1
    call calc_masslvlx(LEVEL, N_r0, N_theta0, mass0, massx)

! write force components for every corner in grid
    do i = 1, N_r0, factor_x
        do j = 1, N_theta0, factor_x
            ! sum up the forces at the point (i, j)
            do ishift = 0, factor_x-1
                shifted_i = i+ishift    ! index to the point shifted from the corner where the force is summed up
                dr_shift = r0(shifted_i)-r0(i)   ! distance of shift in r
                do jshift = 0, factor_x-1
                    shifted_j = j+jshift    ! index to the point shifted from the corner where the force is summed up
                    dtheta_shift = theta0(shifted_j)-theta0(j)   ! distance of shift in theta
                    do iprime = 1, N_rx
                        rprime = rx(iprime)+dr_shift   ! shift the position of the masses
                        ratio = (r0(shifted_i)-.5*dr0(shifted_i))/rprime
                        ratio_squared = ratio*ratio
                        rprime = rprime*rprime    ! iprime squared for later...
                        do jprime = 1, N_thetax
                            ! formula for the gravitational force split into four parts for faster calculation expressed in terms of ratios of r
                            ! F_grav = Sum ( Mass / (r_iprime³*sqrt(1+r_ratio_squared-2*r_ratio*cos))³/² ) * r_iprime | (r_ratio - cos), (sin)
                            vmass = ( (factor_x-ishift) * (factor_x-jshift)  * massx(iprime, jprime)+ &
                                     & ishift  * (factor_x-jshift)         * massx(iprime+1, jprime)+ &
                                     & (factor_x-ishift) * jshift          * massx(iprime, jprime+1)+ &
                                     &          ishift * jshift          * massx(iprime+1, jprime+1)  &
                                     &)*inv_factor2
                            ! variables depending on iprime, jprime have to be shifted according to dr_shift / dtheta_shift
                            ! rx_ratio_squared(i,iprime) ; rx_ratio(i,iprime) ; cos_tablex(j,jprime) ; rx_prime_squared(iprime) ; sin_tablex(j,jprime)
                            delta_theta = theta0(shifted_j)-.5*dtheta0(shifted_j)-thetax(jprime)-dtheta_shift   ! theta - theta_prime
                            cosine = cos(delta_theta)
                            sine = sin(delta_theta)
                            denom_point = 1.+ratio_squared-2.*ratio*cosine    ! denominator
                            denom_point = sqrt(denom_point)*denom_point*rprime
                            force_point = vmass/denom_point
                            force_r(shifted_i, shifted_j) = force_r(shifted_i, shifted_j) + force_point&
                                            &*(ratio-cosine)
                            force_theta(shifted_i, shifted_j) = force_theta(shifted_i, shifted_j) + force_point&
                                                &*sine
                        end do
                    end do
                end do
            end do
        end do
    end do
    
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

END PROGRAM oscillation_mass
