!-------------------------------------------------------------------------
! PROGRAM grav_refined
!    compile: with make (Makefile included), e.g. make grav_refined
!
!    use: ./grav_refined 
!-------------------------------------------------------------------------
PROGRAM grav_refined
    USE grav_parameters
    USE grav_IO
    USE grav_precalcs
    USE grav_mass
    IMPLICIT NONE

! DECLARATIONS
    INTEGER :: i, iprime, j, jprime
    INTEGER :: ishift, jshift, ishiftref, jshiftref
    INTEGER :: iref_low, iref_up, ix_low, ix_up
    INTEGER :: jref_low, jref_up, jx_low, jx_up
    INTEGER :: periodic_low, periodic_up
    INTEGER :: N_rx, N_thetax, N_rxref, N_thetaxref
    INTEGER :: LEVEL, factor_x, half_factor_x
    INTEGER :: LEVELref, factor_xref, half_factor_xref

    REAL(8) :: inv_factor2, inv_factor2ref, force_point, denom_point,&
             & c1, c2, c3, c4, c1ref, c2ref, c3ref, c4ref,&
             & r_corner, theta_corner,&
             & dr_shift, dr_shiftref,&
             & dtheta_shift, dtheta_shiftref,&
             & vmass, rprime, ratio, delta_theta, cosine

    REAL :: start, finish
    CHARACTER(len=1) :: arg
!--------------------------------------------------------------------------------------
! READ THE LEVEL (FROM THE ARGUMENT LINE)
    call getarg(1, arg)
    read (arg, '(i5)') LEVEL

    write(*,*)"Setting up grid surface of dimensions ",N_r0,N_theta0
    write(*,*)"Calculating level "//arg//" with refinement..."

    factor_x = 2**LEVEL
    half_factor_x = factor_x/2
    inv_factor2 = 1./(factor_x*factor_x)
    N_rx = N_r0/factor_x
    N_thetax = N_theta0/factor_x

    LEVELref = LEVEL-1
    factor_xref = 2**LEVELref
    half_factor_xref = factor_xref/2
    inv_factor2ref = 1./(factor_xref*factor_xref)
    N_rxref = N_r0/factor_xref
    N_thetaxref = N_theta0/factor_xref

    allocate(force_r(0:N_r0+factor_x, 0:N_theta0+factor_x), force_theta(0:N_r0+factor_x, 0:N_theta0+factor_x))     ! allocate with ghost cells

!--------------------------------------------------------------------------------------
! OPEN FILES
    open(unit=11, file=TRIM('./data/radial_refined_lvl'//arg//'.data'), action='write')
    open(unit=12, file=TRIM('./data/angular_refined_lvl'//arg//'.data'), action='write')

!--------------------------------------------------------------------------------------
! READ INPUT FILES
    call read_input(N_r0, N_theta0, r0, theta0, sigma0)   ! r0 and theta0 indicate the centers of the grid cells

! READ IN LEVEL X ------------------
    call grid_param_lvlx(LEVEL, N_r0, N_theta0, r0, theta0, sigma0, rx, thetax, sigmax)   ! rx and thetax indicate the centers of the grid cells in LEVEL X
    call grid_param_lvlx(LEVELref, N_r0, N_theta0, r0, theta0, sigma0, rxref, thetaxref, sigmaxref)   ! rxref and thetaxref indicate the centers of the grid cells in LEVEL X-1

!--------------------------------------------------------------------------------------
! PRECALCULATIONS
    call precalcs_lvl0(N_r0, N_theta0, r0, theta0, dr0, dtheta0, r0_squared, r0_prime_squared,&   ! r0_squared indicates the corners, r0_prime_squared indicates the centers
                       & r0_ratio, r0_ratio_squared, cos_table0, sin_table0)                      ! cos and sin with differences of corners to centers
    call precalcs_lvlx(LEVEL, N_r0, N_theta0, r0, theta0, dr0, dtheta0, rx, thetax,&
                                 & drx, dthetax, rx_squared, rx_prime_squared,&
                                 & rx_ratio, rx_ratio_squared,&
                                 & cos_tablex, sin_tablex)                                        ! same as above for LEVEL X, ratios indicate the corners to the centers
    call precalcs_lvlx(LEVELref, N_r0, N_theta0, r0, theta0, dr0, dtheta0, rxref, thetaxref,&
                                 & drxref, dthetaxref, rxref_squared, rxref_prime_squared,&
                                 & rxref_ratio, rxref_ratio_squared,&
                                 & cos_tablexref, sin_tablexref)
    deallocate(r0_squared, r0_prime_squared, rx_squared, rx_prime_squared)
    deallocate(sigmax, sigmaxref)

!--------------------------------------------------------------------------------------
! PROPAGATION starts here
    call CPU_TIME(start)

! fill the <mass table> for LEVEL 0
    call calc_masslvl0(N_r0, N_theta0, sigma0, r0, dr0, dtheta0, mass0)

! fill the <mass table> for LEVEL X
    call calc_masslvlx(LEVEL, N_r0, N_theta0, mass0, massx)
! fill the <mass table> for LEVEL X-1
    call calc_masslvlx(LEVELref, N_r0, N_theta0, mass0, massxref)


!! LIMIT CONVERSION METHOD
!! go through the force grid
    do i = 1, N_r0
        r_corner = r0(i)-.5*dr0(i)
        ! index shifts for interpolation in r
        ishift = MOD(i-1, factor_x)
        ishiftref = MOD(ishift, factor_xref)
        ! shifts in r
        dr_shift = r0(i)-r0(i-ishift)
        dr_shiftref = r0(i)-r0(i-ishiftref)
        ! refined area indices in r
        iref_low = (i-factor_x-ishiftref-1)/factor_xref+1
        iref_up = (i+factor_x-ishiftref-2)/factor_xref+1
        ix_low = (i-factor_x-ishift-1)/factor_x+1
        ix_up = (i+factor_x-ishift-2)/factor_x+1
        do j = 1, N_theta0
            theta_corner = theta0(j)-.5*dtheta0(j)
            ! index shifts for interpolation in theta
            jshift = MOD(j-1, factor_x)
            jshiftref = MOD(jshift, factor_xref)
            ! shifts in theta
            dtheta_shift = (theta0(j)-theta0(j-jshift))
            dtheta_shiftref = theta0(j)-theta0(j-jshiftref)
            ! refined area indices in theta
            jref_low = (j-factor_x-jshiftref-1)/factor_xref+1
            jref_up = (j+factor_x-jshiftref-2)/factor_xref+1
            jx_low = (j-factor_x-jshift-1)/factor_x+1
            jx_up = (j+factor_x-jshift-2)/factor_x+1
            ! index shift for the last two loops if true
            periodic_low = 0
            if (jx_low < 1) then
                periodic_low = 1
            end if
            periodic_up = 0
            if (jx_up > N_thetax) then
                periodic_up = 1
            end if
            ! bilinear interpolation coefficients
            c1 = (factor_x-ishift)*(factor_x-jshift)
            c2 = ishift*(factor_x-jshift)
            c3 = (factor_x-ishift)*jshift
            c4 = ishift*jshift
            c1ref = (factor_xref-ishiftref)*(factor_xref-jshiftref)
            c2ref = ishiftref*(factor_xref-jshiftref)
            c3ref = (factor_xref-ishiftref)*jshiftref
            c4ref = ishiftref*jshiftref
            ! sum up refined area (LEVEL X-1)
            do iprime = iref_low, iref_up
                do jprime = jref_low, jref_up
                    vmass = ( c1ref * massxref(iprime, jprime)    +&
                             &c2ref * massxref(iprime+1, jprime)  +&
                             &c3ref * massxref(iprime, jprime+1)  +&
                             &c4ref * massxref(iprime+1, jprime+1))&
                             & * inv_factor2ref
                    ! calculate shifted values
                    rprime = rxref(iprime)+dr_shiftref
                    ratio = r_corner/rprime
                    delta_theta = theta_corner-thetaxref(jprime)-dtheta_shiftref
                    cosine = cos(delta_theta)
                    denom_point = 1.+ratio*ratio-2.*ratio*cosine
                    denom_point = sqrt(denom_point)*denom_point*rprime*rprime
                    force_point = vmass/denom_point
                    force_r(i,j) = force_r(i,j)+force_point*(ratio-cosine)
                    force_theta(i,j) = force_theta(i,j)+force_point*sin(delta_theta)
                end do
            end do
            ! sum up the rest as Level X in four parts/loops
            do iprime = 0, ix_low-1
                do jprime = 1, N_thetax
                    vmass = ( c1 * massx(iprime, jprime)    +&
                             &c2 * massx(iprime+1, jprime)  +&
                             &c3 * massx(iprime, jprime+1)  +&
                             &c4 * massx(iprime+1, jprime+1))&
                             & * inv_factor2
                    ! calculate shifted values
                    rprime = rx(iprime)+dr_shift
                    ratio = r_corner/rprime
                    delta_theta = theta_corner-thetax(jprime)-dtheta_shift
                    cosine = cos(delta_theta)
                    denom_point = 1.+ratio*ratio-2.*ratio*cosine
                    denom_point = sqrt(denom_point)*denom_point*rprime*rprime
                    force_point = vmass/denom_point
                    force_r(i,j) = force_r(i,j)+force_point*(ratio-cosine)
                    force_theta(i,j) = force_theta(i,j)+force_point*sin(delta_theta)
                end do
            end do
            do iprime = ix_up+1, N_rx
                do jprime = 1, N_thetax
                    vmass = ( c1 * massx(iprime, jprime)    +&
                             &c2 * massx(iprime+1, jprime)  +&
                             &c3 * massx(iprime, jprime+1)  +&
                             &c4 * massx(iprime+1, jprime+1))&
                             & * inv_factor2
                    ! calculate shifted values
                    rprime = rx(iprime)+dr_shift
                    ratio = r_corner/rprime
                    delta_theta = theta_corner-thetax(jprime)-dtheta_shift
                    cosine = cos(delta_theta)
                    denom_point = 1.+ratio*ratio-2.*ratio*cosine
                    denom_point = sqrt(denom_point)*denom_point*rprime*rprime
                    force_point = vmass/denom_point
                    force_r(i,j) = force_r(i,j)+force_point*(ratio-cosine)
                    force_theta(i,j) = force_theta(i,j)+force_point*sin(delta_theta)
                end do
            end do
            do iprime = ix_low, ix_up
                do jprime = 1+periodic_up, jx_low-1
                    vmass = ( c1 * massx(iprime, jprime)    +&
                             &c2 * massx(iprime+1, jprime)  +&
                             &c3 * massx(iprime, jprime+1)  +&
                             &c4 * massx(iprime+1, jprime+1))&
                             & * inv_factor2
                    ! calculate shifted values
                    rprime = rx(iprime)+dr_shift
                    ratio = r_corner/rprime
                    delta_theta = theta_corner-thetax(jprime)-dtheta_shift
                    cosine = cos(delta_theta)
                    denom_point = 1.+ratio*ratio-2.*ratio*cosine
                    denom_point = sqrt(denom_point)*denom_point*rprime*rprime
                    force_point = vmass/denom_point
                    force_r(i,j) = force_r(i,j)+force_point*(ratio-cosine)
                    force_theta(i,j) = force_theta(i,j)+force_point*sin(delta_theta)
                end do
            end do
            do iprime = ix_low, ix_up
                do jprime = jx_up+1, N_thetax-periodic_low
                    vmass = ( c1 * massx(iprime, jprime)    +&
                             &c2 * massx(iprime+1, jprime)  +&
                             &c3 * massx(iprime, jprime+1)  +&
                             &c4 * massx(iprime+1, jprime+1))&
                             & * inv_factor2
                    ! calculate shifted values
                    rprime = rx(iprime)+dr_shift
                    ratio = r_corner/rprime
                    delta_theta = theta_corner-thetax(jprime)-dtheta_shift
                    cosine = cos(delta_theta)
                    denom_point = 1.+ratio*ratio-2.*ratio*cosine
                    denom_point = sqrt(denom_point)*denom_point*rprime*rprime
                    force_point = vmass/denom_point
                    force_r(i,j) = force_r(i,j)+force_point*(ratio-cosine)
                    force_theta(i,j) = force_theta(i,j)+force_point*sin(delta_theta)
                end do
            end do
        end do
    end do

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

END PROGRAM grav_refined
