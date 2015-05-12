!--------------------------------------------------------------------------------------
! PROGRAM OSCILLATION_MASS
!	compile: with make (Makefile included), e.g. make grav_pure_lvlx
!       
!	use:	 ./grav_pure_lvlx [LEVEL], e.g.  ./grav_pure_lvlx 2
!--------------------------------------------------------------------------------------
PROGRAM grav_pure_lvl
    USE grav_parameters
    USE grav_IO
    USE grav_precalcs
    USE grav_mass
    IMPLICIT NONE

! DECLARATIONS
    INTEGER :: i, iprime, j, jprime, ishift, jshift!, isign, jsign
    INTEGER :: N_rx, N_thetax
    INTEGER :: LEVEL, factor_x, half_factor_x

    REAL(8) :: inv_factor2, force_point, denom_point, vmass, dr_shift, dtheta_shift,&
             & ratio, cosine, r_corner, theta_corner, rprime, delta_theta,&
             & c1, c2, c3, c4!, epsilon2=1e-6

    REAL :: start, finish
    CHARACTER(len=1) :: arg
!--------------------------------------------------------------------------------------
! READ THE LEVEL (FROM THE ARGUMENT LINE)
    call getarg(1, arg)
    read (arg, '(i5)') LEVEL

    write(*,*)"Calculating level "//arg//"..."

    factor_x = 2**LEVEL
    half_factor_x = factor_x/2
    inv_factor2 = 1./(factor_x*factor_x)
    N_rx = N_r0/factor_x
    N_thetax = N_theta0/factor_x

    allocate(force_r(0:N_r0+factor_x, 0:N_theta0+factor_x), force_theta(0:N_r0+factor_x, 0:N_theta0+factor_x))     ! allocate with ghost cells

!--------------------------------------------------------------------------------------
! OPEN FILES
    open(unit=11, file=TRIM('./data/radial_pure_lvl'//arg//'.data'), action='write')
    open(unit=12, file=TRIM('./data/angular_pure_lvl'//arg//'.data'), action='write')

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
    deallocate(r0_squared, r0_prime_squared, rx_squared, rx_prime_squared)
    deallocate(sigmax)
!--------------------------------------------------------------------------------------
! PROPAGATION starts here
    call CPU_TIME(start)

! fill the <mass table> for LEVEL 0
    call calc_masslvl0(N_r0, N_theta0, sigma0, r0, dr0, dtheta0, mass0)

! fill the <mass table> for LEVEL 1
    call calc_masslvlx(LEVEL, N_r0, N_theta0, mass0, massx)


!! SHIFT WITH MOD METHOD
!! go through the force grid
    do i = 1, N_r0
        ishift = MOD(i-1, factor_x)                     ! shift index from nearest corner to the left
        dr_shift = r0(i)-r0(i-ishift)                   ! shift in r
        r_corner = r0(i)-.5*dr0(i)                      ! shift to the corner
        do j = 1, N_theta0
            jshift = MOD(j-1, factor_x)                 ! shift index from the nearest corner to the bottom
            dtheta_shift = theta0(j)-theta0(j-jshift)   ! shift in theta
            theta_corner = theta0(j)-.5*dtheta0(j)
            ! factors for interpolation of the shifted mass
            c1 = (factor_x-ishift)*(factor_x-jshift)
            c2 = ishift*(factor_x-jshift)
            c3 = (factor_x-ishift)*jshift
            c4 = ishift*jshift
            do iprime = 0, N_rx
                do jprime = 1, N_thetax
                    ! bilinear interpolation
                    vmass = ( c1 * massx(iprime, jprime)    + &
                             &c2 * massx(iprime+1, jprime)  + &
                             &c3 * massx(iprime, jprime+1)  + &
                             &c4 * massx(iprime+1, jprime+1)) &
                             & * inv_factor2
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


!! EXPANDED INTERPOLATION
!! go through the force grid
!    do i = 1, N_r0
!        ishift = MOD(i-1, factor_x)                     ! shift index from nearest corner to the left
!        dr_shift = r0(i)-r0(i-ishift)                   ! shift in r
!        r_corner = r0(i)-.5*dr0(i)                      ! shift to the corner
!        do j = 1, N_theta0
!            jshift = MOD(j-1, factor_x)                 ! shift index from the nearest corner to the bottom
!            dtheta_shift = theta0(j)-theta0(j-jshift)   ! shift in theta
!            theta_corner = theta0(j)-.5*dtheta0(j)
!            ! factors for interpolation of the shifted mass
!            c1 = (factor_x-ishift)*(factor_x-jshift)
!            c2 = ishift*(factor_x-jshift)/2
!            c3 = (factor_x-ishift)*jshift/2
!            c4 = ishift*jshift/2
!            do iprime = 1, N_rx
!                do jprime = 1, N_thetax
!                    ! bilinear interpolation
!                    vmass = ( c1 * massx(iprime, jprime)    + &
!                             &c2 * massx(iprime+1, jprime)  + &
!                             &c3 * massx(iprime, jprime+1)  + &
!                             &c4 * massx(iprime+1, jprime+1)+ &
!                             &c2 * massx(iprime-1, jprime)  + &
!                             &c3 * massx(iprime, jprime-1)  + &
!                             &c4 * massx(iprime-1, jprime-1)  &
!                             &)* 2 * inv_factor2
!                    rprime = rx(iprime)+dr_shift
!                    ratio = r_corner/rprime
!                    delta_theta = theta_corner-thetax(jprime)-dtheta_shift
!                    cosine = cos(delta_theta)
!                    denom_point = 1.+ratio*ratio-2.*ratio*cosine
!                    denom_point = sqrt(denom_point)*denom_point*rprime*rprime
!                    force_point = vmass/denom_point
!                    force_r(i,j) = force_r(i,j)+force_point*(ratio-cosine)
!                    force_theta(i,j) = force_theta(i,j)+force_point*sin(delta_theta)
!                end do
!            end do
!        end do
!    end do


!! NEGATIVE SHIFT INCLUDED
!! go through the force grid
!    do i = 1, N_r0
!        ishift = MOD(i-1, factor_x)                     ! shift index from nearest corner to the left
!        isign = 1
!        if (ishift > half_factor_x) then
!            ishift = half_factor_x-ishift
!            isign = -1
!        end if
!        dr_shift = r0(i)-r0(i-ishift)                   ! shift in r
!        r_corner = r0(i)-.5*dr0(i)                      ! shift to the corner
!        do j = 1, N_theta0
!            jshift = MOD(j-1, factor_x)                 ! shift index from the nearest corner to the bottom
!            jsign = 1
!            if (jshift > half_factor_x) then
!                jshift = half_factor_x-jshift
!                jsign = -1
!            end if
!            dtheta_shift = theta0(j)-theta0(j-jshift)   ! shift in theta
!            theta_corner = theta0(j)-.5*dtheta0(j)
!            ! factors for interpolation of the shifted mass
!            c1 = (half_factor_x-abs(ishift))*(half_factor_x-abs(jshift))
!            c2 = abs(ishift)*(half_factor_x-abs(jshift))
!            c3 = (half_factor_x-abs(ishift))*abs(jshift)
!            c4 = abs(ishift)*abs(jshift)
!            do iprime = 1, N_rx
!                do jprime = 1, N_thetax
!                    ! bilinear interpolation
!                    vmass = ( c1 * massx(iprime, jprime)    + &
!                             &c2 * massx(iprime+isign, jprime)  + &
!                             &c3 * massx(iprime, jprime+jsign)  + &
!                             &c4 * massx(iprime+isign, jprime+jsign)) &
!                             & * inv_factor2
!                    rprime = rx(iprime)+dr_shift
!                    ratio = r_corner/rprime
!                    delta_theta = theta_corner-thetax(jprime)-dtheta_shift
!                    cosine = cos(delta_theta)
!                    denom_point = 1.+ratio*ratio-2.*ratio*cosine
!                    denom_point = sqrt(denom_point)*denom_point*rprime*rprime
!                    force_point = vmass/denom_point
!                    force_r(i,j) = force_r(i,j)+force_point*(ratio-cosine)
!                    force_theta(i,j) = force_theta(i,j)+force_point*sin(delta_theta)
!                end do
!            end do
!        end do
!    end do


!! SHIFT IMPLEMENTATION WITH LOOPS
!! write force components for every corner in grid
!    do i = 1, N_r0, factor_x
!        do j = 1, N_theta0, factor_x
!            ! sum up the forces at the point (i, j)
!            do ishift = 0, factor_x-1
!                shifted_i = i+ishift    ! index to the point shifted from the corner where the force is summed up
!                dr_shift = r0(shifted_i)-r0(i)   ! distance of shift in r
!                do jshift = 0, factor_x-1
!                    shifted_j = j+jshift    ! index to the point shifted from the corner where the force is summed up
!                    dtheta_shift = theta0(shifted_j)-theta0(j)   ! distance of shift in theta
!                    do iprime = 1, N_rx
!                        rprime = rx(iprime)+dr_shift   ! shift the position of the masses
!                        ratio = (r0(shifted_i)-.5*dr0(shifted_i))/rprime
!                        ratio_squared = ratio*ratio
!                        rprime = rprime*rprime    ! iprime squared for later...
!                        do jprime = 1, N_thetax
!                            ! formula for the gravitational force split into four parts for faster calculation expressed in terms of ratios of r
!                            ! F_grav = Sum ( Mass / (r_iprime³*sqrt(1+r_ratio_squared-2*r_ratio*cos))³/² ) * r_iprime | (r_ratio - cos), (sin)
!                            vmass = ( (factor_x-ishift) * (factor_x-jshift)  * massx(iprime, jprime)+ &
!                                     & ishift  * (factor_x-jshift)         * massx(iprime+1, jprime)+ &
!                                     & (factor_x-ishift) * jshift          * massx(iprime, jprime+1)+ &
!                                     &          ishift * jshift          * massx(iprime+1, jprime+1)  &
!                                     &)*inv_factor2
!                            ! variables depending on iprime, jprime have to be shifted according to dr_shift / dtheta_shift
!                            ! rx_ratio_squared(i,iprime) ; rx_ratio(i,iprime) ; cos_tablex(j,jprime) ; rx_prime_squared(iprime) ; sin_tablex(j,jprime)
!                            delta_theta = theta0(shifted_j)-.5*dtheta0(shifted_j)-thetax(jprime)-dtheta_shift   ! theta - theta_prime
!                            cosine = cos(delta_theta)
!                            sine = sin(delta_theta)
!                            denom_point = 1.+ratio_squared-2.*ratio*cosine    ! denominator
!                            denom_point = sqrt(denom_point)*denom_point*rprime
!                            force_point = vmass/denom_point
!                            force_r(shifted_i, shifted_j) = force_r(shifted_i, shifted_j) + force_point&
!                                            &*(ratio-cosine)
!                            force_theta(shifted_i, shifted_j) = force_theta(shifted_i, shifted_j) + force_point&
!                                                &*sine
!                        end do
!                    end do
!                end do
!            end do
!        end do
!    end do


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

END PROGRAM grav_pure_lvl
