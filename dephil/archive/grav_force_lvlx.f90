!--------------------------------------------------------------------------------------
! PROGRAM GRAV_FORCE_LVLX
!	compile: with make (Makefile included), e.g. make grav_force_lvlx
!       
!	use:	 ./grav_force_lvlx [LEVEL], e.g.  ./grav_force_lvlx 1
!--------------------------------------------------------------------------------------
PROGRAM grav_force_lvlx
    USE grav_parameters
    USE grav_IO
    USE grav_precalcs
    USE grav_mass
    IMPLICIT NONE

! DECLARATIONS
    INTEGER :: i, iprime, j, jprime
    INTEGER :: N_rx, N_thetax
    INTEGER :: LEVEL
    REAL(8)  :: epsilon2 ! this value seems to be the closest to the LEVEL0 solution
                         ! the value for LEVEL 2 seems to be close to 0.000336
    REAL(8) :: force_point, denom_point,&
             & f_rx=0., f_thetax=0.
    REAL :: start, finish

    CHARACTER(len=1) :: arg
!--------------------------------------------------------------------------------------
! READ THE LEVEL (FROM THE ARGUMENT LINE)

    call getarg(1, arg)
    read (arg, '(i5)') LEVEL

    write(*,*)"Calculating level "//arg//"..."

    N_rx = N_r0/2**LEVEL
    N_thetax = N_theta0/2**LEVEL

    if (LEVEL==0) then
        epsilon2 = 0
    else if (LEVEL==1) then
        epsilon2 = .000118
    else if (LEVEL==2) then
        epsilon2 = .000346
    else if (LEVEL==3) then
        epsilon2 = .000902      ! try different values
    else if (LEVEL==4) then
        epsilon2 = .001804
    else
        epsilon2 = .000001      ! try different values
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

    call precalcs_lvlx(LEVEL, N_r0, N_theta0, r0, theta0, dr0, dtheta0, rx, thetax,&
                                 & drx, dthetax, rx_squared, rx_prime_squared,&
                                 & rx_ratio, rx_ratio_squared,&
                                 & cos_tablex, sin_tablex)
    deallocate(r0_squared, r0_prime_squared, rx_squared)

!--------------------------------------------------------------------------------------
! PROPAGATION starts here
    call CPU_TIME(start)

! fill the <mass table> for LEVEL 0
    call calc_masslvl0(N_r0, N_theta0, sigma0, r0, dr0, dtheta0, mass0)

! fill the <mass table> for LEVEL 1
    call calc_masslvlx(LEVEL, N_r0, N_theta0, mass0, massx)

! write force components for every corner in grid
    do i = 1, N_r0
        do j = 1, N_theta0
            ! sum up the forces on the point (i, j)
            do iprime = 1, N_rx
                do jprime = 1, N_thetax
                    ! formula for the gravitational force split into four parts for faster calculation expressed in terms of ratios of r
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
            write(11, '(e20.10)') f_rx
            write(12, '(e20.10)') f_thetax
            f_rx = 0.
            f_thetax = 0.
        end do
    end do

! here the calculation has ended
    call CPU_TIME(finish)
    write(*,*) 'Calculation time =', finish-start, 'seconds'

    close(unit=11)
    close(unit=12)

END PROGRAM grav_force_lvlx
