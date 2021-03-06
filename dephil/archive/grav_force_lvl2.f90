PROGRAM grav_force_lvl2
    IMPLICIT NONE

! DECLARATIONS
    INTEGER :: i, iprime, j, jprime
    INTEGER, PARAMETER :: N_r0=128, N_theta0=256,&
                        & N_r2=N_r0/4, N_theta2=N_theta0/4
    REAL(8), PARAMETER :: epsilon2=0.000346 ! this value seems to be the closest to the LEVEL0 solution

    REAL(8) :: density,&
             & r0_i, r2_i,&
             & theta0_j, theta2_j, diff_theta,&
             & r0_iprime, r2_iprime,&
             & force_point, denom_point,&
             & f_r2=0., f_theta2=0.

    REAL(8), DIMENSION(N_r0) :: r0, dr0, r0_squared, r0_prime_squared
    REAL(8), DIMENSION(N_r2) :: r2, dr2, r2_squared, r2_prime_squared
    REAL(8), DIMENSION(N_theta0) :: theta0, dtheta0
    REAL(8), DIMENSION(N_theta2) :: theta2, dtheta2
    REAL(8), DIMENSION(N_r0, N_r0) :: r0_ratio, r0_ratio_squared
    REAL(8), DIMENSION(N_r0, N_r2) :: r2_ratio, r2_ratio_squared
    REAL(8), DIMENSION(N_r0, N_theta0) :: sigma0, mass0
    REAL(8), DIMENSION(N_r2, N_theta2) :: sigma2, mass2
    REAL(8), DIMENSION(1-N_theta0:N_theta0-1) :: cos_table0, sin_table0
    REAL(8), DIMENSION(N_theta0, N_theta2) :: cos_table2, sin_table2

    REAL :: start, finish
!--------------------------------------------------------------------------------------
! READ IN FILES
    open(unit=8, file='./data/r_project.data', status='old', action='read')
    open(unit=9, file='./data/theta_project.data', status='old', action='read')
    open(unit=10, file='./data/density_project.data', status='old', action='read')
    open(unit=11, file='./data/f_radial_lvl2.data', status='old', action='write')
    open(unit=12, file='./data/f_angular_lvl2.data', status='old', action='write')

! read in LEVEL ZERO ------------------
! read <radii> into 1-D array for LEVEL0
    do i = 1, N_r0
        read(8, '(e20.10)') r0_i
        r0(i) = r0_i
    end do
    close(unit=8)
! read <angles> into 1-D array for LEVEL0
    do j = 1, N_theta0
        read(9, '(e20.10)') theta0_j
        theta0(j) = theta0_j
    end do
    close(unit=9)
! read <density> into 2-D array for LEVEL0
    do i = 1, N_r0
        do j = 1, N_theta0
            read(10, '(e20.10)') density
            sigma0(i,j) = density
        end do
    end do
    close(unit=10)

! read in LEVEL TWO ------------------
! read <radii> into 1-D array for LEVEL2
    do i = 4, N_r0, 4
        r2(i/4) = (r0(i)+r0(i-3))*.5
    end do
! read <angles> into 1-D array for LEVEL2
    do j = 4, N_theta0, 4
        theta2(j/4) = (theta0(j)+theta0(j-3))*.5
    end do
! read <density> into 2-D array for LEVEL2
    do i = 4, N_r0, 4
        do j = 4, N_theta0, 4
            sigma2(i/4, j/4) = sigma0(i, j)+sigma0(i-1, j)+sigma0(i, j-1)+sigma0(i-1, j-1)&
                            &+sigma0(i-2,j)+sigma0(i-2,j-1)+sigma0(i,j-2)+sigma0(i-1,j-2)+sigma0(i-2,j-2)&
                            &+sigma0(i-3,j)+sigma0(i-3,j-1)+sigma0(i-3,j-2)&
                            &+sigma0(i,j-3)+sigma0(i-1,j-3)+sigma0(i-2,j-3)+sigma0(i-3,j-3)
        end do
    end do

!--------------------------------------------------------------------------------------
! PRECALCULATIONS

! different <dr/dtheta> for LEVEL 0 ------------------
    do i = 1, N_r0-1
        dr0(i) = r0(i+1)-r0(i)
    end do
    dr0(N_r0) = dr0(N_r0-1)
    do j = 1, N_theta0-1
        dtheta0(j) = theta0(j+1)-theta0(j)
    end do
    dtheta0(N_theta0) = dtheta0(N_theta0-1)

! different <dr/dtheta> for LEVEL2 ------------------
    do i = 1, N_r2-1
        dr2(i) = r2(i+1)-r2(i)
    end do
    dr2(N_r2) = dr2(N_r2-1)
    do j = 1, N_theta2-1
        dtheta2(j) = theta2(j+1)-theta2(j)
    end do
    dtheta2(N_theta2) = dtheta2(N_theta2-1)

! <r²> in the corners and <r'²> in the centres for LEVEL0 ---------
    do i = 1, N_r0
        r0_i = r0(i)-.5*dr0(i) ! shift r to the corners
        r0_squared(i) = r0_i*r0_i
        r0_prime_squared(i) = r0(i)*r0(i)
    end do

! <r²> in the corners and <r'²> in the centres for LEVEL2 ---------
    do i = 1, N_r2
        r2_i = r2(i)-.5*dr2(i) ! shift r to the corners
        r2_squared(i) = r2_i*r2_i
        r2_prime_squared(i) = r2(i)*r2(i)
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

! <r ratios> and <r ratios squared> for LEVEL2 ------------------
    do i = 1, N_r0
        r0_i = r0(i)-.5*dr0(i) ! shift r to the corners
        do iprime = 1, N_r2
            r2_iprime = r0_i/r2(iprime) ! r_iprime only used as temp
            r2_ratio(i, iprime) = r2_iprime
            r2_ratio_squared(i, iprime) = r2_iprime*r2_iprime
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

! fill the <sin/cos table> for LEVEL2 ------------------
    do j = 1, N_theta0
        theta0_j = theta0(j)-.5*dtheta0(j)
        do jprime = 1, N_theta2
            diff_theta = theta0_j-theta2(jprime)
            cos_table2(j, jprime) = cos(diff_theta)
            sin_table2(j, jprime) = sin(diff_theta)
        end do
    end do

!--------------------------------------------------------------------------------------
! PROPAGATION starts here
    call CPU_TIME(start)

! fill the <mass table> for LEVEL 0
    do iprime = 1, N_r0
        r0_iprime = r0(iprime)*dr0(iprime) ! multiply dr already here to save some operations
        do jprime = 1, N_theta0
            mass0(iprime,jprime) = -sigma0(iprime,jprime)*r0_iprime*dtheta0(jprime)
        end do
    end do

! fill the <mass table> for LEVEL2
    do i = 4, N_r0, 4
        do j = 4, N_theta0, 4
            mass2(i/4, j/4) = mass0(i, j)+mass0(i-1, j)+mass0(i, j-1)+mass0(i-1, j-1)&
                            &+mass0(i-2,j)+mass0(i-2,j-1)+mass0(i,j-2)+mass0(i-1,j-2)+mass0(i-2,j-2)&
                            &+mass0(i-3,j)+mass0(i-3,j-1)+mass0(i-3,j-2)+mass0(i,j-3)+mass0(i-1,j-3)+mass0(i-2,j-3)+mass0(i-3,j-3)

        end do
    end do

    write(*,*)"Calculating level 2..."

! write force components for every corner in grid
    do i = 1, N_r0
        do j = 1, N_theta0
            ! sum up the forces on the point (i, j)
            do iprime = 1, N_r2
                do jprime = 1, N_theta2
                    ! formula for the gravitational force split into four parts for faster calculation
                            ! expressed in terms of ratios of r
                            ! F_grav = Sum ( Mass / (r_iprime*sqrt(1+r_ratio_squared-2*r_ratio*cos))³/² ) * r_iprime | (r_ratio - cos), (sin)

                    !denom_point = InvSqrt(1+r2_ratio_squared(i, iprime)-2*r2_ratio(i, iprime)*cos_table2(j, jprime)+epsilon2)
                    denom_point = sqrt(1+r2_ratio_squared(i, iprime)-2*r2_ratio(i, iprime)*cos_table2(j, jprime)+epsilon2)
                    !denom_point = c_invsqrt64(1+r2_ratio_squared(i, iprime)-2*r2_ratio(i, iprime)*cos_table2(j, jprime)+epsilon2)   ! returns somehow Nan's

                    force_point = mass2(iprime,jprime)/(denom_point*denom_point*denom_point*r2_prime_squared(iprime))

                    f_r2 = f_r2 + force_point&
                            &*(r2_ratio(i, iprime)-cos_table2(j, jprime))
                    f_theta2 = f_theta2 + force_point&
                            &*sin_table2(j, jprime)
                end do
            end do
        write(11, '(e20.10)') f_r2
        write(12, '(e20.10)') f_theta2
        f_r2 = 0.
        f_theta2 = 0.
        end do
    end do

! here the calculation has ended
    call CPU_TIME(finish)
    write(*,*) 'Calculation time =', finish-start, 'seconds'

    close(unit=11)
    close(unit=12)

! difference of the forces
!    open(unit=13, file='./data/f_radial_diff0to2.data', status='old', action='write')
!    open(unit=14, file='./data/f_radial.data', status='old', action='read')
!    open(unit=15, file='./data/f_radial_lvl2.data', status='old', action='read')

!    do i = 1, N_r0*N_theta0
!        read(14, '(e20.10)') f_r
!        read(15, '(e20.10)') f_r2
!        write(13, '(e20.10)') f_r-f_r2
!    end do

!    close(unit=13)
!    close(unit=14)
!    close(unit=15)


END PROGRAM grav_force_lvl2
