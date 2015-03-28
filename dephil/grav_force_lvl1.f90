PROGRAM grav_force_lvl1
    IMPLICIT NONE

! DECLARATIONS
    INTEGER :: i, iprime, j, jprime
    INTEGER, PARAMETER :: N_r0=128, N_theta0=256,&
                        & N_r1=N_r0/2, N_theta1=N_theta0/2

    REAL(8) :: r0_i, theta0_j, density,&
             & r1_i, theta1_j, diff_theta,&
             & r1_iprime,&
             & r1_i_squared, r1_iprime_squared,&
             & r1sum_squared,&
             & r1_iprime_cos,&
             & force_point, denom_point,&
             & f_r1=0., f_theta1=0.

    REAL(8), DIMENSION(N_r0) :: r0
    REAL(8), DIMENSION(N_r1) :: r1, dr1
    REAL(8), DIMENSION(N_theta0) :: theta0
    REAL(8), DIMENSION(N_theta1) :: theta1, dtheta1
    REAL(8), DIMENSION(N_r0, N_theta0) :: sigma0
    REAL(8), DIMENSION(N_r1, N_theta1) :: sigma1, mass1
    REAL(8), DIMENSION(1-N_theta1:N_theta1-1) :: cos_table1, sin_table1

    REAL :: start, finish
!--------------------------------------------------------------------------------------
! READ IN FILES
    open(unit=8, file='./data/r_project.data', status='old', action='read')
    open(unit=9, file='./data/theta_project.data', status='old', action='read')
    open(unit=10, file='./data/density_project.data', status='old', action='read')
    open(unit=11, file='./data/f_radial_lvl1.data', status='old', action='write')
    open(unit=12, file='./data/f_angular_lvl1.data', status='old', action='write')

! read in LEVEL ZERO ------------------
! read radii into 1-D array for lvl0
    do i = 1, N_r0
        read(8, '(e20.10)') r0_i
        r0(i) = r0_i
    end do
    close(unit=8)
! read angles into 1-D array for lvl0
    do j = 1, N_theta0
        read(9, '(e20.10)') theta0_j
        theta0(j) = theta0_j
    end do
    close(unit=9)
! read density into 2-D array for lvl0
    do i = 2, N_r0
        do j = 1, N_theta0
            read(10, '(e20.10)') density
            sigma0(i,j) = density
        end do
    end do
    close(unit=10)

! read in LEVEL ONE ------------------
! read radii into 1-D array for lvl1
    do i = 2, N_r0, 2
        r1(i/2) = (r0(i)+r0(i-1))*.5
    end do
! read angles into 1-D array for lvl1
    do j = 2, N_theta0, 2
        theta1(j/2) = (theta0(j)+theta0(j-1))
    end do
! read density into 2-D array for lvl1
    do i = 2, N_r0, 2
        do j = 2, N_theta0, 2
            sigma1(i/2, j/2) = sigma0(i, j)+sigma0(i-1, j)+sigma0(i, j-1)+sigma0(i-1, j-1)
        end do
    end do

!--------------------------------------------------------------------------------------
! CALCULATION OF THE FORCE COMPONENTS

! PRECALCULATIONS

! different dr/dtheta for LEVEL 1 ------------------
    do i = 1, N_r1-1
        dr1(i) = r1(i+1)-r1(i)
    end do
    dr1(N_r1) = dr1(N_r1-1)
        do j = 1, N_theta1-1
        dtheta1(j) = theta1(j+1)-theta1(j)
    end do
    dtheta1(N_theta1) = dtheta1(N_theta1-1)

! fill the sin/cos table for LEVEL 1 ------------------
    do j = 1, N_theta1
        theta1_j = theta1(j)-.5*dtheta1(j)
        do jprime = 1, N_theta1
            diff_theta = theta1_j-theta1(jprime)
            cos_table1(j-jprime) = cos(diff_theta)
            sin_table1(j-jprime) = sin(diff_theta)
        end do
    end do

! PROPAGATION starts here
    call CPU_TIME(start)

! fill the mass table for LEVEL 1
    do iprime = 1, N_r1
        r1_iprime = r1(iprime)*dr1(iprime) ! multiply dr already here to save some operations
        do jprime = 1,N_theta1
            mass1(iprime,jprime) = -sigma1(iprime,jprime)*r1_iprime*dtheta1(jprime)
        end do
    end do


! write force components for every corner in grid ------------------------------------
    do i = 1, N_r1
        r1_i = r1(i)-.5*dr1(i) ! shift r to the corners
        r1_i_squared = r1_i*r1_i   ! better than calculating it in formula
        do j = 1, N_theta1
            theta1_j = theta1(j)-.5*dtheta1(j) ! shift theta to the corners
            ! sum up the forces on the point (i, j)
            do iprime = 1, N_r1
                r1_iprime = r1(iprime)
                r1_iprime_squared = r1_iprime*r1_iprime   ! slight performance improvement
                r1sum_squared = r1_i_squared+r1_iprime_squared   ! slight performance improvement
                do jprime = 1, N_theta1
                    ! formula for the gravitational force split into four parts for faster calculation
                    r1_iprime_cos = r1_iprime*cos_table1(j-jprime)

                    !denom_point1 = InvSqrt(rsum_squared-2.*r_i*r_iprime_cos)
                    denom_point = 1./sqrt(r1sum_squared-2.*r1_i*r1_iprime_cos)
                    !denom_point = c_invsqrt64(rsum_squared-2.*r_i*r_iprime_cos)   ! returns somehow Nan's

                    denom_point = denom_point*denom_point*denom_point
                    force_point = mass1(iprime,jprime)*denom_point

                    f_r1 = f_r1 + force_point&
                            &*(r1_i-r1_iprime_cos)
                    f_theta1 = f_theta1 + force_point&
                            &*r1_iprime*sin_table1(j-jprime)
                end do
            end do
        write(11, '(e20.10)') f_r1
        write(12, '(e20.10)') f_theta1
        f_r1 = 0.
        f_theta1 = 0.
        end do
    end do

! here the calculation has ended
    call CPU_TIME(finish)
    write(*,*) 'Calculation time =', finish-start, 'seconds'

    close(unit=11)
    close(unit=12)

END PROGRAM grav_force_lvl1
