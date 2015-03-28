! Program to compute the force of a density map
!
!	to compile use make (Makefile included)
!
PROGRAM grav_force
    IMPLICIT NONE

! DECLARATIONS
    INTEGER :: i, iprime, j, jprime
    INTEGER, PARAMETER :: N_r=128, N_theta=256
    REAL :: start, finish
    REAL(8) :: density,&
             & r_i, r_i_squared,&
             & theta_j, diff_theta,&
             & r_iprime, r_iprime_squared,&
             & force_point, denom_point,&
             & f_r=0., f_theta=0.,&
             & r_iprime_cos, rsum_squared

    REAL(8), DIMENSION(N_r) :: r, dr
    REAL(8), DIMENSION(N_theta) :: theta, dtheta
    REAL(8), DIMENSION(N_r, N_theta) :: sigma, mass
    REAL(8), DIMENSION(1-N_theta:N_theta-1) :: cos_table, sin_table
    !REAL(8) :: InvSqrt
    !REAL(8), EXTERNAL :: c_invsqrt64

!--------------------------------------------------------------------------------------
! READ IN FILES
! Attention: define the path to the file
    open(unit=8, file='/home/dephil/grav_project/data/r_project.data', status='old', action='read')
    open(unit=9, file='./data/theta_project.data', status='old', action='read')
    open(unit=10, file='./data/density_project.data', status='old', action='read')
    open(unit=11, file='./data/f_radial.data', status='old', action='write')
    open(unit=12, file='./data/f_angular.data', status='old', action='write')

! read radii into 1-D array
    do i = 1, N_r
        read(8, '(e20.10)') r_i
        r(i) = r_i
    end do
    close(unit=8)
! read angles into 1-D array
    do j = 1, N_theta
        read(9, '(e20.10)') theta_j
        theta(j) = theta_j
    end do
    close(unit=9)
! read density into 2-D array
    do i = 1, N_r
        do j = 1, N_theta
            read(10, '(e20.10)') density
            sigma(i,j) = density
        end do
    end do
    close(unit=10)

!--------------------------------------------------------------------------------------
! CALCULATION OF THE FORCE COMPONENTS

! PRECALCULATIONS

! an array for all the different dr/dtheta ! now constant, but later useful
    do i = 1, N_r-1
        dr(i) = r(i+1)-r(i)
    end do
    dr(N_r) = dr(N_r-1)
    do j = 1, N_theta-1
        dtheta(j) = theta(j+1)-theta(j)
    end do
    dtheta(N_theta) = dtheta(N_theta-1)

! fill the sin/cos table
    do j = 1,N_theta
        theta_j = theta(j)-.5*dtheta(j)
        do jprime = 1,N_theta
            diff_theta = theta_j-theta(jprime)
            cos_table(j-jprime) = cos(diff_theta)
            sin_table(j-jprime) = sin(diff_theta)
        end do
    end do

! PROPAGATION starts here
    call CPU_TIME(start)
! fill the mass table
    do iprime = 1, N_r
        r_iprime = r(iprime)*dr(iprime) ! multiply dr already here to save some operations
        do jprime = 1,N_theta
            mass(iprime,jprime) = -sigma(iprime,jprime)*r_iprime*dtheta(jprime)
        end do
    end do

! write force components for every corner in grid
    do i = 1, N_r
        r_i = r(i)-.5*dr(i) ! shift r to the corners
        r_i_squared = r_i*r_i   ! better than calculating it in formula
        do j = 1, N_theta
            theta_j = theta(j)-.5*dtheta(j) ! shift theta to the corners
            ! sum up the forces on the point (i, j)
            do iprime = 1, N_r
                r_iprime = r(iprime)
                r_iprime_squared = r_iprime*r_iprime   ! slight performance improvement
                rsum_squared = r_i_squared+r_iprime_squared   ! slight performance improvement
                do jprime = 1, N_theta
                    ! formula for the gravitational force split into four parts for faster calculation
                    r_iprime_cos = r_iprime*cos_table(j-jprime)
                    
                    !denom_point = InvSqrt(rsum_squared-2.*r_i*r_iprime_cos)
                    denom_point = 1./sqrt(rsum_squared-2.*r_i*r_iprime_cos)
                    !denom_point = c_invsqrt64(rsum_squared-2.*r_i*r_iprime_cos)   ! returns somehow Nan's
                    
                    denom_point = denom_point*denom_point*denom_point
                    force_point = mass(iprime,jprime)*denom_point

                    f_r = f_r + force_point&
                            &*(r_i-r_iprime_cos)
                    f_theta = f_theta + force_point&
                            &*r_iprime*sin_table(j-jprime)
                end do
            end do
        write(11, '(e20.10)') f_r
        write(12, '(e20.10)') f_theta
        f_r = 0.
        f_theta = 0.
        end do
    end do

! here the calculation has ended
    call CPU_TIME(finish)
    write(*,*) 'Calculation time =', finish-start, 'seconds'

    close(unit=11)
    close(unit=12)

END PROGRAM grav_force

!--------------------------------------------------------------------------------------
!  OPTIMAL TIME STANDS AT:                10.1942997     seconds     -> current standard
! 
!                         without mass    255.950150     seconds     (only slight performance improvements by avoiding too much array calls)
!                            with mass    152.530518     seconds     (half mass calc operations; inner loop)
!     with sine and cosine inside loop    119.338303     seconds     (one cos calculation eliminated; inner loop)
!                     with denominator    57.8149910     seconds     (denominator using sqrt*sqrt*sqrt instead **1.5)
!                  with sin/cos tables    17.7189388     seconds     (calculating cos/sin in external loop)
!                           mass table    15.7216377     seconds     (calculating the mass in external loop)
!                    with r_iprime_cos    14.8603277     seconds     (saves a few calculations in the inner loop since it appears twice there)
! include optimizer compiling flag -O1    10.1942997     seconds     (remembered that gnu compilers have no optimization default)
!    with custom inverse sqrt function    84.1394730     seconds     (sadly in Fortran much slower than intrinsic 1./sqrt, but the error is actually not so bad)
!    with Quakes inverse sqrt function    returns somehow Nan's      (linking the C program in the compilation)
!--------------------------------------------------------------------------------------
! IDEAS:
!       a Fortran version of Quakes inverse square root!
!           -> http://en.wikipedia.org/wiki/Fast_inverse_square_root
!           -> try magic numbers: 0x5fe6ec85e7de30da    or    0x5fe6eb50c7b537a9
!
!--------------------------------------------------------------------------------------

REAL(8) FUNCTION InvSqrt (x)
    IMPLICIT NONE

    TYPE casting
        DOUBLE PRECISION :: x
    END TYPE casting

    REAL(8), INTENT(in) :: x

    ! casting
    TYPE(casting), TARGET :: pointerTo

    ! Encode data as an array of integers
    INTEGER(8), DIMENSION(:), ALLOCATABLE :: enc
    INTEGER(8) :: length
    INTEGER(8) :: magic_number = 6910469410427058089
    REAL(8) :: xhalf

    xhalf = .5*x

    ! transfer to heap
    pointerTo%x = x
    ! encode a memory section from a type to other
    length = size(transfer(pointerTo, enc))
    allocate(enc(length))
    ! encoded to integer
    enc = transfer(pointerTo, enc)  ! evil floating point bit level hacking
    enc(1) = magic_number - rshift(enc(1),1)  ! wtf! (for int64: 0x5fe6eb50c7b537a9 = 6910469410427058089)
    ! decode
    pointerTo = transfer(enc, pointerTo)
    ! dealloc
    deallocate(enc)

    InvSqrt = pointerTo%x*(1.5 - xhalf*pointerTo%x*pointerTo%x)
END FUNCTION InvSqrt
