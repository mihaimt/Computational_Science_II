program gravity
  !  optimised with -O2
  !  2 trials
!      Time:      3.6949229999999997      s
!      Time:      3.6917869999999997      s
          
  implicit none             ! all variables must be defined
  real(8) t_init, t_end     ! global timing
  real, parameter :: G = 1  ! normalised to 1; 6.67384e-11 !(m*m*m)/(kg*s*s)
  real, parameter :: pi = 3.1415926538
  integer, parameter :: N_r = 128
  integer, parameter :: N_theta = 256
  real(8), parameter :: epsilonCos = 1d-8
  integer, parameter :: level = 1
  integer, parameter :: level_mult = 2 ** level
  real(8), parameter :: eps = 1e-6
  ! run-time variables
  real(8) :: a,r,r2, theta, r_step, theta_step, r_p, theta_p, force_mag_temp, surface_area, d_rd_theta, force_denominator

  ! use zero-based indexing for arrays to make modulo math work better below.
  real(8),dimension(0:(N_r*N_theta)-1) :: force_r, force_theta, force_mag, sigma, mass
  real(8) :: r_prime(0:N_r-1), theta_prime(0:N_theta-1)
  real(8) :: u_r, u_theta ! the unit vectors

  ! iterator variables
  integer :: out_i, in_i, theta_i

  !cache and lookup variables
  real(8) :: cosvar, sinvar, cosCache(0:N_theta-1), sinCache(0:N_theta-1)
 
  !size of the inner loop, assume start at 0
  integer, parameter :: inner_loop_size = ((N_r / level_mult) * (N_theta / level_mult)) - 1

  ! temp variables
  integer :: r_lookup

  write (*,*) 'Program Start ', ' ... '
  write (*,*) 'Data read ...'
  call readFile("r_project.data", r_prime)
  call readFile("theta_project.data", theta_prime)
  call readFile("density_project.data", sigma)
  write (*,*) 'Data read done.'

  ! grid is regular, all step sizes are equal
  r_step = r_prime(1) - r_prime(0)                
  theta_step = theta_prime(1) - theta_prime(0)  
  d_rd_theta = r_step * theta_step ! precalculate dtheta*dr
   
  ! pre-compute the cache for all theta differences that occur, for N_theta/4 values of N_theta
  do theta_i=0, (N_theta)-1
    cosCache(theta_i) = cos((theta_i+ (0.5 * level_mult))*theta_step)
    sinCache(theta_i) = sin((theta_i+ (0.5 * level_mult))*theta_step)
  end do

  call CPU_TIME(t_init)

  ! Compute the mass at level_0
  do out_i=0, (N_r*N_theta)-1
    r = r_prime(out_i/N_theta)
    mass(out_i) = -sigma(out_i)*r*d_rd_theta
  end do

  call increaseLevel(mass, N_r, N_theta, level)

  ! out_i: index to output grid
  do out_i=0, (N_r*N_theta)-1          
                                                                    !avoid singularities
    r=r_prime(out_i/N_theta)-(r_step/2.0)                           !r_prime - half_step
    theta = theta_prime(MODULO(out_i, N_theta))-(theta_step/2.0)    !same as above, for theta
    r2=r*r                                                          !putting into memory causes slowdown
    force_r(out_i)=0
    force_theta(out_i)=0
    do in_i=0, inner_loop_size
      ! current r-prime is the average of the min and max r-prime of the level
      r_lookup = in_i / (N_theta/level_mult)
      r_p = (r_prime(r_lookup * level_mult) + &
             r_prime((r_lookup+1) * level_mult - 1)) * 0.5
      ! current theta-prime
      theta_p = (theta_prime(MODULO(in_i*level_mult, N_theta)) + &
                 theta_prime(MODULO((in_i+1)*level_mult - 1, N_theta))) * 0.5
      
      cosvar = cosCache(MODULO(out_i-in_i*level_mult, N_theta))  
      sinvar = sinCache(MODULO(out_i-in_i*level_mult, N_theta))

      force_denominator = (r2+r_p*r_p-(2*r*r_p*cosvar))
      force_denominator = sqrt(force_denominator)*force_denominator + eps

      force_mag_temp = mass(in_i)/force_denominator
      force_r(out_i) = (r - r_p*cosvar) * force_mag_temp + force_r(out_i)
      force_theta(out_i) = (r_p*sinvar) * force_mag_temp + force_theta(out_i)
    end do
    force_mag(out_i) = sqrt(force_r(out_i)*force_r(out_i) + force_theta(out_i)*force_theta(out_i)) 
  end do

  call CPU_TIME(t_end)

  call writeFile("force_r.data",force_r)
  call writeFile("force_theta.data",force_theta)
  call writeFile("force_mag.data", force_mag)

  print *, "Time: ", t_end - t_init , "s"


contains 

  subroutine increaseLevel(mass, N_r, N_theta, level)
    real(8) :: mass(:)          ! In place replacement of array. Replace Level L with L+level
    integer :: N_r, N_theta, level  ! N_r and N_theta are sizes of input data
    integer :: in_N_r, in_N_t, out_N_r, out_N_t ! size of the array at old and new level
    integer :: i, r, t          ! iterator variables
    in_N_r = N_r
    in_N_t = N_theta

    do i = 1, level
      out_N_r = in_N_r/2
      out_N_t = in_N_t/2    
        
      do r = 0, out_N_r-1
        do t = 0, out_N_t-1
          mass(t+out_N_t*r) = mass(2*t +     (in_N_t*2*r)) + &
                              mass(2*t + 1 + (in_N_t*2*r)) + &
                              mass(2*t +     (in_N_t*(2*r + 1))) + &
                              mass(2*t + 1 + (in_N_t*(2*r + 1)))
        end do
      end do
        
      in_N_r = out_N_r
      in_N_t = out_N_t
    end do
  end subroutine increaseLevel

  subroutine readFile(filename, x)
    character(*) :: filename
    real(8) :: x(:)

    integer :: i,n
    real(8) :: temp
    n=size(x)
    open(unit=8, file=filename, status='old', action='read')
    do i=1, n
        read(8, '(e20.10)') temp
        x(i) = temp
    end do
    close(8)
  end subroutine readFile

  subroutine writeFile(filename, x)
    character(*) :: filename
    real(8) :: x(:)

    integer :: i,n
    real(8) :: temp
    n=size(x)
    open(unit=8, file=filename, status='old', action='write')
    do i=1, n
        write(8, '(e20.10)') x(i)
    end do
    close(8)
  end subroutine writeFile

end program gravity
