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
  integer, parameter :: level = 4
  integer, parameter :: level_mult = 2 ** level
  real(8), parameter :: eps = 1e-6
  ! run-time variables
  real(8) :: a,r,r2, theta, r_step, t_step, r_p, theta_p, force_mag_temp, surface_area, d_rd_theta, force_denominator

  ! use zero-based indexing for arrays to make modulo math work better below.
  real(8),dimension(0:(N_r*N_theta)-1) :: force_r, force_theta, force_mag, sigma
  ! Store masses of all levels in the same array with concat. Ex: [level0]+[level1]+[level2]+...
  ! limit of geometric series is the size of array: 4/3 of original size 
  real(8),dimension(0:(N_r*N_theta) * 4 / 3) :: mass
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
  COMPLEX(8) :: force

  ! lookup array
  integer,dimension(0:20) :: level_offset_lookup  ! conseravtively oversized



  write (*,*) 'Program Start ', ' ... '
  write (*,*) 'Data read ...'
  call readFile("r_project.data", r_prime)
  call readFile("theta_project.data", theta_prime)
  call readFile("density_project.data", sigma)
  write (*,*) 'Data read done.'

  ! grid is regular, all step sizes are equal
  r_step = r_prime(1) - r_prime(0)                
  t_step = theta_prime(1) - theta_prime(0)  
  d_rd_theta = r_step * t_step ! precalculate dtheta*dr
   
  ! pre-compute the cache for all theta differences that occur, for N_theta/4 values of N_theta
  do theta_i=0, (N_theta)-1
    cosCache(theta_i) = cos((theta_i+ (0.5 * level_mult))*t_step)
    sinCache(theta_i) = sin((theta_i+ (0.5 * level_mult))*t_step)
  end do

  call computeLevelOffsets(level_offset_lookup)
 
  call CPU_TIME(t_init)

  ! Compute the mass at level_0
  do out_i=0, (N_r*N_theta)-1
    r = r_prime(out_i/N_theta)
    mass(out_i) = sigma(out_i)*r*d_rd_theta
  end do

  call computeMassAtHigherLevels(N_r, N_theta, level)
  call writeMasses(level)

  ! out_i: index to output grid
  do out_i=0, (N_r*N_theta)-1          
                                                                    !avoid singularities
    r=r_prime(out_i/N_theta)-(r_step/2.0)                           !r_prime - half_step
    theta = theta_prime(MODULO(out_i, N_theta))-(t_step/2.0)    !same as above, for theta
    r2=r*r                                                          !putting into memory causes slowdown
    force_r(out_i)=0
    force_theta(out_i)=0
    do in_i=0, inner_loop_size
      force = computeForcesAtLevel(level, level_mult, r, theta, in_i)
      force_r(out_i) = force_r(out_i) + REALPART(force)
      force_theta(out_i) = force_theta(out_i) + IMAGPART(force)
    end do
    force_mag(out_i) = sqrt(force_r(out_i)*force_r(out_i) + force_theta(out_i)*force_theta(out_i)) 
  end do

  call CPU_TIME(t_end)

  call writeFile("force_r.data",force_r)
  call writeFile("force_theta.data",force_theta)
  call writeFile("force_mag.data", force_mag)

  print *, "Time: ", t_end - t_init , "s"


contains 

  subroutine computeLevelOffsets(offsets)
    integer :: offsets(:)
    integer :: i,n,s
    n=size(offsets)
    s=N_r*N_theta
    offsets(1) = 0
    do i=2, n
      offsets(i) = offsets(i-1) + s
      s = s/4
    end do
  end subroutine computeLevelOffsets


  ! using the COMPLEX type as a container for the f_r and f_t return value
  recursive function computeForcesAtLevel(level, level_mult, r, theta, in_i) result(force)
    implicit none
    ! input variables
    integer, intent(in) :: level, level_mult, in_i
    real(8), intent(in) :: r, theta
    complex(8) :: force
    real(8) :: level_cutoff
    integer :: r_lookup, in_r, in_t
    real(8) :: f_mag, f_r, f_t
    !
    level_cutoff = level_mult * 2.0
    r_lookup = in_i / (N_theta/level_mult)
    r_p = (r_prime(r_lookup * level_mult) + &
             r_prime((r_lookup+1) * level_mult - 1)) * 0.5
    theta_p = (theta_prime(MODULO(in_i*level_mult, N_theta)) + &
               theta_prime(MODULO((in_i+1)*level_mult - 1, N_theta))) * 0.5
    if (level > 1 .AND. (abs(r-r_p) > level_cutoff * r_step .OR. abs(theta-theta_p) > level_cutoff * t_step)) THEN
      in_r = in_i / (N_theta / level_mult)
      in_t = MODULO(in_i, N_theta / level_mult)
      force =         computeForcesAtLevel(level-1,level_mult/2, r, theta, (in_r  ) * 4 + in_t * 2)
      force = force + computeForcesAtLevel(level-1,level_mult/2, r, theta, (in_r  ) * 4 + in_t * 2 + 1)
      force = force + computeForcesAtLevel(level-1,level_mult/2, r, theta, (in_r+1) * 4 + in_t * 2)
      force = force + computeForcesAtLevel(level-1,level_mult/2, r, theta, (in_r+1) * 4 + in_t * 2 + 1)
    ELSE
      !cosvar = cosCache(MODULO(in_i*level_mult-out_i, N_theta))  
      !sinvar = sinCache(MODULO(in_i*level_mult-out_i, N_theta))
      cosvar = cos(theta_p - theta)
      sinvar = sin(theta_p - theta)

      force_denominator = (r*r+r_p*r_p-(2*r*r_p*cosvar))
      force_denominator = sqrt(force_denominator)*force_denominator + eps

      ! TODO: fix mass lookup
      f_mag = mass(level_offset_lookup(level) + in_i)/force_denominator
      f_r = (r - r_p*cosvar) * f_mag
      f_t = (r_p*sinvar) * f_mag
      force = COMPLEX(f_r, f_t)
    END IF
  end function computeForcesAtLevel

  subroutine computeMassAtHigherLevels(N_r, N_theta, level)
    !real(8) :: mass(:)                          ! In place replacement of array. Replace Level L with L+level
    integer :: N_r, N_theta, level              ! N_r and N_theta are sizes of input data
    integer :: in_N_r, in_N_t, out_N_r, out_N_t ! size of the array at old and new level
    integer :: i, r, t                          ! iterator variables
    integer :: offset_in, offset_out            ! offsets
    in_N_r = N_r
    in_N_t = N_theta

    offset_out = N_r * N_theta
    offset_in = 0

    do i = 1, level
      out_N_r = in_N_r/2
      out_N_t = in_N_t/2
      write(*,*) i, offset_out, offset_in, out_N_r, out_N_t
        
      do r = 0, out_N_r-1
        do t = 0, out_N_t-1
          !write(*,*) r, t, 2*t+out_N_t*r*2, 4*t +     (in_N_t*2*r)
          mass(t+out_N_t*r+offset_out) = mass(offset_in+2*t +     (in_N_t*2*r)) + &
                                         mass(offset_in+2*t + 1 + (in_N_t*2*r)) + &
                                         mass(offset_in+2*t +     (in_N_t*(2*r + 1))) + &
                                         mass(offset_in+2*t + 1 + (in_N_t*(2*r + 1)))
        end do
      end do
      offset_in = offset_out
      offset_out = offset_out + out_N_r * out_N_t
        
      in_N_r = out_N_r
      in_N_t = out_N_t
    end do
  end subroutine computeMassAtHigherLevels

  subroutine writeMasses(level)
    integer :: level
    integer :: i
    integer :: begin, end, cur_size
    character(len=128) :: filename
    begin = 1
    cur_size = N_r * N_theta
    end = cur_size
    do i=0, level
      write(filename, "(A5,I0.2,A4)") "mass_", i, ".txt"
      print *, filename
      call writeFilePartly(filename, mass, begin, end)
      begin = end + 1
      cur_size = cur_size / 4
      end = end + cur_size
    end do
  end subroutine writeMasses

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

  ! write subset of the array
  subroutine writeFilePartly(filename, x, begin, end)
    character(*) :: filename
    real(8) :: x(:)
    integer :: begin, end
    logical :: exists

    integer :: i,n
    real(8) :: temp
    inquire(file=filename, exist=exists)
    ! merge is the same as the terneray operator from C: (exists ? 'new' : 'old')
    open(unit=8, file=filename, status=merge('old', 'new', exists), action='write')
    do i=begin, end
        write(8, '(e20.10)') x(i)
    end do
    close(8)
  end subroutine writeFilePartly

  ! write the full array.
  subroutine writeFile(filename, x)
    character(*) :: filename
    real(8) :: x(:)
    call writeFilePartly(filename, x, 1, size(x))
  end subroutine writeFile

end program gravity
