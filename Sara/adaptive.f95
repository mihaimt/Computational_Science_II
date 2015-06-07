program gravity

  !  optimised with -O2
  !  2 trials
          
  implicit none             ! all variables must be defined
  real(8) t_init, t_end     ! global timing
  real, parameter :: G = 1  ! normalised to 1; 6.67384e-11 !(m*m*m)/(kg*s*s)
  real, parameter :: pi = 3.1415926538
  integer, parameter :: N_r = 128
  integer, parameter :: N_theta = 256
  integer, parameter :: min_level = 0
  integer, parameter :: max_level = 6
  integer, parameter :: level_mult = 2 ** max_level
  real(8), parameter :: eps = 0
  ! run-time variables
  real(8) :: a,r,r2, theta, r_step, t_step, r_p, theta_p, force_mag_temp, surface_area, d_rd_theta, force_denominator

  ! use zero-based indexing for arrays to make modulo math work better below.
  real(8),dimension(0:(N_r*N_theta)-1) :: force_r, force_theta, force_mag, sigma, debug_value
  ! Store masses of all levels in the same array with concat. Ex: [level0]+[level1]+[level2]+...
  ! limit of geometric series is the size of array: 4/3 of original size 
  real(8),dimension(0:(N_r*N_theta) * 4 / 3) :: mass
  real(8) :: r_prime(0:N_r-1), theta_prime(0:N_theta-1)
  real(8) :: u_r, u_theta ! the unit vectors

  ! iterator variables
  integer :: out_i, in_i, theta_i

 
  !size of the inner loop, assume start at 0
  integer, parameter :: inner_loop_size = ((N_r / level_mult) * (N_theta / level_mult)) - 1

  ! temp variables
  integer :: r_lookup
  COMPLEX(8) :: force

  ! lookup array
  integer,dimension(0:20) :: level_offset_lookup  ! conservatively oversized

  call readFile("r_project.data", r_prime)
  call readFile("theta_project.data", theta_prime)
  call readFile("density_project.data", sigma)
  print *, "Simulation start..."

  ! grid is regular, all step sizes are equal
  r_step = r_prime(1) - r_prime(0)
  t_step = theta_prime(1) - theta_prime(0)  
  d_rd_theta = r_step * t_step ! precalculate dtheta*dr
   
  call computeLevelOffsets(level_offset_lookup)
 
  call CPU_TIME(t_init)

  ! Compute the mass at level_0
  do out_i=0, (N_r*N_theta)-1
    r = r_prime(out_i/N_theta)
    mass(out_i) = -sigma(out_i)*r*d_rd_theta
  end do

  call computeMassAtHigherLevels(N_r, N_theta, max_level)
  !call writeMasses(max_level)

  ! out_i: index to output grid
  do out_i=0, (N_r*N_theta)-1
    ! The output grid is half a theta and r step shifted.
    r=r_prime(out_i/N_theta)-(r_step/2.0)
    theta = theta_prime(MODULO(out_i, N_theta))-(t_step/2.0)
    force_r(out_i)=0
    force_theta(out_i)=0
    do in_i=0, inner_loop_size
      force = computeForcesAtLevel(max_level, level_mult, r, theta, in_i, out_i == 12401)
      force_r(out_i) = force_r(out_i) + REALPART(force)
      force_theta(out_i) = force_theta(out_i) + IMAGPART(force)
    end do
    force_mag(out_i) = sqrt(force_r(out_i)*force_r(out_i) + &
                            force_theta(out_i)*force_theta(out_i)) 
  end do

  call CPU_TIME(t_end)

  call writeFile("force_r.data",force_r)
  call writeFile("force_theta.data",force_theta)
  call writeFile("force_mag.data", force_mag)
  !call writeFile("debug_level.data", debug_value)

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

  ! Used to see which level contributes to a specific pixel.
  recursive subroutine showDebug(in_i, in_level, read_level)
    integer, intent(in) :: in_i, in_level, read_level
    integer :: level_mult, in_r, in_t, in_j
    level_mult = 2 ** in_level
    if (in_level > 0) then
      ! This uses the same recursion logic as computeForcesAtLevel to propagate the read level
      ! down to every individual grid cell (=pixel)
      in_r = in_i / (N_theta / level_mult)
      in_t = MODULO(in_i, N_theta / level_mult)
      in_j = (in_r * 2) * (2*N_theta / level_mult) + in_t * 2
      call showDebug(in_j, in_level-1, read_level)
      call showDebug(in_j + 1, in_level-1, read_level)
      in_j = (in_r * 2 + 1) * (2*N_theta / level_mult) + in_t * 2
      call showDebug(in_j, in_level-1, read_level)
      call showDebug(in_j + 1, in_level-1, read_level)
    else
      debug_value(in_i) = dble(read_level)
    end if
  end subroutine showDebug

  ! This method computes the force between a given (r,theta) and the mass at (in_i,
  ! which defines r' and theta') at a given level. If this level is too big/close to
  ! the current (r,theta) then the function recurses one level lower until the
  ! the distance constraint is no longer violated or the lowest level has been
  ! reached.
  ! HACK: Using the COMPLEX type as a container for the f_r and f_t return value
  recursive function computeForcesAtLevel(level, level_mult, r, theta, in_i, debug_this) result(force)
    implicit none
    ! input variables
    integer, intent(in) :: level, level_mult, in_i
    real(8), intent(in) :: r, theta
    complex(8) :: force
    real(8) :: level_cutoff, theta_diff
    integer :: r_lookup, in_r, in_t, in_j
    real(8) :: f_mag, f_r, f_t
    real(8) :: cosvar, sinvar
    logical :: debug_this ! Used to display what levels are used for this specific point.
    ! Grid cell distance that allows this level to be used
    !  Arbitrarily chosen value, try trial and error to find best rate of trade-off:
    !  low error, good speed. try: [.5, 3] ?
    level_cutoff = level_mult * 3.0
    r_lookup = in_i / (N_theta/level_mult)
    r_p = (r_prime(r_lookup * level_mult) + &
           r_prime((r_lookup+1) * level_mult - 1)) * 0.5
    theta_p = (theta_prime(MODULO(in_i*level_mult, N_theta)) + &
               theta_prime(MODULO((in_i+1)*level_mult - 1, N_theta))) * 0.5

    ! properly handle the cyclical nature of the theta dimension.
    theta_diff = abs(theta - theta_p)
    theta_diff = min(theta_diff, abs(theta_diff - (2 * pi)))
    IF (level > min_level .AND. &
        (abs(r-r_p) < level_cutoff * r_step .AND. &
         theta_diff < level_cutoff * t_step)) THEN
      in_r = in_i / (N_theta / level_mult)
      in_t = MODULO(in_i, N_theta / level_mult)
      in_j = (in_r * 2) * (2*N_theta / level_mult) + in_t * 2
      force =         computeForcesAtLevel(level-1,level_mult/2, r, theta, in_j, debug_this)
      force = force + computeForcesAtLevel(level-1,level_mult/2, r, theta, in_j + 1, debug_this)
      in_j = (in_r * 2 + 1) * (2*N_theta / level_mult) + in_t * 2
      force = force + computeForcesAtLevel(level-1,level_mult/2, r, theta, in_j, debug_this)
      force = force + computeForcesAtLevel(level-1,level_mult/2, r, theta, in_j + 1, debug_this)
    ELSE
      ! TODO: Implement the cos cache again, using cos/sin is super slow.
      cosvar = cos(theta - theta_p)
      sinvar = sin(theta - theta_p)

      force_denominator = (r*r+r_p*r_p-(2*r*r_p*cosvar))
      force_denominator = sqrt(force_denominator)*force_denominator + eps

      f_mag = mass(level_offset_lookup(level) + in_i)/force_denominator
      f_r = (r - r_p*cosvar) * f_mag
      f_t = (r_p*sinvar) * f_mag
      force = COMPLEX(f_r, f_t)

      ! Only turn on for debugging
      !IF (debug_this) THEN
      !  call showDebug(in_i, level, level)
      !END IF
    END IF
  end function computeForcesAtLevel

  ! All masses are in the same array, the first level is at index [0, N_r * N_theta),
  ! the second level is at [N_r * N_theta, (N_r * N_theta) + (N_r/2) * (N_theta/2) ),
  ! and so forth. Every higher level is only 1/4th of the size of the previous one.
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
