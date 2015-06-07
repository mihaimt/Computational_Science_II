PROGRAM gravity
  IMPLICIT NONE             ! all variables must be defined
  REAL(8) t_init, t_end     ! global timing
  REAL, PARAMETER :: G = 1  ! normalised to 1; 6.67384e-11 !(m*m*m)/(kg*s*s)
  REAL, PARAMETER :: pi = 3.1415926538
  INTEGER, PARAMETER :: N_r = 128
  INTEGER, PARAMETER :: N_theta = 256
  INTEGER :: min_level = 0
  INTEGER :: max_level = 6
  INTEGER :: level_mult 
  REAL(8), PARAMETER :: eps = 1e-6
  ! run-time variables
  REAL(8) :: a,r,r2, theta, r_step, t_step, force_mag_temp, surface_area, d_rd_theta, force_denominator
  REAL(8) :: mass_diff_threshold, mass_diff_percentage, mass_diff_total
  LOGICAL :: com_correction

  ! use zero-based indexing for arrays to make modulo math work better below.
  REAL(8),DIMENSION(0:(N_r*N_theta)-1) :: force_r, force_theta, force_mag, sigma, debug_value
  ! Store masses of all levels in the same array with concat. Ex: [level0]+[level1]+[level2]+...
  ! limit of geometric series is the size of array: 4/3 of original size 
  REAL(8),DIMENSION(0:(N_r*N_theta) * 4 / 3) :: mass, mass_min, mass_max, mass_offset_r, mass_offset_t
  REAL(8) :: r_prime(0:N_r-1), theta_prime(0:N_theta-1)

  ! iterator variables
  INTEGER :: out_i, in_i, theta_i


  !size of the inner loop, assume start at 0
  INTEGER :: inner_loop_size

  ! temp variables
  INTEGER :: r_lookup
  COMPLEX(8) :: force

  CHARACTER(len=64) :: arg, output_prefix

  ! lookup array
  INTEGER,DIMENSION(0:20) :: level_offset_lookup  ! conservatively oversized

  debug_value(:) = -1.0

  ! The first command line parameter will be the max level
  IF (iargc() > 0) THEN
    CALL getarg(1, arg)
    READ(arg, '(i10)') max_level
  END IF
  ! The second command line parameter will be the min level
  IF (iargc() > 1) THEN
    CALL getarg(2, arg)
    READ(arg, '(i10)') min_level
  END IF
  PRINT *, "MaxLevel:", max_level
  PRINT *, "MinLevel:", min_level
  IF (min_level > max_level) THEN
    PRINT *, "MaxLevel must be larger or equal to MinLevel"
    CALL EXIT(1)
  END IF

  ! The third command line parameter will be the mass difference threshold
  IF (iargc() > 2) THEN
    CALL getarg(3, arg)
    READ(arg, *) mass_diff_percentage
    IF (mass_diff_percentage > 0.0) THEN
      PRINT *, "Adaptive mass subsampling turned on", mass_diff_percentage * 100, "%"
    END IF
  ELSE
    mass_diff_percentage = 0.0
  END IF

  com_correction = .FALSE.
  IF (iargc() > 3) THEN
    CALL getarg(4, arg)
    IF (arg == "com") THEN
      com_correction = .TRUE.
      print *, "Center of Mass correction turned on"
    END IF
  END IF

  CALL readFile("r_project.data", r_prime)
  CALL readFile("theta_project.data", theta_prime)
  ! The fifth command line parameter will be the input file
  IF (iargc() > 4) THEN
    CALL getarg(5, arg)
    CALL readFile(arg, sigma)
  ELSE
    CALL readFile("density_project.data", sigma)
  END IF

  ! The sixth command line parameter will be the output prefix
  IF (iargc() > 5) THEN
    CALL getarg(6, output_prefix)
  END IF

  PRINT *, "Simulation start..."

  ! grid is regular, all step sizes are equal
  r_step = r_prime(1) - r_prime(0)
  t_step = theta_prime(1) - theta_prime(0)  
  d_rd_theta = r_step * t_step ! precalculate dtheta*dr
   
  CALL computeLevelOffsets(level_offset_lookup)
 
  CALL CPU_TIME(t_init)


  CALL computeMassAtLevel0()
  CALL computeMassAtHigherLevels(N_r, N_theta, max_level)
  CALL writeMasses(max_level)

  level_mult = 2 ** max_level
  inner_loop_size = ((N_r / level_mult) * (N_theta / level_mult)) - 1

  ! out_i: index to output grid
  DO out_i=0, (N_r*N_theta)-1
    IF (MODULO(out_i, N_r) == N_r - 1) THEN
      ! This costs about 0.03 seconds, totally worth it.
      call progress(100 * out_i / (N_r * N_theta - 1))
    END IF
    ! The output grid is half a theta and r step shifted.
    r=r_prime(out_i/N_theta)-(r_step/2.0)
    theta = theta_prime(MODULO(out_i, N_theta))-(t_step/2.0)
    force_r(out_i)=0
    force_theta(out_i)=0
    DO in_i=0, inner_loop_size
      force = computeForcesAtLevel(max_level, level_mult, r, theta, in_i, out_i == 21098)
      force_r(out_i) = force_r(out_i) + REALPART(force)
      force_theta(out_i) = force_theta(out_i) + IMAGPART(force)
    END DO
    force_mag(out_i) = sqrt(force_r(out_i)*force_r(out_i) + &
                            force_theta(out_i)*force_theta(out_i)) 
  END DO

  CALL CPU_TIME(t_end)

  CALL writeFile("force_r.data",force_r)
  CALL writeFile("force_theta.data",force_theta)
  CALL writeFile("force_mag.data", force_mag)
  CALL writeFile("level.data", debug_value)

  PRINT *, "Time: ", t_end - t_init , "s"

contains 

  SUBROUTINE computeLevelOffsets(offsets)
    INTEGER :: offsets(:)
    INTEGER :: i,n,s
    n=size(offsets)
    s=N_r*N_theta
    offsets(1) = 0
    DO i=2, n
      offsets(i) = offsets(i-1) + s
      s = s/4
    END DO
  END SUBROUTINE computeLevelOffsets

  ! Used to see which level contributes to a specific pixel.
  RECURSIVE SUBROUTINE showDebug(in_i, in_level, read_level)
    INTEGER, INTENT(IN) :: in_i, in_level, read_level
    INTEGER :: level_mult, in_r, in_t, in_j
    level_mult = 2 ** in_level
    IF (in_level > 0) THEN
      ! This uses the same recursion logic as computeForcesAtLevel to propagate the read level
      ! down to every individual grid cell (=pixel)
      in_r = in_i / (N_theta / level_mult)
      in_t = MODULO(in_i, N_theta / level_mult)
      in_j = (in_r * 2) * (2*N_theta / level_mult) + in_t * 2
      CALL showDebug(in_j, in_level-1, read_level)
      CALL showDebug(in_j + 1, in_level-1, read_level)
      in_j = (in_r * 2 + 1) * (2*N_theta / level_mult) + in_t * 2
      CALL showDebug(in_j, in_level-1, read_level)
      CALL showDebug(in_j + 1, in_level-1, read_level)
    ELSE
      debug_value(in_i) = dble(read_level)
    END IF
  END SUBROUTINE showDebug

  ! This method computes the force between a given (r,theta) and the mass at (in_i,
  ! which defines r' and theta') at a given level. If this level is too big/close to
  ! the current (r,theta) then the function recurses one level lower until the
  ! the distance constraint is no longer violated or the lowest level has been
  ! reached.
  ! HACK: Using the COMPLEX type as a container for the f_r and f_t return value
  RECURSIVE FUNCTION computeForcesAtLevel(level, level_mult, r, theta, in_i, debug_this) result(force)
    IMPLICIT NONE
    ! input variables
    INTEGER, INTENT(IN) :: level, level_mult, in_i
    REAL(8), INTENT(IN) :: r, theta
    complex(8) :: force
    REAL(8) :: level_cutoff, theta_diff
    INTEGER :: r_lookup, in_r, in_t, in_j, mass_lookup
    REAL(8) :: r_p, t_p, f_mag, f_r, f_t, mass_diff
    REAL(8) :: cosvar, sinvar
    LOGICAL :: recurseFurther
    LOGICAL :: debug_this ! Used to display what levels are used for this specific point.
    ! Grid cell distance that allows this level to be used
    !  Arbitrarily chosen value, try trial and error to find best rate of trade-off:
    !  low error, good speed. try: [.5, 3] ?
    level_cutoff = level_mult * 3.0
    r_lookup = in_i / (N_theta/level_mult)
    r_p = (r_prime(r_lookup * level_mult) + &
           r_prime((r_lookup+1) * level_mult - 1)) * 0.5
    t_p = (theta_prime(MODULO(in_i*level_mult, N_theta)) + &
           theta_prime(MODULO((in_i+1)*level_mult - 1, N_theta))) * 0.5

    mass_lookup = level_offset_lookup(level) + in_i
    recurseFurther = .FALSE.
    IF (level > min_level) THEN
      ! properly handle the cyclical nature of the theta dimension.
      theta_diff = abs(theta - t_p)
      theta_diff = min(theta_diff, abs(theta_diff - (2 * pi)))
      recurseFurther = abs(r-r_p) < level_cutoff * r_step .AND. theta_diff < level_cutoff * t_step

      IF (mass_diff_percentage > 0) THEN
        ! the index into the mass arrays.
        mass_diff = mass_max(mass_lookup) - mass_min(mass_lookup)
        ! masses are negative (-sigma) so min is the important one
        ! this value must be tuned (or precomputed?)
        recurseFurther = recurseFurther .OR. mass_diff > mass_diff_threshold
      END IF
    END IF
    IF (recurseFurther) THEN
      in_r = in_i / (N_theta / level_mult)
      in_t = MODULO(in_i, N_theta / level_mult)
      in_j = (in_r * 2) * (2*N_theta / level_mult) + in_t * 2
      force =         computeForcesAtLevel(level-1,level_mult/2, r, theta, in_j, debug_this)
      force = force + computeForcesAtLevel(level-1,level_mult/2, r, theta, in_j + 1, debug_this)
      in_j = (in_r * 2 + 1) * (2*N_theta / level_mult) + in_t * 2
      force = force + computeForcesAtLevel(level-1,level_mult/2, r, theta, in_j, debug_this)
      force = force + computeForcesAtLevel(level-1,level_mult/2, r, theta, in_j + 1, debug_this)
    ELSE
      IF (com_correction) THEN
        r_p = r_p + mass_offset_r(mass_lookup)
        t_p = t_p + mass_offset_t(mass_lookup)
      END IF
      ! The cosine cache might still be a very nice optimization, unfortunately with the center of mass
      ! correction all theta differences might be unique, at which point a cache is completely useless.
      cosvar = cos(theta - t_p)
      sinvar = sin(theta - t_p)

      force_denominator = (r*r+r_p*r_p-(2*r*r_p*cosvar))
      force_denominator = sqrt(force_denominator)*force_denominator + eps

      f_mag = mass(mass_lookup) / force_denominator
      f_r = (r - r_p*cosvar) * f_mag
      f_t = (r_p*sinvar) * f_mag
      force = COMPLEX(f_r, f_t)

      ! Only turn on for debugging
      IF (debug_this) THEN
        CALL showDebug(in_i, level, level)
      END IF
    END IF
  END FUNCTION computeForcesAtLevel

  SUBROUTINE computeMassAtLevel0()
    ! Compute the mass at level_0
    INTEGER :: out_i, in_i, theta_i
    REAL(8) :: smallest_mass, largest_mass, cur_mass
    smallest_mass = 1e20
    largest_mass = -1e20
    DO out_i=0, (N_r*N_theta)-1
      r = r_prime(out_i/N_theta)
      cur_mass = -sigma(out_i) * r * d_rd_theta
      smallest_mass = min(smallest_mass, cur_mass)
      largest_mass = max(largest_mass, cur_mass)
      mass(out_i) = cur_mass
      mass_min(out_i) = cur_mass
      mass_max(out_i) = cur_mass
      mass_offset_r(out_i) = 0.0
      mass_offset_t(out_i) = 0.0
    END DO

    mass_diff_total = largest_mass - smallest_mass
    IF (mass_diff_percentage > 0) THEN
      mass_diff_threshold = mass_diff_percentage * mass_diff_total
    ELSE
      mass_diff_threshold = 0
    END IF
  END SUBROUTINE computeMassAtLevel0

  ! All masses are in the same array, the first level is at index [0, N_r * N_theta),
  ! the second level is at [N_r * N_theta, (N_r * N_theta) + (N_r/2) * (N_theta/2) ),
  ! and so forth. Every higher level is only 1/4th of the size of the previous one.
  SUBROUTINE computeMassAtHigherLevels(N_r, N_theta, level)
    INTEGER :: N_r, N_theta, level              ! N_r and N_theta are sizes of input data
    INTEGER :: in_N_r, in_N_t, out_N_r, out_N_t ! size of the array at old and new level
    INTEGER :: i, r, t                          ! iterator variables
    INTEGER :: i1, i2, i3, i4, o
    INTEGER :: offset_in, offset_out            ! offsets
    REAL(8) :: min_1, min_2, max_1, max_2
    REAL(8) :: half_step_r, half_step_t
    REAL(8) :: center_r, center_t
    in_N_r = N_r
    in_N_t = N_theta

    offset_out = N_r * N_theta
    offset_in = 0

    half_step_r = r_step / 2.0
    half_step_t = t_step / 2.0

    DO i = 1, level
      out_N_r = in_N_r/2
      out_N_t = in_N_t/2
        
      DO r = 0, out_N_r-1
        DO t = 0, out_N_t-1
          i1 = (offset_in + 2 * t    ) + (in_N_t * (2 * r    ))
          i2 = (offset_in + 2 * t + 1) + (in_N_t * (2 * r    ))
          i3 = (offset_in + 2 * t    ) + (in_N_t * (2 * r + 1))
          i4 = (offset_in + 2 * t + 1) + (in_N_t * (2 * r + 1))
          o = offset_out + t + out_N_t * r
          mass(o) = mass(i1) + mass(i2) + mass(i3) + mass(i4)
          
          ! Compute the spread between the lowest and highest mass of all the level 0 cells that are 
          ! contained by this one.
          min_1 = min(mass_min(i1), mass_min(i2))
          min_2 = min(mass_min(i3), mass_min(i4))
          mass_min(o) = min(min_1, min_2)
          max_1 = max(mass_max(i1), mass_max(i2))
          max_2 = max(mass_max(i3), mass_max(i4))
          mass_max(o) = max(max_1, max_2)

          ! Compute the location of the center of mass of the aggregated cell.

          ! Careful here, if the mass(o) is too small, dividing by it can yield very large results that will
          ! over-correct and lead to numerical instability.
          IF (abs(mass(o)) > 0.01 * mass_diff_total) THEN
            ! This summation doesn't account for polar coordinates and is therefore at best an
            ! approximation. The error will increase in higher levels as each cell covers a larger arc.
            center_r =            (mass_offset_r(i1) - half_step_r) * mass(i1)
            center_r = center_r + (mass_offset_r(i2) - half_step_r) * mass(i2)
            center_r = center_r + (mass_offset_r(i3) + half_step_r) * mass(i3)
            center_r = center_r + (mass_offset_r(i4) + half_step_r) * mass(i4)
            
            ! For angles it's supposed to work better, but is probably succeptible to the same mistake.
            center_t =            (mass_offset_t(i1) - half_step_t) * mass(i1)
            center_t = center_t + (mass_offset_t(i2) + half_step_t) * mass(i2)
            center_t = center_t + (mass_offset_t(i3) - half_step_t) * mass(i3)
            center_t = center_t + (mass_offset_t(i4) + half_step_t) * mass(i4)
            mass_offset_r(o) = center_r / mass(o)
            mass_offset_t(o) = center_t / mass(o)
          ELSE
            mass_offset_r(o) = 0
            mass_offset_t(o) = 0
          END IF
        END DO
      END DO
      offset_in = offset_out
      offset_out = offset_out + out_N_r * out_N_t
        
      in_N_r = out_N_r
      in_N_t = out_N_t
      half_step_r = half_step_r * 2.0
      half_step_t = half_step_t * 2.0
    END DO
  END SUBROUTINE computeMassAtHigherLevels

  ! Write the masses of all levels computed to separate files.
  SUBROUTINE writeMasses(level)
    INTEGER :: level
    INTEGER :: i
    INTEGER :: begin, end, cur_size
    CHARACTER(len=128) :: filename
    begin = 1
    cur_size = N_r * N_theta
    end = cur_size
    DO i=0, level
      write(filename, "(A11,I0.2,A4)") "mass_", i, ".txt"
      !PRINT *, filename
      CALL writeFilePartly(filename, mass, begin, end)
      begin = end + 1
      cur_size = cur_size / 4
      end = end + cur_size
    END DO
  END SUBROUTINE writeMasses

  SUBROUTINE readFile(filename, x)
    CHARACTER(*) :: filename
    REAL(8) :: x(:)

    INTEGER :: i,n
    REAL(8) :: temp
    n=size(x)
    OPEN(unit=8, file=filename, status='old', action='read')
    DO i=1, n
        READ(8, '(e20.10)') temp
        x(i) = temp
    END DO
    CLOSE(8)
  END SUBROUTINE readFile

  ! write subset of the array
  SUBROUTINE writeFilePartly(filename, x, begin, end)
    CHARACTER(*) :: filename
    REAL(8) :: x(:)
    INTEGER :: begin, end
    logical :: exists

    CHARACTER(128) :: full_filename

    INTEGER :: i,n
    REAL(8) :: temp

    full_filename = trim(output_prefix)//trim(adjustl(filename))
    inquire(file=full_filename, exist=exists)
    ! merge is the same as the terneray operator from C: (exists ? 'new' : 'old')
    OPEN(unit=8, file=full_filename, status=merge('old', 'new', exists), action='write')
    DO i=begin, end
        write(8, '(e20.10)') x(i)
    END DO
    CLOSE(8)
  END SUBROUTINE writeFilePartly

  ! write the full array.
  SUBROUTINE writeFile(filename, x)
    CHARACTER(*) :: filename
    REAL(8) :: x(:)
    CALL writeFilePartly(filename, x, 1, size(x))
  END SUBROUTINE writeFile

  
  ! a simple progress bar.
  SUBROUTINE progress(j)
    IMPLICIT NONE
    INTEGER(KIND=4)::j,k
    CHARACTER(LEN=57)::bar="???% |                                                  |"
    WRITE(unit=bar(1:3),fmt="(i3)") j
    IF (j > 2) THEN
      bar(7 : 6 + (j/2)) = "**************************************************"
    END IF
    ! print the progress bar.
    WRITE(unit=6,fmt="(a1,a57)",advance="no") char(13), bar
    IF (j /= 100) THEN
      FLUSH(unit=6)
    ELSE
      WRITE(unit=6,fmt=*)
    ENDIF
    RETURN
  END SUBROUTINE progress
  
END PROGRAM gravity
