program ex1
  !  optimised 
  !  using cosine and sin caching of variables without symmetry
  !  caching any available values outside the loops
  !  f95 -O3 ex7.f95 -o ex7 && ./ex7
  !  5 trials
  !      Time:    12.281551000000000      s
  !      Time:    12.228482000000000      s
  !      Time:    12.225080000000000      s
  !      Time:    12.246321000000000      s
  !      Time:    12.264409000000001      s

  ! implicit none // all variables must be defined
  real(8) t_init, t_end
  real, parameter :: G = 1 !normalised to 1; 6.67384e-11 !(m*m*m)/(kg*s*s)
  real, parameter :: pi = 3.1415926538
  integer, parameter :: N_r = 128
  integer, parameter :: N_theta = 256
  real(8) :: epsilonCos = 1d-8
  real(8) :: a,r,r2, theta, r_step, theta_step, r_p, theta_p, force_mag_temp, surface_area, d_rd_theta
  ! use zero-based indexing for arrays to make modulo math work better below.
  real(8),dimension(0:(N_r*N_theta)-1) :: force_r, force_theta, force_mag, sigma
  real(8) :: r_prime(0:N_r-1), theta_prime(0:N_theta-1)
  real(8) :: u_r, u_theta ! the unit vectors
  integer :: out_i, in_i, theta_i
  real(8) :: cosvar, sinvar, cosCache(0:N_theta-1), sinCache(0:N_theta-1)
  real(8) :: i_table, t_table, i_force, t_force, i_trig, t_trig
  
  write (*,*) 'Program Start ', ' ... '
  write (*,*) 'Data read ...'
  call readFile("r_project.data", r_prime)
  call readFile("theta_project.data", theta_prime)
  call readFile("density_project.data", sigma)
  write (*,*) 'Data read done.'
   
  call CPU_TIME(t_init)
  r_step = r_prime(1) - r_prime(0)              !  
  theta_step = theta_prime(1) - theta_prime(0)  !

  !Old and new cosine Caching
  !i_table = start() 
  do theta_i=0, (N_theta)-1
    cosCache(theta_i) = cos((theta_i+0.5)*theta_step) 
    sinCache(theta_i) = sin((theta_i+0.5)*theta_step)
  end do

  !do theta_i=0, (N_theta/4)-1
  !  cosCache(theta_i) = cos((theta_i+0.5)*theta_step) 
  !end do

  !old and new cosine caching timers
  !t_table = stop(i_table)
  !print *, "Table Time: ", t_table, "s"

  ! Checking for errors within our epsilon threshold
  !do theta_i=0, N_theta-1
  !  IF (abs(cosCache(theta_i) - cosFast(theta_i)) > epsilonCos) THEN
  !     print *, "ERROR COS, " , theta_i, cosCache(theta_i), cosFast(theta_i), "..."
  !  END IF
  !  IF (abs(sinCache(theta_i) - sinFast(theta_i)) > epsilonCos) THEN
  !     print *, "ERROR SIN, " , theta_i, sinCache(theta_i), sinFast(theta_i), abs(sinCache(theta_i)-sinFast(theta_i)), "..."
  !  END IF
  !end do

  !t_force = 0
  !t_trig = 0
  
  d_rd_theta = r_step * theta_step
  do out_i=0, (N_r*N_theta)-1           ! 
  !do out_i=0, 1000
    r=r_prime(out_i/N_theta)-(r_step/2.0)
    r2=r*r !putting into memory causes slowdown
    theta = theta_prime(MODULO(out_i, N_theta))-(theta_step/2.0)
    force_r(out_i)=0
    force_theta(out_i)=0
    do in_i=0, (N_r*N_theta)-1 
      !i_trig = start()
      r_p = r_prime(in_i/N_theta)               !current r-prime
      theta_p = theta_prime(MODULO(in_i, N_theta))  !current theta-prime
      !theta_p = theta_prime(IAND(in_i, N_theta-1)) !no significant change
      !surface_area = theta_step*r_p*r_step
      !d_rd_theta = r_step * theta_step ! slower outside loop???!?!?WHY?!?
      !cosvar = cos(theta-theta_p)
      !sinvar = sin(theta-theta_p)
      !cosvar = cosCache(MODULO(in_i-out_i, N_theta))
      !sinvar = sinCache(MODULO(in_i-out_i, N_theta))
      !cosvar = cosFast(in_i-out_i)
      !sinvar = sinFast(in_i-out_i)
      cosvar = cosCache(IAND(in_i-out_i, N_theta-1))  ! no significant change
      sinvar = sinCache(IAND(in_i-out_i, N_theta-1))
      !t_trig = t_trig + stop(i_trig)

      !error margins Check
      !IF (abs(cosvar - cosCache( MODULO(in_i-out_i, N_theta) ) ) > 1d-8) THE
      !  write(*,*) 'Does not equal', out_i, in_i, cosvar, cosCache( MODULO(in_i-out_i, N_theta))
      !  stop
      !END IF
      
      !i_force = start()
      force_denominator = (r2+r_p*r_p-(2*r*r_p*cosvar))
      !force_denominator = sqrt(force_denominator)*force_denominator     !pow(2/3)
      force_denominator = sqrt(force_denominator*force_denominator*force_denominator)
      force_mag_temp = (G*sigma(in_i)/(force_denominator)) * d_rd_theta
      force_r(out_i) = (r - r_p*cosvar) * force_mag_temp + force_r(out_i)
      force_theta(out_i) = (r_p*sinvar) * force_mag_temp + force_theta(out_i)
      !t_force = t_force + stop(i_force)
    end do
    force_mag(out_i) = sqrt(force_r(out_i)*force_r(out_i) + force_theta(out_i)*force_theta(out_i)) 
  end do

  call CPU_TIME(t_end)

  call writeFile("force_r.data",force_r)
  call writeFile("force_theta.data",force_theta)
  call writeFile("force_mag.data", force_mag)

  print *, "Time: ", t_end - t_init , "s"
  !print *, "Trigonometry Time: ", t_trig, "s"
  !print *, "Force Time: ", t_force , "s"

  !call testCpuTime(100000)
  !call testCosine(100000, 500d0)
  !call testSine(100000, 500d0)
  !call testMultiplications(100000,2)
  !call testPowers(100000,2.0)
  !call testMultiplications(100000,3)
  !call testPowers(100000,3.0)
  !call testMultiplications(100000,4)
  !call testPowers(100000,4.0)
  !call testSqrt(100000, 100000d0)
  !call testPowers(100000,0.5)
  !call testIfStatements(10000)

contains 

  real(8) function cosFast(x)
    integer :: x, y
    y = MODULO(x, N_theta)
    IF (y < N_theta/4) THEN
        cosFast = cosCache(y)
    ELSE IF (y < N_theta/2) THEN
        cosFast = -cosCache((N_theta/2-y) -1)
    ELSE IF (y < 3*N_theta/4) THEN
        cosFast = -cosCache(y - N_theta/2 )
    ELSE 
        cosFast = cosCache(N_theta-y -1)
    END IF
  end function cosFast

  real(8) function sinFast(x)
    integer :: x
    sinFast = cosFast(x - N_theta/4)
  end function sinFast

  real(8) function start()
    call CPU_TIME(start)
    return
  end function start

  real(8) function stop(started)
    real(8) :: started
    real(8) :: current
    call CPU_TIME(current)
    stop = current - started
    return
  end function stop

  subroutine testIfStatements(n)
    integer :: i, n, x
    
    call CPU_TIME (t_init)
    do i=1,n
      IF (n > i) THEN
        x=i
      ELSE IF (n < i) THEN
        x=n
      END IF
    end do
    call CPU_TIME(t_end)
    print *, "If-Statement Test - Time: ", (t_end - t_init)/n , "s"

  end subroutine

  subroutine testCpuTime(n)
    integer :: i, n
    real(8) :: x, y, t_init, t_end
    y = 0
    call CPU_TIME (t_init)
    do i=1,n
        x = start()
        y = y + stop(x)
        !call CPU_TIME(x)
    end do
    call CPU_TIME(t_end)
    print *, "CpuTime Test - Time: ", (t_end - t_init)/n , "s"
  end subroutine

  subroutine testCosine(n, y)
    integer :: i, n
    real(8) :: x, y, t_init, t_end
    call CPU_TIME (t_init)
    do i=1,n
        x=cos(y)
    end do
    call CPU_TIME(t_end)
    print *, "Cosine Test - Time: ", (t_end - t_init)/n , "s"
    call CPU_TIME (t_init)
    do i=1,n
        x=cosFast(i)
    end do
    call CPU_TIME(t_end)
    print *, "CosFast Test - Time: ", (t_end - t_init)/n , "s"
  end subroutine

  subroutine testSine(n, y)
    integer :: i, n
    real(8) :: x, y, t_init, t_end
    call CPU_TIME (t_init)
    do i=1,n
        x=sin(y)
    end do
    call CPU_TIME(t_end)
    print *, "Sine Test - Time: ", (t_end - t_init)/n , "s"
    call CPU_TIME (t_init)
    do i=1,n
        x=sinFast(i)
    end do
    call CPU_TIME(t_end)
    print *, "SinFast Test - Time: ", (t_end - t_init)/n , "s"
  end subroutine

  subroutine testMultiplications(n, m)
    integer :: i, n, m
    real(8) :: x, t_init, t_end
    IF (m==2) THEN
      call CPU_TIME (t_init)
      do i=1,n
        x=n*n
      end do
      call CPU_TIME(t_end)
    ELSE IF (m==3) THEN
      call CPU_TIME (t_init)
      do i=1,n
        x=n*n*n
      end do
      call CPU_TIME(t_end)
    ELSE IF (m==4) THEN
      call CPU_TIME (t_init)
      do i=1,n
        x=n*n*n*n
      end do
      call CPU_TIME(t_end)
    ELSE IF (m==5) THEN
      call CPU_TIME (t_init)
        do i=1,n
          x=n*n*n*n
        end do
      call CPU_TIME(t_end)
    ELSE
      print *, "ERROR for Multiplication test"
    END IF 
    print *, "Multiplications ", m, "* - Time: ", (t_end - t_init)/n , "s"
  end subroutine

  subroutine testSqrt(n, m)
    integer :: i, n
    real(8) :: x, m
    call CPU_TIME (t_init)
    do i=1,n
      x=sqrt(m)
    end do
    call CPU_TIME(t_end)
    print *, "Square Root * - Time: ", (t_end - t_init)/n , "s"
  end subroutine

  subroutine testPowers(n, m)
    integer :: i, n
    real(4) :: m
    real(8) :: x, t_init, t_end
    call CPU_TIME (t_init)
      do i=1,n
        x=n**m
      end do
    call CPU_TIME(t_end)
    print *, "Powers of ", n, " to the ",m, "* - Time: ", (t_end - t_init)/n , "s"
  end subroutine




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


end program ex1
