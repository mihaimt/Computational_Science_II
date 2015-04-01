program ex1

  implicit none             ! all variables must be defined
  real(8) t_init, t_end     ! global timing
  integer, parameter :: N_r = 128
  integer, parameter :: N_theta = 256
  integer :: theta_i
  real(8) :: r_prime(0:N_r-1), theta_prime(0:N_theta-1)
  real(8),dimension(0:(N_r*N_theta)-1) :: sigma
  real(8) :: theta_step
  real(8) :: cosCache(0:N_theta-1), sinCache(0:N_theta-1)

  
  write (*,*) 'Program Start ', ' ... '
  write (*,*) 'Data read ...'
  call CPU_TIME(t_init)
  call readFile("r_project.data", r_prime)
  call readFile("theta_project.data", theta_prime)
  call readFile("density_project.data", sigma)
  call CPU_TIME(t_end)
  write (*,*) 'Data read done.'
  print *, "Time: ", t_end - t_init , "s"

  theta_step = theta_prime(1) - theta_prime(0)
  ! pre-compute the cache for all theta differences that occur, for N_theta/4 values of N_theta
  do theta_i=0, (N_theta/4)-1
    cosCache(theta_i) = cos((theta_i+0.5)*theta_step)
  end do


  call testCpuTime(1000000)
  call testCosine(1000000, 500d0)
  call testSine(1000000, 500d0)
  call testMultiplications(1000000,2)
  call testPowers(1000000,2.0)
  call testMultiplications(1000000,3)
  call testPowers(1000000,3.0)
  call testMultiplications(1000000,4)
  call testPowers(1000000,4.0)
  call testSqrt(1000000, 100000d0)
  call testPowers(1000000,0.5)
  call testIfStatements(10000000)

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
