program Gravity
  ! Open data files and read to array for r, theta and density
  integer , parameter ::  Nr = 128 , Nth =256
  integer , parameter :: Ndim = Nr * Nth
  real(8), dimension(Nr) :: r
  real(8), dimension(Nth) :: th
  real(8), dimension(Nr,Nth) :: dens
  real(8), dimension(Nr,Nth) :: Mass
  real(8), dimension(Nth,Nth) :: cos_tab, sin_tab

  ! Set parameters and constant for force calculations    
  real(8) :: r0, r02, rp, theta_i, r_prime, Fm, dr, dth, dA, Den, r_squared, cosinerp, density
  real(8) :: Fr0 = 0.0 , Fth0 = 0.0 , FMag = 0.0                                  
  integer :: d0index, tha, thb, radius, angle
  real :: start, finish

  open(unit=1, file="../Data/r_project.data", action="read")
  open(unit=2, file="../Data/theta_project.data", action="read")
  open(unit=3, file="../Data/density_project.data", action="read")
  open(unit=100, file="forces.txt", action="write", status="replace")
  write(100, *) "r                    theta                       F_radial                      F_azimuthal" 
  read(1, *), r
  read(2, *), th

  do i = 1, Nr
     do j = 1, Nth
        read(3, *) density
        dens(i,j) = density
     end do
  end do

  dr  = r(2) -r(1)    ! spacing between grid centers
  dth = th(2)-th(1)
  dA = dr*dth

  ! Fill sin and cos tables
  do tha = 1, Nth
     theta_i = th(tha) - dth/2
     do thb = 1, Nth
        diff_theta = theta_i - th(thb)
        cos_tab(tha,thb) = COS(diff_theta)
        sin_tab(tha,thb) = SIN(diff_theta)
     end do
  end do

  call cpu_time(start)
  ! Calculate masses for each patch
  do radius = 1, Nr
     r_prime = r(radius)
     do angle = 1, Nth
        Mass(radius,angle) = -dA*dens(radius,angle)*r_prime
     end do
  end do

  ! Loop over all points r0 and theta0
  do i = 1, Nr
     r0 = r(i) - dr / 2
     r02 = r0*r0
     do j = 1, Nth
        th0 = th(j) - dth / 2
        ! Loop over all points of r' and theta'
        do k = 1, Nr
           rp = r(k)
           r_squared = r02 + rp*rp
           do m = 1, Nth
              cosinerp = rp*cos_tab(j,m)
              Den = r_squared - 2*r0*cosinerp
              Fm = Mass(k,m)/ (Sqrt(Den)*Den)
              Fth0  = Fth0  + Fm * rp * sin_tab(j,m)   ! azimuthal component
              Fr0   = Fr0 + Fm * ( r0 - cosinerp)      ! radial component
           end do
        end do
        
        ! Write out forces for corner point
        write(100, *) r0,th0,Fr0,Fth0
        Fr0 = 0
        Fth0 = 0
     end do
  end do

  ! Report time used for computation
  call cpu_time(finish)
  print *,"Time (seconds) = ",finish-start

  close(100)
end program Gravity 

! Best time Time (seconds) =    10.2020264  with -O2
