MODULE grav_parameters
    ! LEVEL 0 parameters
    INTEGER, PARAMETER :: N_r0=128, N_theta0=256
    REAL(8), DIMENSION(N_r0) :: r0, dr0
    REAL(8), DIMENSION(N_theta0) :: theta0, dtheta0
    REAL(8), DIMENSION(N_r0, N_r0) :: r0_ratio, r0_ratio_squared
    REAL(8), DIMENSION(N_r0, N_theta0) :: mass0
    REAL(8), DIMENSION(1-N_theta0:N_theta0-1) :: cos_table0, sin_table0
    REAL(8), DIMENSION(:), ALLOCATABLE :: r0_squared, r0_prime_squared
    REAL(8), DIMENSION(:, :), ALLOCATABLE :: sigma0
    ! LEVEL X parameters
    REAL(8), DIMENSION(:), ALLOCATABLE :: rx, drx, rx_squared, rx_prime_squared
    REAL(8), DIMENSION(:), ALLOCATABLE :: thetax, dthetax
    REAL(8), DIMENSION(:, :), ALLOCATABLE :: rx_ratio, rx_ratio_squared
    REAL(8), DIMENSION(:, :), ALLOCATABLE :: sigmax, massx
    REAL(8), DIMENSION(:, :), ALLOCATABLE :: cos_tablex, sin_tablex
    REAL(8), DIMENSION(:, :), ALLOCATABLE :: force_r, force_theta
END MODULE grav_parameters
