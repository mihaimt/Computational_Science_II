MODULE grav_parameters
    ! LEVEL 0 parameters
    INTEGER, PARAMETER :: N_r0=128, N_theta0=256
    REAL(8), DIMENSION(:), ALLOCATABLE :: r0, dr0
    REAL(8), DIMENSION(:), ALLOCATABLE :: theta0, dtheta0
    REAL(8), DIMENSION(0:N_r0+1, 0:N_r0+1) :: r0_ratio, r0_ratio_squared
    REAL(8), DIMENSION(:, :), ALLOCATABLE :: mass0
    REAL(8), DIMENSION(1-N_theta0:N_theta0-1) :: cos_table0, sin_table0
    REAL(8), DIMENSION(:), ALLOCATABLE :: r0_squared, r0_prime_squared
    REAL(8), DIMENSION(:, :), ALLOCATABLE :: sigma0
    ! FORCE GRID
    REAL(8), DIMENSION(:, :), ALLOCATABLE :: force_r, force_theta
    ! LEVEL X parameters
    REAL(8), DIMENSION(:), ALLOCATABLE :: rx, drx, rx_squared, rx_prime_squared
    REAL(8), DIMENSION(:), ALLOCATABLE :: thetax, dthetax
    REAL(8), DIMENSION(:, :), ALLOCATABLE :: rx_ratio, rx_ratio_squared
    REAL(8), DIMENSION(:, :), ALLOCATABLE :: sigmax, massx
    REAL(8), DIMENSION(:, :), ALLOCATABLE :: cos_tablex, sin_tablex
    ! LEVEL X-1 parameters
    REAL(8), DIMENSION(:), ALLOCATABLE :: rxref, drxref, rxref_squared, rxref_prime_squared
    REAL(8), DIMENSION(:), ALLOCATABLE :: thetaxref, dthetaxref
    REAL(8), DIMENSION(:, :), ALLOCATABLE :: rxref_ratio, rxref_ratio_squared
    REAL(8), DIMENSION(:, :), ALLOCATABLE :: sigmaxref, massxref
    REAL(8), DIMENSION(:, :), ALLOCATABLE :: cos_tablexref, sin_tablexref
END MODULE grav_parameters
