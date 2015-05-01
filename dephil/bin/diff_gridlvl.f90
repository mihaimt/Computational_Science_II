!--------------------------------------------------------------------------
! PROGRAM TO GET THE DIFFERENCES OF THE DATA 
!        FROM grav_force_lvlx, oscillation_force, oscillation_mass
!    compile by hand...
!
!    use:      ./diff_grid_lvl [FILENAME1].data [FILENAME2].data
!--------------------------------------------------------------------------
PROGRAM diff_gridlvl
    IMPLICIT NONE
    ! DECLARATIONS
    INTEGER :: N_r0=128, N_theta0=256    ! dimensions of the grid
    INTEGER :: i                         ! iterators
    REAL(8) :: f_1, f_2
    CHARACTER(len=50) :: arg1, arg2, arg3
    
    ! GET FILES FROM TERMINAL LINE
    call getarg(1, arg1)
    call getarg(2, arg2)
    call getarg(3, arg3)

    ! OPEN FILES
    open(unit=8, file=arg1, action='read')
    open(unit=9, file=arg2, action='read')
    open(unit=10, file=arg3, action='write')

    do i=1, N_r0*N_theta0
        read(8, '(e20.10)') f_1
        read(9, '(e20.10)') f_2
        write(10, '(e20.10)') f_2-f_1
    end do

    close(unit=8)
    close(unit=9)
    close(unit=10)

END PROGRAM diff_gridlvl
