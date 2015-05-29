!-----------------------------------------------------------------------
! MODULE module_grav_mass
!
!     subroutine calc_masslvl0       calculates the mass from density sigma0
!     subroutine calc_masslvlx       calculates the mass from mass0 to a given level N_lvl
!
!-----------------------------------------------------------------------
MODULE grav_mass

    CONTAINS

        SUBROUTINE calc_masslvl0(N_r0, N_theta0, sigma0, r0, dr0, dtheta0, mass0)
            IMPLICIT NONE
            ! DECLARATIONS
            INTEGER, intent(in) :: N_r0, N_theta0
            REAL(8), DIMENSION(:), ALLOCATABLE, intent(in) :: r0, dr0
            REAL(8), DIMENSION(:), ALLOCATABLE, intent(in) :: dtheta0
            REAL(8), DIMENSION(N_r0, N_theta0), intent(in) :: sigma0
            REAL(8), DIMENSION(:,:), ALLOCATABLE, intent(out) :: mass0
            INTEGER :: iprime, jprime
            REAL(8) :: r0dr0_iprime

            allocate(mass0(-1:N_r0+2,-1:N_theta0+2))
            ! fill the <mass table> for LEVEL 0 in the center of the cells
            do iprime = 1, N_r0
                r0dr0_iprime = r0(iprime)*dr0(iprime) ! multiply dr already here to save some operations
                do jprime = 1, N_theta0
                    mass0(iprime,jprime) = -sigma0(iprime,jprime)*r0dr0_iprime*dtheta0(jprime)
                end do
            end do
            ! assign boundary conditions / ghost cells
            mass0(0, :) = 0
            mass0(-1, :) = 0
            mass0(:, 0) = mass0(:, N_theta0)      ! periodic
            mass0(:, -1) = mass0(:, N_theta0-1)   ! periodic
            mass0(N_r0+1, :) = 0
            mass0(N_r0+2, :) = 0
            mass0(:, N_theta0+1) = mass0(:, 1)    ! periodic
            mass0(:, N_theta0+2) = mass0(:, 2)    ! periodic
        END SUBROUTINE calc_masslvl0


        SUBROUTINE calc_masslvlx(N_lvl, N_r0, N_theta0, mass0, masslvl)
            IMPLICIT NONE
            ! DECLARATIONS
            INTEGER, intent(in) :: N_lvl
            INTEGER, intent(in) :: N_r0, N_theta0
            REAL(8), DIMENSION(:,:), ALLOCATABLE, intent(in) :: mass0
            REAL(8), DIMENSION(:,:), ALLOCATABLE, intent(out) :: masslvl

            INTEGER :: i, j, permute_i, permute_j     ! iterators
            INTEGER :: factor_lvl                     ! = 2^N_lvl
            INTEGER :: N_rlvl, N_thetalvl             ! dimensions of parameters in Level N_lvl
            REAL(8) :: temp=0

            factor_lvl=2**N_lvl
            N_rlvl=N_r0/factor_lvl
            N_thetalvl=N_theta0/factor_lvl

            allocate(masslvl(-1:N_rlvl+2, -1:N_thetalvl+2))

            ! fill the <mass table> for LEVEL X
            do i = factor_lvl, N_r0, factor_lvl
                do j = factor_lvl, N_theta0, factor_lvl
                    do permute_i = 0, factor_lvl-1
                        do permute_j = 0, factor_lvl-1
                            temp = temp + mass0(i-permute_i, j-permute_j)
                        end do
                    end do
                masslvl(i/factor_lvl, j/factor_lvl) = temp
                temp = 0
                end do
            end do
            ! assign boundary conditions / ghost cells
            masslvl(0, :) = 0
            masslvl(-1, :) = 0
            masslvl(:, 0) = masslvl(:, N_thetalvl)      ! periodic
            masslvl(:, -1) = masslvl(:, N_thetalvl-1)   ! periodic
            masslvl(N_rlvl+1, :) = 0
            masslvl(N_rlvl+2, :) = 0
            masslvl(:, N_thetalvl+1) = masslvl(:, 1)    ! periodic
            masslvl(:, N_thetalvl+2) = masslvl(:, 2)    ! periodic
        END SUBROUTINE

END MODULE grav_mass
