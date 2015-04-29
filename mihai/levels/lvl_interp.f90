program gravity 
!Variable section
implicit none
INTEGER i, j, k, l

INTEGER, parameter :: nr = 128, nt = 256, lvl = 1

REAL a, G, dr, dth, rad, den ,den_inv, comm, prod, nothing 
double precision t_init, t_end, timed
REAL, dimension(nr) :: radius, radius_corn, radius_corn_2_inv
REAL, dimension(nt) :: theta, cos_center, cos_corner, sin_center, sin_corner
REAL, dimension(nr,nt) :: density, acc_r, acc_t
REAL, dimension(nr,nt) ::  mass

REAL, dimension(nt/2**lvl,nt) :: cos_diff, sin_diff
REAL, dimension(nr/2**lvl,nr) :: ratio, ratio_f
REAL, dimension(1:2) :: acc

INTEGER :: nr_l, nt_l
REAL, dimension(:) , allocatable:: radius_l
REAL, dimension(:) , allocatable:: theta_l
REAL, dimension(:,:), allocatable :: mass_l


REAL, parameter ::  eps=0.0001

call readd(nr, nt, radius, theta, density)

call mass_lvl_0(nr, nt, radius, theta, density, mass)

call corner_r(nr, radius, radius_corn, radius_corn_2_inv)





!IF (lvl == 0) THEN
!   radius_l = radius
!   theta_l  = theta
!   mass_l   = mass
!ELSE IF (lvl == 1) THEN
   call raise_lvl(lvl, nr, nt, radius, theta, mass, radius_l, theta_l, mass_l, nr_l, nt_l)
!ELSE 
!   call raise_lvl(nr, nt, radius, theta, mass, radius_1, theta_1, mass_1, nr_1, nt_1)
!   do i = 1, lvl 
!        i_nr = nr_l
!        i_nt = nt_l
!        raise_lvl(nr, nt, radius, theta, mass, radius_l, theta_l, mass_l, nr_l, nt_l)
!   enddo
!END IF

!call trig_diff(nt, theta, nt_l, theta_l, cos_diff, sin_diff)

!call ration(nr, radius_corn, nr_l, radius_l, ratio, ratio_f)

call cpu_time(t_init)

call gravity_calculator(lvl, nr,nt,nr_l, nt_l,  &
mass_l, radius_corn_2_inv, eps, acc_r, acc_t)

!do i = 1, nr
!        do j = 1, nt
        !acceleration calculations for a single cell
!                acc(1) = 0
!                acc(2) = 0
                
!                do k = 1,nr_l
!                        do l = 1, nt_l
                                !reducing double operations
!                                prod    = ratio(k,i)*2*cos_diff(l,j)

!                                rad     = ratio_f(k,i) - prod + eps 
                                       
!                                den     = rad * sqrt(rad)
!                                comm    = mass_l(k,l)*radius_corn_2_inv(i)/den                                
                                !acc calculations

!                                acc(1) = acc(1) + comm* & 
!                                        (ratio(k,i) - cos_diff(l,j))
                                        
!                                acc(2) = acc(2) + comm* &
 !                                       sin_diff(l,j)
 !                       enddo
 !               enddo
                        
!                acc_r(i,j) = acc(1)
!                acc_t(i,j) = acc(2)
                    
!        enddo
!enddo

call writee(nr, nt, acc_r, acc_t)

call cpu_time(t_end)

OPEN(UNIT = 7, FILE='/home/ics/mihai/git/Computational_Science_II_Open/time_lvl_l.data')

write(7,*) t_end - t_init, t_init, t_end

close(7)


CONTAINS

        SUBROUTINE readd(nr, nt, radius, theta, density)

        IMPLICIT NONE

        integer, intent(in) :: nr, nt
        
        real, dimension(nr), intent(out) :: radius
        real, dimension(nt), intent(out) :: theta
        real, dimension(nr, nt), intent(out) :: density

        real :: a
        integer :: i, j

 
        OPEN(UNIT = 2, FILE ='/home/ics/mihai/git/Computational_Science_II/Data/r_project.data')
        OPEN(UNIT = 3, FILE ='/home/ics/mihai/git/Computational_Science_II/Data/theta_project.data')
        OPEN(UNIT = 4, FILE ='/home/ics/mihai/git/Computational_Science_II/Data/density_project.data')

        do i = 1, nr
                read(2,*) a
                radius(i) = a

        enddo

        do i = 1, nt
                read(3,*) a
                theta(i) = a
        enddo


        do i = 1, nr
                do j = 1, nt
                        read(4,*) a
                        density(i,j) = a
                enddo
        enddo

        close(2)
        close(3)
        close(4)

        END SUBROUTINE readd


        SUBROUTINE writee(nr, nt, acc_r, acc_t)

        IMPLICIT NONE

        integer, intent(in) :: nr, nt

        real, dimension(nr, nt), intent(in) :: acc_r, acc_t
               
        integer :: i, j

        OPEN(UNIT = 5,FILE='/home/ics/mihai/git/Computational_Science_II_Open/acc_r_project.data')
        OPEN(UNIT = 6,FILE='/home/ics/mihai/git/Computational_Science_II_Open/acc_t_project.data')
     
        do i = 1, nr
                do j = 1, nt

                        write(5,*) acc_r(i,j)
                        write(6,*) acc_t(i,j)
                enddo
        enddo

        close(5)
        close(6)

        END SUBROUTINE writee

!PRECALCULATIONS' SUBROUTINES


        SUBROUTINE mass_lvl_0(nr, nt, radius, theta, density, mass)

        IMPLICIT NONE

        integer, intent(in) :: nr, nt 
        real, dimension(nr), intent(in) :: radius
        real, dimension(nt), intent(in) :: theta
        real, dimension(nr, nt), intent(in) :: density


        real, dimension(nr, nt), intent(out) :: mass

        integer :: k, l
        real    :: dr, dt

        dr = radius(2)- radius(1)
        dt = theta(2) - theta(1)

        do k = 1, nr
                do l = 1, nt
                        mass(k,l) = density(k,l)*radius(k)*dr*dt
                enddo
        enddo

        END SUBROUTINE mass_lvl_0


        SUBROUTINE raise_lvl(lvl, nr, nt, radius, theta, mass, radius_l, theta_l, mass_l, nr_l, nt_l) 
        
        IMPLICIT NONE

        integer, intent(in) :: nr, nt, lvl
        real, dimension(nr), intent(in) :: radius
        real, dimension(nt), intent(in) :: theta
        real, dimension(nr, nt), intent(in) :: mass


        integer, intent(out) :: nr_l, nt_l
        real, dimension(:), allocatable, intent(out) :: radius_l
        real, dimension(:), allocatable, intent(out) :: theta_l
        real, dimension(:, :), allocatable, intent(out) :: mass_l
            
        integer :: k, l, j, i, a, b, leng


        leng = 2**lvl
        nr_l = nr/2**lvl
        nt_l = nt/2**lvl

        allocate(mass_l(nr_l, nt_l))
        allocate(radius_l(nr_l))
        allocate(theta_l(nt_l))


        do l = leng, nt, leng
                i = l/leng
                theta_l(i) = 0
                do j = l-leng+1, l 
                        theta_l(i) = theta_l(i) + theta(j)
                enddo
                theta_l(i) = theta_l(i)/leng
        enddo


        do l = leng, nr, leng
                i = l/leng
                radius_l(i) = 0
                do j = l-leng+1, l
                        radius_l(i) = radius_l(i) + radius(j)
                enddo
                radius_l(i) = radius_l(i)/leng
        enddo


        do l = leng, nr, leng
                i = l/leng
                do k = leng, nt, leng
                        j = k/leng
                        mass_l(i,j) = 0

                        do a = l-leng+1, l
                                do b = k-leng+1, k
                                        mass_l(i,j) = mass_l(i,j) + mass(a,b)
                                enddo
                        enddo
                enddo
        enddo
                        

        END SUBROUTINE raise_lvl

!VIRTUAL GRIDS

        SUBROUTINE virtual_grids(lvl, nr_l, nt_l, radius_l, theta_l, mass_l, radius_vg, theta_vg, mass_vg)

        integer, intent(in) :: nt_l, nr_l, lvl
        real, dimension(nr_l), intent(in) :: radius_l
        real, dimension(nt_l), intent(in) :: theta_l
        real, dimension(nr_l,nt_l), intent(in) :: mass_l
        
        real, dimension(:,:,:,:), allocatable, intent(out) :: mass_vg
        real, dimension(:,:,:), allocatable, intent(out) :: radius_vg
        real, dimension(:,:,:), allocatable, intent(out) :: theta_vg        

        integer :: nvr, nvt, aa, bb, i, j, leng
        real :: dr, dr_mic, dt, dt_mic, c11,c22,c12,c21
        nvr = 2**lvl
        nvt = 2**lvl

        leng = 2**lvl

        allocate(mass_vg(nvr, nvt, nr_l+2, nt_l))
        allocate(radius_vg(nvr, nvt, nr_l+2))
        allocate(theta_vg(nvr, nvt, nt_l)) 

        dr = radius_l(2)-radius_l(1)
        dr_mic = dr/leng

        do aa = 1, nvr
                do bb = 1, nvt
                        radius_vg(aa, bb, 1) = radius_l(1)-dr/2 + (aa-1)*dr_mic
                       
                        do i = 2, nr_l+1
                        radius_vg(aa,bb,i) = radius_l(i+1)+(aa-1)*dr_mic
                        enddo
                        radius_vg(aa, bb, nr_l+2) = radius(nr_l) + dr/2 +(aa-1)*dr_mic 
                enddo
        enddo



        dt = theta_l(2)-theta_l(1)
        dt_mic = dt/leng

        do aa = 1, nvr
                do bb = 1, nvt
                        
                        do j = 1, nt_l
                                theta_vg(aa,bb,j) = theta_l(j)+(aa-1)*dt_mic
                        enddo
                        
                enddo
        enddo






        do aa = 1, nvr
                do bb = 1, nvt
                        c11 = (nvr-aa+1)*(nvt-bb+1)
                        c21 = (nvr-aa+1)*(bb-1)
                        c12 = (aa-1)*(nvt-bb+1)
                        c22 = (aa-1)*(bb-1)
                        do i = 2, nr_l+1 
                                do j = 2, nt_l-1
                                        mass_vg(aa, bb, i,j) = c11*mass_l(i-1, j-1) +&
                                                               c12*mass_l(i-1, j) +&
                                                               c21*mass_l(i, j-1) + &
                                                               c22*mass_l(i,j)               
                                enddo
                        enddo                

                        do i = 2, nr_l + 1
                                mass_vg(aa, bb, i,1) = c11*mass_l(i-1,nt_l) +&
                                                               c12*mass_l(i-1,1) +&
                                                               c21*mass_l(i,nt_l) + &
                                                               c22*mass_l(i,1)
                                mass_vg(aa, bb, i,nt_l) = c11*mass_l(i-1,i) +&
                                                               c12*mass_l(i-1,nt_l)+&
                                                               c21*mass_l(i,i)+ &
                                                               c22*mass_l(i,nt_l)   
                                
                        enddo

                        do j = 2, nt_l
                        
                                  mass_vg(aa, bb, 1,j) = c21*mass_l(1,j-1) + &
                                                               c22*mass_l(1,j)
                                  mass_vg(aa, bb, nr_l+2, j) = c11*mass_l(nr_l,j-1) +&
                                                               c12*mass_l(nr_l,j)

                        enddo


                        mass_vg(aa, bb, 1,1) = c11*mass_l(1,1) + c12*mass_l(1,2)
                        mass_vg(aa, bb, nr_l+2, 1) = c11*mass_l(nr_l,1) + c12*mass_l(nr_l, 2)

        enddo
        enddo

        ENDSUBROUTINE virtual_grids


        SUBROUTINE trig_diff(nt, theta, nt_l, theta_l, cos_diff, sin_diff)

        integer, intent(in) :: nt, nt_l
        real, dimension(nt), intent(in) :: theta
        real, dimension(nt_l), intent(in) :: theta_l

        real, dimension(nt_l, nt), intent(out) :: cos_diff, sin_diff

        real :: dt
        integer :: l, j

        dt = theta(3)-theta(2)

        do l = 1, nt_l
                do j = 1, nt
                        cos_diff(l,j) = cos(theta_l(l)-theta(j) + dt/2.)
                        sin_diff(l,j) = sin(theta_l(l)-theta(j) + dt/2.)
                enddo
        enddo


        END SUBROUTINE trig_diff

        SUBROUTINE corner_r(nr, radius, radius_corn, radius_corn_2_inv)

        integer, intent(in) :: nr
        real, dimension(nr), intent(in) :: radius
        real, dimension(nr), intent(out) :: radius_corn, radius_corn_2_inv
        real :: dr
        integer :: i

        dr = radius(3) -radius(2)
        do i = 1, nr
                radius_corn(i) = radius(i) + dr/2
                radius_corn_2_inv(i) = 1./(radius_corn(i)*radius_corn(i))
        enddo

        END SUBROUTINE corner_r


        SUBROUTINE ration(nr, radius_corn, nr_l, radius_l, ratio, ratio_f)

        integer, intent(in) :: nr, nr_l
        real, dimension(nr), intent(in) :: radius_corn
        real, dimension(nr_l), intent(in) :: radius_l
        real, dimension(nr_l, nr), intent(out) :: ratio, ratio_f
        
        integer :: k,i


        do k = 1,nr_l
                do i = 1, nr
                        ratio(k,i) = radius_l(k)/radius_corn(i)
                        ratio_f(k,i) = ratio(k,i)*ratio(k,i) + 1
                enddo
        enddo
        
        END SUBROUTINE ration

!GRAVITY calculator



        SUBROUTINE gravity_calculator(lvl, nr,nt,nr_l, nt_l, mass_l, radius_corn_2_inv, eps, acc_r, acc_t)

        integer, intent(in) :: nr, nt, nr_l, nt_l, lvl
        real, dimension(nr_l,nt_l) , intent(in) :: mass_l
        real, dimension(nr_l, nr) :: ratio, ratio_f
        real, dimension(nt_l, nt)  :: cos_diff
        real, dimension(nr), intent(in) :: radius_corn_2_inv
        
        real, intent(in) :: eps

        real, dimension(nr,nt), intent(out) :: acc_r, acc_t

        real, dimension(:, :, :), allocatable :: radius_vg        
        real, dimension(:, :, :), allocatable:: theta_vg
        real, dimension(:,:, :, :), allocatable :: mass_vg

        real, dimension(nt_l):: theta_l1
        real, dimension(nr_l+2):: radius_l1

        integer:: i, j, aa, bb, k, l

        call virtual_grids(lvl, nr_l, nt_l, radius_l, theta_l, mass_l, radius_vg,theta_vg, mass_vg)




        do aa = 1, 2**lvl

                do bb = 1, 2**lvl
                        
                        !placing the grids
                                do j = 1, nt_l
                                        theta_l1(j) = theta_vg(aa,bb, j)
                                enddo
                                call trig_diff(nt, theta, nt_l, theta_l1, cos_diff, sin_diff)

                                do i = 1, nr_l+2
                                        radius_l1(i) = radius_vg(aa, bb, i)
                                enddo
                                call ration(nr, radius_corn, nr_l+2, radius_l1, ratio, ratio_f)


                                do i = 1, nr_l
                                        do j = 1, nt_l         
                                        acc(1) = 0
                                        acc(2) = 0
                                       


                                do k = 1,nr_l+2
                                      do l = 1, nt_l
                                                 !reducing double operations
                                prod    = ratio(k,(i-1)*2**lvl+aa)*2*cos_diff(l,(j-1)*2**lvl+bb)

                                rad     = ratio_f(k,(i-1)*2**lvl+aa) - prod + eps

                                den     = rad * sqrt(rad)
                                comm    = mass_vg(aa,bb,k,l)*radius_corn_2_inv((i-1)*2**lvl+aa)/den
                                !acc calculations

                                acc(1) = acc(1) + comm* &
                                        (ratio(k,(i-1)*2**lvl+aa) - cos_diff(l,(j-1)*2**lvl+bb))

                                acc(2) = acc(2) + comm* &
                                        sin_diff(l,(j-1)*2**lvl+bb)
                                  enddo
                                enddo

                                acc_r((i-1)*2**lvl+aa,(j-1)*2**lvl+bb) = acc(1)
                                acc_t((i-1)*2**lvl+aa,(j-1)*2**lvl+bb) = acc(2)

                                 enddo
                                enddo

        enddo
        enddo
        




        END SUBROUTINE gravity_calculator







end program gravity
