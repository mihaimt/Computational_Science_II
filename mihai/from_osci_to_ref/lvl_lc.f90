program gravity 
!Variable section
implicit none
INTEGER ::i, j, k, l, p, q

INTEGER, parameter :: nr = 128, nt = 256, lvl = 1

REAL :: a, G, dr, dt, rad, den ,den_inv, comm, prod, nothing, c00, c01, c10, c11, massa, radiuss, thetaa 
REAL :: un_ratio, un_cos, un_sin, un_ratio_f

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

REAL, parameter ::  eps=0.001

call readd(nr, nt, radius, theta, density)

call mass_lvl_0(nr, nt, radius, theta, density, mass)

call corner_r(nr, radius, radius_corn, radius_corn_2_inv)

call raise_lvl(lvl, nr, nt, radius, theta, mass, radius_l, theta_l, mass_l, nr_l, nt_l)

call cpu_time(t_init)


dt =   theta(3)-theta(2)
dr = radius(3) -radius(2)


do i = 1, nr
        do j = 1, nt
        !acceleration calculations for a single cell
                acc(1) = 0
                acc(2) = 0


                
                do k = 1, nr_l-1
                        do l = 1, nt_l-1

                                if ((i<2*k+1) .and. (i>2*k-3) .and. (j<2*l+1) .and. (j>2*l-3)) then
                                                                                              
                                do p = 2*k-3, 2*k+1
                                do q = 2*l-3, 2*l+1

                                un_cos = cos(theta(q)-theta(j) + dt/2.)
                                un_sin = sin(theta(q)-theta(j) + dt/2.)

                                radiuss = radius(p)

                                un_ratio = radiuss/radius_corn(i)
                                un_ratio_f = un_ratio*un_ratio + 1




                                !reducing double operations
                                prod    = un_ratio*2*un_cos

                                rad     = un_ratio_f - prod + eps

                                den     = rad * sqrt(rad)

                                massa   = mass(p,q)



                                comm    = massa*radius_corn_2_inv(i)/den
                                !acc calculations

                                acc(1) = acc(1) + comm* &
                                        (un_ratio - un_cos)

                                acc(2) = acc(2) + comm* &
                                        un_sin

                                enddo
                                enddo



                                else

                                if ((MODULO(i,2) ==0) .and. (MODULO(j,2) ==0)) then
                                
                                thetaa = theta_l(l) + dt

                                cos_diff(l,j) = cos(thetaa-theta(j) + dt/2.)
                                sin_diff(l,j) = sin(thetaa-theta(j) + dt/2.)

                                radiuss= radius_l(k)+dr

                                ratio(k,i) = radiuss/radius_corn(i)
                                ratio_f(k,i) = ratio(k,i)*ratio(k,i) + 1




                                !reducing double operations
                                prod    = ratio(k,i)*2*cos_diff(l,j)

                                rad     = ratio_f(k,i) - prod + eps 
                                       
                                den     = rad * sqrt(rad)
                                comm    = mass_l(k,l)*radius_corn_2_inv(i)/den                                
                                !acc calculations

                                acc(1) = acc(1) + comm* & 
                                        (ratio(k,i) - cos_diff(l,j))
                                        
                                acc(2) = acc(2) + comm* &
                                        sin_diff(l,j)

                                endif





                                if ((MODULO(i,2) ==0) .and. (MODULO(j,2) ==1)) then

                                c00 = 0
                                c01 = 0.5
                                c10 = 0
                                c11 = 0.5


                                cos_diff(l,j) = cos(theta_l(l)-theta(j) + dt/2.)
                                sin_diff(l,j) = sin(theta_l(l)-theta(j) + dt/2.)


                                radiuss = radius_l(k)+dr

                                ratio(k,i) = radiuss/radius_corn(i)
                                ratio_f(k,i) = ratio(k,i)*ratio(k,i) + 1




                                !reducing double operations
                                prod    = ratio(k,i)*2*cos_diff(l,j)

                                rad     = ratio_f(k,i) - prod + eps

                                den     = rad * sqrt(rad)


                                massa   = mass_l(k,l)*c00+&
                                          mass_l(k,l+1)*c01+&
                                          mass_l(k+1,l)*c10+&
                                          mass_l(k+1,l+1)*c11




                                comm    = massa*radius_corn_2_inv(i)/den
                                !acc calculations

                                acc(1) = acc(1) + comm* &
                                        (ratio(k,i) - cos_diff(l,j))

                                acc(2) = acc(2) + comm* &
                                        sin_diff(l,j)

                                endif



                                if ((MODULO(i,2) ==1) .and. (MODULO(j,2) ==0)) then

                                c00 = 0.5
                                c01 = 0
                                c10 = 0.5
                                c11 = 0

                                thetaa = theta_l(l) + dt

                                cos_diff(l,j) = cos(thetaa-theta(j) + dt/2.)
                                sin_diff(l,j) = sin(thetaa-theta(j) + dt/2.)

                                radiuss = radius_l(k)

                                ratio(k,i) = radiuss/radius_corn(i)
                                ratio_f(k,i) = ratio(k,i)*ratio(k,i) + 1




                                !reducing double operations
                                prod    = ratio(k,i)*2*cos_diff(l,j)

                                rad     = ratio_f(k,i) - prod + eps

                                den     = rad * sqrt(rad)

                                massa   = mass_l(k,l)*c00+&
                                          mass_l(k,l+1)*c01+&
                                          mass_l(k+1,l)*c10+&
                                          mass_l(k+1,l+1)*c11




                                comm    = massa*radius_corn_2_inv(i)/den
                                !acc calculations

                                acc(1) = acc(1) + comm* &
                                        (ratio(k,i) - cos_diff(l,j))

                                acc(2) = acc(2) + comm* &
                                        sin_diff(l,j)

                                endif




                                if ((MODULO(i,2) ==1) .and. (MODULO(j,2) ==1))then
                                c00 = 0.25
                                c01 = 0.25
                                c10 = 0.25
                                c11 = 0.25



                                cos_diff(l,j) = cos(theta_l(l)-theta(j) + dt/2.)
                                sin_diff(l,j) = sin(theta_l(l)-theta(j) + dt/2.)

                                radiuss = radius_l(k)

                                ratio(k,i) = radiuss/radius_corn(i)
                                ratio_f(k,i) = ratio(k,i)*ratio(k,i) + 1




                                !reducing double operations
                                prod    = ratio(k,i)*2*cos_diff(l,j)

                                rad     = ratio_f(k,i) - prod + eps

                                den     = rad * sqrt(rad)

                                massa   = mass_l(k,l)*c00+&
                                          mass_l(k,l+1)*c01+&
                                          mass_l(k+1,l)*c10+&
                                          mass_l(k+1,l+1)*c11



                                comm    = massa*radius_corn_2_inv(i)/den
                                !acc calculations

                                acc(1) = acc(1) + comm* &
                                        (ratio(k,i) - cos_diff(l,j))

                                acc(2) = acc(2) + comm* &
                                        sin_diff(l,j)

                                endif

                                !level 0 part
                                
                              
                                endif













                        enddo
                enddo
                        
                acc_r(i,j) = acc(1)
                acc_t(i,j) = acc(2)
                    
        enddo
enddo

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

end program gravity
