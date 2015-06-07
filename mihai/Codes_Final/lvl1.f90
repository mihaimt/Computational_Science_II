program gravity 
!Variable section
implicit none
INTEGER i, j, k, l

INTEGER, parameter :: nr = 128, nt = 256, lvl = 1, ng = 2

REAL ag, a, G, dr, dt, rad, den ,den_inv, comm, prod, nothing, c00, c01, c10, c11, massa, radiuss, thetaa 
double precision t_init, t_end, timed
REAL, dimension(nr) :: radius
REAL, dimension(nt) :: theta, cos_center, cos_corner, sin_center, sin_corner
REAL, dimension(nr,nt) :: density, acc_r, acc_t
REAL, dimension(:,:), allocatable ::   acc_rg, acc_tg
REAL, dimension(nr,nt) ::  mass
REAL :: cos_diff, sin_diff
REAL:: ratio, ratio_f
REAL, dimension(1:2) :: acc

INTEGER :: nr_l, nt_l,  nrg
REAL, dimension(:) , allocatable:: radius_l, radius_g
REAL, dimension(:) , allocatable:: theta_l
REAL, dimension(:,:), allocatable :: mass_l, mass_g

REAL, parameter ::  eps=0.0001
REAL theta_ad, radius_ad, valu, radius_corn, radius_corn_2_inv, tl1, tl2
REAL c1, c2, c3, c4
logical don

call readd(nr, nt, radius, theta, density)

call mass_lvl_0(nr, nt, radius, theta, density, mass)

call add_ghost_cells(nr, nt, radius, mass, ng, nrg, radius_g, mass_g)

call raise_lvl(0,0,lvl, nrg, nt, radius_g, theta, mass_g, radius_l, theta_l, mass_l, nr_l, nt_l)

call cpu_time(t_init)


dt =   theta(3)-theta(2)
dr = radius(3) -radius(2)

write(*,*) "Done"


allocate(acc_rg(nrg, nt))
allocate(acc_tg(nrg, nt))

do i = 1, nrg
        do j = 1, nt

                acc_rg(i,j) = 0
                acc_tg(i,j) = 0
        enddo
enddo



!i%2=0, j%2=0  grid
do i = ng, nrg-ng, 2
        do j = 2, nt, 2
        !acceleration calculations for a single cell
                acc(1) = 0
                acc(2) = 0

                radius_corn = radius_g(i)+(radius_g(5)-radius_g(4))/2
                radius_corn_2_inv = 1./(radius_corn*radius_corn)                                

                do k = 1, nr_l

                       radiuss = radius_l(k)
                       ratio = radiuss/radius_corn
                       ratio_f = ratio*ratio + 1

                        do l = 1, nt_l
                 
                                cos_diff = cos(theta_l(l)-theta(j) - dt/2.)
                                sin_diff = sin(theta_l(l)-theta(j) - dt/2.)
                                !reducing double operations
                                prod    = ratio*2*cos_diff

                                rad     = ratio_f - prod + eps

                                den     = rad * sqrt(rad)

                                massa   = mass_l(k,l)
                                         
                                comm    = massa*radius_corn_2_inv/den
                                !acc calculations

                                acc(1) = acc(1) + comm* &
                                        (ratio - cos_diff)

                                acc(2) = acc(2) + comm* &
                                        sin_diff


                        enddo
                enddo 
                acc_rg(i,j) = acc(1)
                acc_tg(i,j) = acc(2)
                    
        enddo
enddo


!i%2=1, j%2=1  grid


c1 = 0.25
c2 = 0.25
c3 = 0.25
c4 = 0.25

do i = ng+1, nrg-ng-1, 2
        do j = 1, nt-1, 2
        !acceleration calculations for a single cell
                acc(1) = 0
                acc(2) = 0

                radius_corn = radius_g(i)+(radius_g(5)-radius_g(4))/2
                radius_corn_2_inv = 1./(radius_corn*radius_corn)
!		write(*,*) i,j, radius_corn, radius_corn_2_inv
                do k = 1, nr_l

                       radiuss = radius_l(k)-dr
                       ratio = radiuss/radius_corn
                       ratio_f = ratio*ratio + 1
!			write(*,*) k, radiuss, ratio, ratio_f
                        do l = 1, nt_l

                                cos_diff = cos(theta_l(l)-dt-theta(j) - dt/2.)
                                sin_diff = sin(theta_l(l)-dt-theta(j) - dt/2.)
!				write(*,*) i,j,k,l,cos_diff, sin_diff
                                !reducing double operations
                                prod    = ratio*2*cos_diff

                                rad     = ratio_f - prod + eps

                                den     = rad * sqrt(rad)

				if (l == 1) then
				massa   = mass_l(k,l)*c1+&
                                          mass_l(k-1,l)*c2+&
                                          mass_l(k,nt_l)*c3+&
                                          mass_l(k-1,nt_l)*c4 
				else


                                massa   = mass_l(k,l)*c1+&
					  mass_l(k-1,l)*c2+&
					  mass_l(k,l-1)*c3+&
					  mass_l(k-1,l-1)*c4


				endif
                                comm    = massa*radius_corn_2_inv/den

!				write(*,*) i,j,k,l,massa
                                !acc calculations

                                acc(1) = acc(1) + comm* &
                                        (ratio - cos_diff)

                                acc(2) = acc(2) + comm* &
                                        sin_diff


                        enddo
                enddo
                acc_rg(i,j) = acc(1)
                acc_tg(i,j) = acc(2)

        enddo
enddo


!i%2=1, j%2=0  grid
c1 = 0.5
c2 = 0.5
c3 = 0
c4 = 0

do i = ng+1, nrg-ng-1, 2
        do j = 2, nt, 2
        !acceleration calculations for a single cell
                acc(1) = 0
                acc(2) = 0

                radius_corn = radius_g(i)+(radius_g(5)-radius_g(4))/2
                radius_corn_2_inv = 1./(radius_corn*radius_corn)
!               write(*,*) i,j, radius_corn, radius_corn_2_inv
                do k = 1, nr_l

                       radiuss = radius_l(k)-dr
                       ratio = radiuss/radius_corn
                       ratio_f = ratio*ratio + 1
!                       write(*,*) k, radiuss, ratio, ratio_f
                        do l = 1, nt_l

                                cos_diff = cos(theta_l(l)-theta(j) - dt/2.)
                                sin_diff = sin(theta_l(l)-theta(j) - dt/2.)
!                               write(*,*) i,j,k,l,cos_diff, sin_diff
                                !reducing double operations
                                prod    = ratio*2*cos_diff

                                rad     = ratio_f - prod + eps

                                den     = rad * sqrt(rad)

                                if (l == 1) then
                                massa   = mass_l(k,l)*c1+&
                                          mass_l(k-1,l)*c2+&
                                          mass_l(k,nt_l)*c3+&
                                          mass_l(k-1,nt_l)*c4
                                else


                                massa   = mass_l(k,l)*c1+&
                                          mass_l(k-1,l)*c2+&
                                          mass_l(k,l-1)*c3+&
                                          mass_l(k-1,l-1)*c4


                                endif
                                comm    = massa*radius_corn_2_inv/den

!                               write(*,*) i,j,k,l,massa
                                !acc calculations

                                acc(1) = acc(1) + comm* &
                                        (ratio - cos_diff)

                                acc(2) = acc(2) + comm* &
                                        sin_diff


                        enddo
                enddo
                acc_rg(i,j) = acc(1)
                acc_tg(i,j) = acc(2)

        enddo
enddo


!i%2=0, j%2=1  grid
c1 = 0.5
c2 = 0
c3 = 0.5
c4 = 0

do i = ng, nrg-ng, 2
        do j = 1, nt-2+1, 2
        !acceleration calculations for a single cell
                acc(1) = 0
                acc(2) = 0

                radius_corn = radius_g(i)+(radius_g(5)-radius_g(4))/2
                radius_corn_2_inv = 1./(radius_corn*radius_corn)
!               write(*,*) i,j, radius_corn, radius_corn_2_inv
                do k = 1, nr_l

                       radiuss = radius_l(k)
                       ratio = radiuss/radius_corn
                       ratio_f = ratio*ratio + 1
!                       write(*,*) k, radiuss, ratio, ratio_f
                        do l = 1, nt_l

                                cos_diff = cos(theta_l(l)-dt-theta(j) - dt/2.)
                                sin_diff = sin(theta_l(l)-dt-theta(j) - dt/2.)
!                               write(*,*) i,j,k,l,cos_diff, sin_diff
                                !reducing double operations
                                prod    = ratio*2*cos_diff

                                rad     = ratio_f - prod + eps

                                den     = rad * sqrt(rad)

                                if (l == 1) then
                                massa   = mass_l(k,l)*c1+&
                                          mass_l(k-1,l)*c2+&
                                          mass_l(k,nt_l)*c3+&
                                          mass_l(k-1,nt_l)*c4
                                else


                                massa   = mass_l(k,l)*c1+&
                                          mass_l(k-1,l)*c2+&
                                          mass_l(k,l-1)*c3+&
                                          mass_l(k-1,l-1)*c4


                                endif
                                comm    = massa*radius_corn_2_inv/den

!                               write(*,*) i,j,k,l,massa
                                !acc calculations

                                acc(1) = acc(1) + comm* &
                                        (ratio - cos_diff)

                                acc(2) = acc(2) + comm* &
                                        sin_diff


                        enddo
                enddo
                acc_rg(i,j) = acc(1)
                acc_tg(i,j) = acc(2)

        enddo
enddo













write(*,*) "Done"

call cpu_time(t_end)

timed = t_end - t_init
write(*,*) "Done"
write(*,*) timed

call writee(nrg, nt, acc_rg, acc_tg)

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

        OPEN(UNIT = 5,FILE='/home/ics/mihai/git/Computational_Science_II_Open/acc_r_noref.data')
        OPEN(UNIT = 6,FILE='/home/ics/mihai/git/Computational_Science_II_Open/acc_t_noref.data')
     
        do i = 3, nr-2
                do j = 1, nt

                        write(5,*) acc_r(i,j)
                        write(6,*) acc_t(i,j)
                enddo
        enddo

        close(5)
        close(6)

        END SUBROUTINE writee


        SUBROUTINE writel(nr, nt, mass_l, radius_l, theta_l)

        IMPLICIT NONE

        integer, intent(in) :: nr, nt

        real, dimension(nr, nt), intent(in) :: mass_l
        real, dimension(nr), intent(in) :: radius_l
        real, dimension(nt), intent(in) :: theta_l

        integer :: i, j

        OPEN(UNIT = 10,FILE='/home/ics/mihai/git/Computational_Science_II_Open/New_attempt/mass_l.data')


        do i = 1, nr
                do j = 1, nt

                        write(10,*) mass_l(i,j)
                       
                enddo
        enddo

        close(10)

        OPEN(UNIT = 11,FILE='/home/ics/mihai/git/Computational_Science_II_Open/New_attempt/radius_l.data')


        do i = 1, nr
                

                        write(11,*) radius_l(i)

                
        enddo

        close(11)


        OPEN(UNIT = 12,FILE='/home/ics/mihai/git/Computational_Science_II_Open/New_attempt/theta_l.data')


        do i = 1, nt
                

                        write(12,*) theta_l(i)

                
        enddo

        close(12)





        END SUBROUTINE writel








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


        SUBROUTINE add_ghost_cells(nr, nt, radius, mass, ng, nrg, radius_g, mass_g)
        
        IMPLICIT NONE
        
        integer, intent(in) :: nr, nt, ng
        real, dimension(nr), intent(in) :: radius
        real, dimension(nr, nt), intent(in) :: mass

        integer, intent(out) :: nrg
        real, dimension(:), allocatable,  intent(out) :: radius_g
        real, dimension(:,:), allocatable,  intent(out)  :: mass_g

        integer :: k, l, i

        nrg= nr + 2*ng

        allocate(radius_g(nrg))
        allocate(mass_g(nrg,nt))



        do k = 1, nrg !+ 2*ng
  
                radius_g(k) = 0.0d0
  
  
                        do l = 1, nt
                                mass_g(k,l) = 0
                        enddo
        enddo



        do k = 1, nr
                i = k+ng
                radius_g(i) = radius(k)
                do l = 1, nt
                        mass_g(i, l) = mass(k,l)

                enddo
        enddo


        do k = 1, ng
                radius_g(k) = radius_g(ng+1) - (ng-k+1)*(radius_g(ng+3)-radius_g(ng+2))
                radius_g(nr+ng+k) = radius_g(nr+ng) + k*(radius_g(ng+3)-radius_g(ng+2))
        enddo
        

 
        
        END SUBROUTINE add_ghost_cells


        SUBROUTINE raise_lvl(di,dj,lvl, nr, nt, radius, theta, mass, radius_l, theta_l, mass_l, nr_l, nt_l) 
        
        IMPLICIT NONE

        integer, intent(in) :: nr, nt, lvl, di, dj
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

        radius_l(1) = radius_l(2)- (radius_l(4) - radius_l(3))
          
        radius_l(nr_l) = radius_l(nr_l-1)  +(radius_l(4) - radius_l(3))


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


end program gravity
