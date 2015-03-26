program gravity 
!Variable section
implicit none
INTEGER i, j, k, l, nr, nt
REAL a, G, dr, dth, rad, den ,den_inv, comm, prod
double precision t_init, t_end, timed
REAL, dimension(1:128) :: radius, radius_corn, radius_corn_2_inv
REAL, dimension(1:256) :: theta, cos_center, cos_corner, sin_center, sin_corner
REAL, dimension(1:128,1:256) :: density, acc_r, acc_t, mass
REAL, dimension(1:256,1:256) :: cos_diff
REAL, dimension(1:128,1:128) :: ratio, ratio_f
REAL, dimension(1:4) :: acc
REAL, dimension(1:128,1:128,1:256) :: proj_1, proj_2, proj_3, proj_4
REAL, dimension(1:128,1:256,1:128) :: surf

nr = 128
nt = 256


!Reading Section

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
OPEN(UNIT = 5, FILE='/home/ics/mihai/git/Computational_Science_II_Open/acc_r_project.data')
OPEN(UNIT = 6, FILE='/home/ics/mihai/git/Computational_Science_II_Open/acc_t_project.data')
!Calculations section



call cpu_time(t_init)

G   = 1
dr  = radius(2)-radius(1) 
dth = theta(2) - theta(1)


!Precalculating the trig. functions
do i = 1, nt
        cos_center(i) = cos(theta(i))
        cos_corner(i) = cos(theta(i)+dth/2)
        sin_center(i) = sin(theta(i))
        sin_corner(i) = sin(theta(i)+dth/2)
enddo
do l = 1, nt
        do j = 1, nt
                cos_diff(l,j) = 2*cos(theta(l)-theta(j) + dth/2)
        enddo
enddo 


!Precalculating the corner radii
do i = 1, nr
        radius_corn(i) = radius(i) + dr/2  
        radius_corn_2_inv(i) = 1./(radius_corn(i)*radius_corn(i))      
enddo

!Precalculating the mass 

do k = 1, nr
        do l = 1, nt
                mass(k,l) = density(k,l)*radius(k)*dr*dth
        enddo
enddo

!Precalculating the radius ratio functions

do k = 1,nr
        do i = 1, nr
                ratio(k,i) = radius(k)/radius_corn(i)
                ratio_f(k,i) = ratio(k,i)*ratio(k,i) + 1   
        enddo
enddo

! 3 do precalculation
do k = 1, nr
        do i = 1, nr
                do l = 1, nt
                        surf(k,l,i) = mass(k,l)*radius_corn_2_inv(i)
              
              
                        !projection for 1 -- nt/2
                        proj_1(k,i,l) = ratio(k,i)*cos_center(l)
                        proj_2(k,i,l) = ratio(k,i)*sin_center(l)
                        
                enddo

        enddo
enddo


do i = 1, nr
        do j = 1, nt/2

        !acceleration calculations for a single cell
        !projection Cartesian
                acc(1) = 0
                acc(2) = 0
                acc(3) = 0
                acc(4) = 0
                do k = 1,nr
                        do l = 1, nt
                                !reducing double operations
                                prod    = ratio(k,i)*cos_diff(l,j)

                                rad     = ratio_f(k,i) - prod 
                                       
                                den     = rad * sqrt(rad)
                                comm    = surf(k,l,i)/den                                
                                !acc calculations

                                acc(1) = acc(1) + comm* & 
                                        (proj_1(k,i,l) - & 
                                        cos_corner(j)) 
                                        
                                acc(2) = acc(2) + comm* &
                                        (proj_2(k,i,l) - &
                                        sin_corner(j))
                                !exploiting symmetry in theta along j 
                                rad     = ratio_f(k,i) + prod
                                den     = rad * sqrt(rad)
                                comm    = surf(k,l,i)/den

                                acc(3) = acc(3) + comm* &
                                        (proj_1(k,i,l) + &
                                        cos_corner(j))

                                acc(4) = acc(4) + comm* &
                                        (proj_2(k,i,l) + &
                                        sin_corner(j))
                                        
                        enddo
                enddo
                        
        !projection Polar
                acc_r(i,j) = acc(1)*cos_center(j) + acc(2)*sin_center(j)
                acc_t(i,j) = acc(2)*cos_center(j) - acc(1)*sin_center(j)
        !the other pi

                acc_r(i,j+nt/2) = -acc(3)*cos_center(j) - acc(4)*sin_center(j)
                acc_t(i,j+nt/2) = -acc(4)*cos_center(j) + acc(3)*sin_center(j)

     
        enddo
enddo

do i = 1, nr
        do j = 1, nt

                write(5,*) acc_r(i,j)
                write(6,*) acc_t(i,j)
        enddo
enddo



call cpu_time(t_end)




close(5)
close(6)

OPEN(UNIT = 7, FILE='/home/ics/mihai/git/Computational_Science_II_Open/time_optimiz_2.data')

write(7,*) t_end - t_init, t_init, t_end

close(7)

end program gravity
