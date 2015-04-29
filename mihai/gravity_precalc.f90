program gravity 
!Variable section
implicit none
INTEGER i, j, k, l, nr, nt
REAL a, G, dr, dth
double precision t_init, t_end, timed
REAL, dimension(1:128) :: radius, radius_corn
REAL, dimension(1:256) :: theta, cos_center, cos_corner, sin_center, sin_corner
REAL, dimension(1:128,1:256) :: density, acc_r, acc_t, mass
REAL, dimension(1:256,1:256) :: cos_diff
REAL, dimension(1:2) :: acc


!Reading Section

OPEN(UNIT = 2, FILE ='/home/ics/mihai/git/Computational_Science_II/Data/r_project.data')
OPEN(UNIT = 3, FILE ='/home/ics/mihai/git/Computational_Science_II/Data/theta_project.data')
OPEN(UNIT = 4, FILE ='/home/ics/mihai/git/Computational_Science_II/Data/density_project.data')



do i = 1, 128
        read(2,*) a
        radius(i) = a

enddo




do i = 1, 256
        read(3,*) a
        theta(i) = a
enddo


do i = 1, 128
        do j = 1, 256
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

nr = 128
nt = 256


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
                cos_diff(l,j) = cos(theta(l)-theta(j) + dth/2)
        enddo
enddo 


!Precalculating the corner radii
do i = 1, nr
        radius_corn(i) = radius(i) + dr/2        
enddo

!Precalculating the mass 

do k = 1, nr
        do l = 1, nt
                mass(k,l) = density(k,l)*radius(k)*dr*dth
        enddo
enddo


do i = 1, nr
        do j = 1, nt

        !acceleration calculations for a single cell
        !projection Cartesian
          acc(1) = 0
          acc(2) = 0
          do k = 1,nr
             do l = 1, nt
                acc(1) = acc(1) + mass(k,l)* & 
                        (radius(k)*cos_center(l) - & 
                        radius_corn(i)*cos_corner(j)) & 
                        /(radius(k)**2 + radius_corn(i)**2 -2*radius(k) & 
                        *radius_corn(i)*cos_diff(l,j))**(1.5) 
         
               acc(2) = acc(2) + mass(k,l)* &
                        (radius(k)*sin_center(l) - &
                         radius_corn(i)*sin_corner(j)) &
                         /(radius(k)**2 + radius_corn(i)**2 -2*radius(k) &
                         *radius_corn(i)*cos_diff(l,j))**(1.5)
              enddo
           enddo
                        
        !projection Polar
           acc_r(i,j) = acc(1)*cos_center(j) + acc(2)*sin_center(j)
           acc_t(i,j) = acc(2)*cos_center(j) - acc(1)*sin_center(j)
           write(5,*) acc_r(i,j)
           write(6,*) acc_t(i,j)
        enddo
enddo

call cpu_time(t_end)




close(5)
close(6)

OPEN(UNIT = 7, FILE='/home/ics/mihai/git/Computational_Science_II_Open/time_precalc.data')
write(7,*) t_end - t_init, t_init, t_end

close(7)

end program gravity
