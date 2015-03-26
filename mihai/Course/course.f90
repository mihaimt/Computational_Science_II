program gravity 
!Variable section
implicit none
INTEGER i, j, k, l, nr, nt
REAL a, G, dr, dth, rad, den ,den_inv, comm, prod
double precision t_init, t_end, timed
REAL, dimension(1:128) :: radius, radius_corn, radius_corn_2_inv
REAL, dimension(1:64) :: radius_c
REAL, dimension(1:256) :: theta, cos_center, cos_corner, sin_center, sin_corner
REAL, dimension(1:128) :: theta_c
REAL, dimension(1:128,1:256) :: density, acc_r, acc_t, mass
REAL, dimension(1:64,1:128) :: density_c
REAL, dimension(1:256,1:256) :: cos_diff, sin_diff
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


OPEN(UNIT = 5, FILE ='/home/ics/mihai/git/Computational_Science_II/Data/r_project_c2.data')
OPEN(UNIT = 6, FILE ='/home/ics/mihai/git/Computational_Science_II/Data/theta_project_c2.data')
OPEN(UNIT = 7, FILE ='/home/ics/mihai/git/Computational_Science_II/Data/density_project_c2.data')


do i = 1, nt/2
        theta_c(i) = (theta(2*i-1) + theta(2*i))/2.
        write(6,*) theta_c(i)
enddo 


do i = 1, nr/2
        radius_c(i) = (radius(2*i-1)+radius(2*i))/2.
        write(5,*) radius_c(i)
        do j = 1, nt/2
                density_c(i,j) = (density(2*i-1, 2*j-1)+ density(2*i-1, 2*j) + density(2*i, 2*j-1) + density(2*i, 2*j))/4.
                write(7,*) density_c(i,j)
        enddo
enddo


close(5)
close(6)
close(7)





end program gravity
