program gravity 
!Variable section
implicit none
INTEGER i, j, k, l
REAL a, G, dr, dth
double precision t_init, t_end, timed
REAL, dimension(1:128) :: radius
REAL, dimension(1:256) :: theta
REAL, dimension(1:128,1:256) :: density, acc_r, acc_t
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
OPEN(UNIT = 5, FILE='/home/ics/mihai/git/Computational_Science_II_Open/force_r_project.data')
OPEN(UNIT = 6, FILE='/home/ics/mihai/git/Computational_Science_II_Open/force_t_project.data')
!Calculations section

call cpu_time(t_init)

G   = 1
dr  = radius(2)-radius(1) 
dth = theta(2) - theta(1)
do i = 1, Nr
      do j = 1, Nt

        
          acc(1) = 0
          acc(2) = 0
          do k = 1,Nr
             do l = 1, Nt
               acc(1) = acc(1) + density(k,l)*dr*dth*radius(k)* &
                        (radius(i) +dr/2 -radius(k)*cos(theta(j)-theta(l)+dth/2))&
                        /(radius(k)**2 + (radius(i)+dr/2)**2 -2*radius(k) & 
                        *(radius(i)+dr/2)*cos(theta(l)-theta(j)+dth/2))**(1.5)
              
               acc(2) = acc(2) + density(k,l)*dr*dth*radius(k)* &
                        sin(theta(j)-theta(l)+dth/2) &
                        /(radius(k)**2 + (radius(i)+dr/2)**2 -2*radius(k) &
                        *(radius(i)+dr/2)*cos(theta(l)-theta(j)+dth/2))**(1.5)
             enddo
           enddo
                        
        
           acc_r(i,j) = acc(1)
           acc_t(i,j) = acc(2)
           write(5,*) acc_r(i,j)
           write(6,*) acc_t(i,j)
        enddo
enddo

call cpu_time(t_end)




close(5)
close(6)

OPEN(UNIT = 7, FILE='/home/ics/mihai/git/Computational_Science_II_Open/time.data')
write(7,*) t_end - t_init, t_init, t_end

close(7)

end program gravity
