PROGRAM timer_test
integer i, a
double precision t1, t2

do n= 6,10



call cpu_time(t1)
do i = 1, 100000000
 !      a = 1.2342342/2.4345
         a = sqrt(0.789)       
!         b = 2.4235**3
!         c = sqrt(b) 
enddo 
call cpu_time(t2)

print*, n, t2-t1
enddo 

endprogram timer_test
