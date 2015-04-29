#!/user/bin/gnuplot

reset
set terminal png
set output 'Radial.png'
set title "Radial Force as a function of radius and azimuth"
set xlabel 'Radius [r0]'
set ylabel 'Azimuth'
set xrange [0.7:1.3]
set yrange [0.0:6.2]
set palette defined ( 0 "#000090",\
                      1 "#000fff",\
                      2 "#0090ff",\
                      3 "#0fffee",\
                      4 "#90ff70",\
                      5 "#ffee00",\
                      6 "#ff7000",\
                      7 "#ee0000",\
                      8 "#7f0000")
#set cblabel "Radial Force"
#set cbrange [-1:1]

plot 'forces.txt' using 1:2:3 w image

set output 'Azimuthal.png'
set title "Azimuthal Force as a function of radius and azimuth"   
set xlabel 'Radius [r0]' 
set ylabel 'Azimuth'
set xrange [0.67:1.33] 
set yrange [0.0:6.28] 
set palette defined ( 0 "#000090",\
                      1 "#000fff",\
                      2 "#0090ff",\
                      3 "#0fffee",\
                      4 "#90ff70",\
                      5 "#ffee00",\
                      6 "#ff7000",\
                      7 "#ee0000",\
                      8 "#7f0000")
#set cbrange [-1:1] 
plot 'forces.txt' using 1:2:4 w image 


