#!/usr/bin/gnuplot

set autoscale
set grid
set xtic auto
set ytic auto

# Plotting one by one
set terminal png size 800,600 enhanced font "Helvetica,14"

set xlabel "Time (s)"
set ylabel "CO2_g"
set output "CO2_g.png"
plot 'CO2_g.txt' w l lt 1 lc rgb "red" lw 2.1 title 'CO2_g'

set xlabel "Time (s)"
set ylabel "CO2_l"
set output "CO2_l.png"
plot 'CO2_l.txt' w l lt 1 lc rgb "red" lw 2.1 title 'CO2_l'

set xlabel "Time (s)"
set ylabel "L_g"
set output "L_g.png"
plot 'L_g.txt' w l lt 1 lc rgb "green" lw 2.1 title 'L_g'

set xlabel "Time (s)"
set ylabel "L_l"
set output "L_l.png"
plot 'L_l.txt' w l lt 1 lc rgb "green" lw 2.1 title 'L_l'

set xlabel "Time (s)"
set ylabel "m0"
set output "m0.png"
plot 'm0.txt' w l lt 1 lc rgb "black" lw 2.1 title 'm0'

set xlabel "Time (s)"
set ylabel "m1"
set output "m1.png"
plot 'm1.txt' w l lt 1 lc rgb "black" lw 2.1 title 'm1'

set xlabel "Time (s)"
set ylabel "m2"
set output "m2.png"
plot 'm2.txt' w l lt 1 lc rgb "black" lw 2.1 title 'm2'

set xlabel "Time (s)"
set ylabel "m3"
set output "m3.png"
plot 'm3.txt' w l lt 1 lc rgb "black" lw 2.1 title 'm3'

set xlabel "Time (s)"
set ylabel "M0"
set output "M0.png"
plot 'M0.txt' w l lt 1 lc rgb "black" lw 2.1 title 'M0'

set xlabel "Time (s)"
set ylabel "M1"
set output "M1.png"
plot 'M1.txt' w l lt 1 lc rgb "black" lw 2.1 title 'M1'

set xlabel "Time (s)"
set ylabel "M2"
set output "M2.png"
plot 'M2.txt' w l lt 1 lc rgb "black" lw 2.1 title 'M2'

set xlabel "Time (s)"
set ylabel "M3"
set output "M3.png"
plot 'M3.txt' w l lt 1 lc rgb "black" lw 2.1 title 'M3'

set xlabel "Time (s)"
set ylabel "rho-bubble"
set output "rho_bubble.png"
plot 'rho_bubble.txt' w l lt 1 lc rgb "yellow" lw 2.1 title 'rho-bubble'

set xlabel "Time (s)"
set ylabel "rho-foam"
set output "rho_foam.png"
plot 'rho_foam.txt' w l lt 1 lc rgb "black" lw 2.1 title 'rho-foam'

set xlabel "Time (s)"
set ylabel "T, (K)"
set output "T.png"
plot 'T.txt' w l lt 1 lc rgb "black" lw 2.1 title 'T'

set xlabel "Time (s)"
set ylabel "XOH"
set output "XOH.png"
plot 'XOH.txt' w l lt 1 lc rgb "black" lw 2.1 title 'XOH'

set xlabel "Time (s)"
set ylabel "XW"
set output "XW.png"
plot 'XW.txt' w l lt 1 lc rgb "black" lw 2.1 title 'XW'

pause 1  

# This will re-read the same files and can be used on the fly
# reread

# lt: linetype, lc: linecolor, ls:linestyle, lw: linewidth, w l: with line
# linetype: -1=black 1=red 2=grn 3=blue 4=purple 5=aqua 6=brn 7=orange 8=light-brn
# for postscipt -1=normal, 1=grey, 2=dashed, 3=hashed, 4=dot, 5=dot-dash
