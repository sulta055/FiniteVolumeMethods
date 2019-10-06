#set terminal png

set terminal gif animate delay 10 size 1280, 680
set output 'eno-weno.gif'

set key font ",10"

nt=400
nx=100

do for [i=0:nt-1] {
  set title "Time Step = #".(i+1)
  #set title "Piecewise Polynomial Reconstruction"
  set xlabel "x"
  set ylabel "v(x)"
  set yrange [-1:1]
  set pointsize 0.6
  plot "output_adv.txt" every ::i*nx+1::nx+(i*nx)-1 using 1:2 with points pointtype 4 lc rgb "red" title "WENO5",\
 "output_vl.txt" every ::i*nx+1::nx+(i*nx)-1 using 1:2 with linespoints pointtype 4 lc rgb "blue" title "Van Leer"

}


unset output
