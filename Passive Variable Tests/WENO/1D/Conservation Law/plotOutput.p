#set terminal png

set terminal gif animate delay 10 size 1280, 680
set output 'eno-weno.gif'

set key font ",10"

nx=100
nt=200

do for [i=0:nt-1] {
  set title "Time Step = #".(i+1)
  #set title "Piecewise Polynomial Reconstruction"
  set xlabel "x"
  set ylabel "v(x)"
  set pointsize 0.6
  plot "output_consv.txt" every ::i*nx+1::nx+(i*nx)-1 using 1:2 with linespoints pointtype 6 lc rgb "blue" title "WENO5" 
}


unset output
