set xrange [0:10]
set xtics 1, 1
set grid lt 0 lw 1
set xlabel font "Helvetica,30"
set xlabel offset 0,0
set xlabel "L"
set bmargin 5
set ylabel font "Helvetica,30"
set ylabel offset 0,0
set ylabel "|ff|"
set lmargin 15
set key font "Helvetica,20"
set tics font "Helvetica,10"
set size 1.0, 1.0 

plot "expandedM0.dat" w l title "M=0" lw 3, \
     "expandedM1.dat" w l title "M=1" lw 3, \
     "expandedM2.dat" w l title "M=2" lw 3, \
     "expandedM3.dat" w l title "M=3" lw 3, \
     "expandedM4.dat" w l title "M=4" lw 3, \
     "expandedM5.dat" w l title "M=5" lw 3
pause -1
