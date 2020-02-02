set xrange [0:10]
set xtics 1, 1
set grid
set xlabel font "Helvetica,30"
set xlabel offset 0,0
set xlabel "L"
#set bmargin 5
set ylabel font "Helvetica,30"
set ylabel offset 0,0
set ylabel "M"
#set lmargin 15
set zlabel font "Helvetica,30"
set zlabel offset 0,0
#set zlabel "|ff|"
set nokey
set tics font "Helvetica,10"

splot "expanded3D.dat" w l lw 2
pause -1
