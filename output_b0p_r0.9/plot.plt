set xrange [0:10]
set xtics 1, 1
set grid
#plot "expandedM0.dat" w l, "expandedM1.dat" w l, "expandedM2.dat" w l, "expandedM3.dat" w l, "expandedM4.dat" w l, "expandedM5.dat" w l
splot "expanded3D.dat" w l
