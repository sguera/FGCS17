reset
set terminal svg noenhanced
set output '/users/staff/math/maffia/workspace_prova/Jacobi2D/experiments/20171222_111857/graph.svg' 
set auto x
set style data histogram
set style histogram errorbars gap 2 lw 1
set style fill solid border -1
set boxwidth 0.9
set yrange [0:*]
set key autotitle columnheader
set xlabel "Implemented Methods"
set key title "Parameters (PAR_I PAR_J PAR_ITER)"
set ylabel "Execution time (ms)"
set title "Performance Comparison of Project: Jacobi2D\nNumber of Threads: 8"
set key top left
set grid
plot '/users/staff/math/maffia/workspace_prova/Jacobi2D/experiments/20171222_111857/results.dat' using 2:3:4:xtic(1) w hist,\
	'' using 5:6:7:xtic(1)  w hist,\
	'' using 8:9:10:xtic(1)  w hist,\
	'' using 11:12:13:xtic(1)  w hist
