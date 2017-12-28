reset
set terminal svg noenhanced
set output '/users/staff/ifi/guerrera/Desktop/workspace_eb/wave_FGCS17/experiments/20171222_092118/graph.svg' 
set auto x
set style data histogram
set style histogram errorbars gap 2 lw 1
set style fill solid border -1
set boxwidth 0.9
set yrange [0:30]
set key autotitle columnheader
set xlabel "Number of Threads"
set key title "Implemented Methods"
set ylabel "GFlop/s"
set title "Performance Comparison of Project: wave_FGCS17\nParameters (X_MAX Y_MAX Z_MAX): 200 200 200"
set key top left
set grid
plot '/users/staff/ifi/guerrera/Desktop/workspace_eb/wave_FGCS17/experiments/20171222_092118/results.dat' using 2:3:4:xtic(1) w hist,\
	'' using 5:6:7:xtic(1)  w hist,\
	'' using 8:9:10:xtic(1)  w hist,\
	'' using 11:12:13:xtic(1)  w hist
