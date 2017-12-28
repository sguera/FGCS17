#
# Makefile for Patus stencil benchmark
#
# Note: $(PATUS_*) variables will be automatically replaced by the
# required runtime files by Patus.
#
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S), Darwin)
	CC = /usr/local/bin/gcc
	GNUPLOT = /Applications/Gnuplot.app/Contents/Resources/bin/gnuplot
else
	CC = gcc
	GNUPLOT = gnuplot
endif

CFLAGS = -O3 -g -msse2 -fopenmp -I.


bench: kernel.o driver.o ssemathlib.o timer.o util.o
	$(CC) $(CFLAGS) -o $@ $+ -lm

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<
	
tune:
ifndef X_MAX
	@echo For tuning, the value of X_MAX must be specified. E.g. 'make tune X_MAX=100'
	@false
endif
ifndef Y_MAX
	@echo For tuning, the value of Y_MAX must be specified. E.g. 'make tune Y_MAX=100'
	@false
endif
ifndef Z_MAX
	@echo For tuning, the value of Z_MAX must be specified. E.g. 'make tune Z_MAX=100'
	@false
endif

	java -jar /home/hpc/ihpc/ihpc07/PROVA/software/easybuild/software/PATUS/0.1.4-Java-1.7.0_79/bin/patus.jar autotune bench $(X_MAX) $(Y_MAX) $(Z_MAX) @WAVE_CB_X=\(\(\$$0-4\)/4\),\(\(\$$0-4\)/2\),\(\$$0-4\) @WAVE_CB_Y=4:4:\(\$$1-4\) @WAVE_CB_Z=4:4:\(\$$2-4\) @WAVE_CHUNK=1,2,4 @WAVE__UNROLL_P3=0:1:8  | tee `date +tune_%F_%H.%M.%S.txt`
	
plot:
	for f in *.data; do echo "set term png; set view map; splot \"$$f\" matrix with image" | $(GNUPLOT) > $${f%.data}.png; done
	
clean:
	rm -rf *.o bench
