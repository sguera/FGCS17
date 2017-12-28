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
ifndef WIDTH
	@echo For tuning, the value of WIDTH must be specified. E.g. 'make tune WIDTH=100'
	@false
endif
ifndef HEIGHT
	@echo For tuning, the value of HEIGHT must be specified. E.g. 'make tune HEIGHT=100'
	@false
endif

	java -jar /home/hpc/ihpc/ihpc07/PROVA/software/easybuild/software/PATUS/0.1.4-Java-1.7.0_79/bin/patus.jar autotune bench $(WIDTH) $(HEIGHT) @BLUR_CB_X=\(\(\$$0-4\)/4\),\(\(\$$0-4\)/2\),\(\$$0-4\) @BLUR_CB_Y=4:4:\(\$$1-4\) @BLUR_CHUNK=1,2,4 @BLUR__UNROLL_P3=0:1:4  | tee `date +tune_%F_%H.%M.%S.txt`
	
plot:
	for f in *.data; do echo "set term png; set view map; splot \"$$f\" matrix with image" | $(GNUPLOT) > $${f%.data}.png; done
	
clean:
	rm -rf *.o bench
