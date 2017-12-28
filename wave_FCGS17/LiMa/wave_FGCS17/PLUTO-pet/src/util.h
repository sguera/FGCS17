#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <ctype.h>

void get_time(struct timespec* t) {
	clock_gettime(CLOCK_MONOTONIC, t);
}
void get_clockres(struct timespec* t) {
	clock_getres(CLOCK_MONOTONIC, t);
}

struct timespec t1, t2, dt;

void start_timing() {

	get_clockres(&t1);
	printf("Timer resolution is %lu nano seconds.\n",t1.tv_nsec);

	get_time(&t1);
}

double stop_timing() {

	get_time(&t2);

	if ((t2.tv_nsec - t1.tv_nsec) < 0) {
		dt.tv_sec = t2.tv_sec - t1.tv_sec - 1;
		dt.tv_nsec = 1000000000 - t1.tv_nsec + t2.tv_nsec;
	} else {
		dt.tv_sec = t2.tv_sec - t1.tv_sec;
		dt.tv_nsec = t2.tv_nsec - t1.tv_nsec;
	}

	return (dt.tv_sec + (double)(dt.tv_nsec) * 0.000000001);

}

void printMetric(char metric[], double val) {

        printf("Performance in %s: %f\n", metric, val);
}

