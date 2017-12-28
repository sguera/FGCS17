/*
 * Gaussian filter
 *
 * Danilo Guerrera: danilo.guerrera@unibas.ch
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <math.h>
#include "dim_input.h"

#define IDX(j,i) ((i)+height*(j))

/**
 * Get current time in seconds.
 */
double seconds ()
{
  struct timeval tv;
  gettimeofday (&tv, NULL);
  return ((double) tv.tv_sec) + 1e-6 * tv.tv_usec;
}

int malloc_error (const char* err)
{
  fprintf (stderr, "Failed to allocate the field %s.\n", err);
  return EXIT_FAILURE;
}


int main(int argc, char * argv[])
{
  int t, i, j, k;
  double time = 0.0;
  double nFlops, GFlops;
  long width = WIDTH;
  long height = HEIGHT;
  
  float* u_0 = (float *) malloc(width * height * sizeof(float));
  if (u_0 == NULL)
  {
      return malloc_error ("u_0");
  }
  float* u_1 = (float *) malloc(width * height * sizeof(float));
  if (u_1 == NULL)
  {
      return malloc_error ("u_1");
  }
  
  /* Initialization */
  const float SIGMA = 3.0;
  const float f0 = 1 / (2 * SIGMA * SIGMA);
  const float s0 = exp ( 0 * f0);
  const float s1 = exp (-1 * f0);
  const float s2 = exp (-2 * f0);
  const float s4 = exp (-4 * f0);
  const float s5 = exp (-5 * f0);
  const float s8 = exp (-8 * f0);
  const float f = 1 / (s0 + 4 * (s1 + s2 + s4 + s8) + 8 * s5);
  const int T_MAX = 50;
  /* Define our arrays */
  
  //memset (u, 0, 3 * WIDTH * height * Z_MAX * sizeof (float));

  #pragma omp parallel for private(j,i)// schedule(static,1)
  for (j = 2; j < height - 2; j++)
  {
    for (i = 2; i < width - 2; i++)
    {
      u_0[IDX(j,i)] = 1;
    }
  }
  

#ifdef LIKWID_PERFMON
  #pragma omp parallel
  {
    LIKWID_MARKER_START("Compute");
  }
#endif

  double fTimeStart = seconds();
  
  #pragma omp parallel for private(j,i)
  for (t = 0; t < T_MAX; t++) {
    for (j = 2; j < height - 2; j++)
    {
      for (i = 2; i < width - 2; i++)
      {
        u_1[IDX(j,i)] = f * (
                          s0 * u_0[IDX(j,i)] +
                          s1 * (u_0[IDX(j,i-1)] + u_0[IDX(j,i+1)] + u_0[IDX(j-1,i)] + u_0[IDX(j+1,i)]) +
                          s2 * (u_0[IDX(j-1,i-1)] + u_0[IDX(j-1,i+1)] + u_0[IDX(j+1,i-1)] + u_0[IDX(j+1,i+1)]) +
                          s4 * (u_0[IDX(j,i-2)] + u_0[IDX(j,i+2)] + u_0[IDX(j-2,i)] + u_0[IDX(j+2,i)]) +
                          s5 * (
                                u_0[IDX(j-1,i-2)] + u_0[IDX(j-2,i-1)] + u_0[IDX(j-2,i+1)] + u_0[IDX(j-1,i+2)] +
                                u_0[IDX(j+1,i-2)] + u_0[IDX(j+2,i-1)] + u_0[IDX(j+2,i+1)] + u_0[IDX(j+1,i+2)]
                                )+
                          s8 * (u_0[IDX(j-2,i-2)] + u_0[IDX(j-2,i+2)] + u_0[IDX(j+2,i-2)] + u_0[IDX(j+2,i+2)])
                        );
      }
    }
  }
  
  double fTimeEnd = seconds ();

#ifdef LIKWID_PERFMON
  #pragma omp parallel
  {
    LIKWID_MARKER_STOP("Compute");
  }
  LIKWID_MARKER_CLOSE;
#endif
  
  /* print statistics */
  double fNumFlops = (double) (width-4) * (double) (height-4) * T_MAX * 31.0;
  printf ("FLOPs in stencil code:      %e\n", fNumFlops);
  printf ("Time spent in stencil code: %f\n", fTimeEnd - fTimeStart);
  printf ("Performance in GFlop/s:     %f\n", fNumFlops / (1e9 * (fTimeEnd - fTimeStart)));

  return EXIT_SUCCESS;
}
