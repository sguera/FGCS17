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
  
  long width = WIDTH;
  long height = HEIGHT;

  float ***u = (float ***) malloc(sizeof(float**)*2);
  u[0] = (float **) malloc(sizeof(float*)*height);
  u[1] = (float **) malloc(sizeof(float*)*height);
  for(i=0; i<height; i++){
    u[0][i] = (float*) malloc(sizeof(float*)*width);
    u[1][i] = (float*) malloc(sizeof(float*)*width);
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
      u[0][j][i] = 1;
    }
  }
  
  double fTimeStart = seconds();
  
#pragma scop
  for (t = 0; t < T_MAX; t++)
  {
    for (j = 2; j < height - 2; j++)
    {
      for (i = 2; i < width - 2; i++)
      {
        u[(t+1)%2][j][i] = f * (
                                s0 * u[t%2][j][i] +
                                s1 * (u[t%2][j][i-1] + u[t%2][j][i+1] + u[t%2][j-1][i] + u[t%2][j+1][i]) +
                                s2 * (u[t%2][j-1][i-1] + u[t%2][j-1][i+1] + u[t%2][j+1][i-1] + u[t%2][j+1][i+1]) +
                                s4 * (u[t%2][j][i-2] + u[t%2][j][i+2] + u[t%2][j-2][i] + u[t%2][j+2][i]) +
                                s5 * (u[t%2][j-1][i-2] + u[t%2][j-2][i-1] + u[t%2][j-2][i+1] + u[t%2][j-1][i+2] +
                                      u[t%2][j+1][i-2] + u[t%2][j+2][i-1] + u[t%2][j+2][i+1] + u[t%2][j+1][i+2]
                                      )+
                                s8 * (u[t%2][j-2][i-2] + u[t%2][j-2][i+2] + u[t%2][j+2][i-2] + u[t%2][j+2][i+2])
          );
      }
    }
  }
#pragma endscop
  
  double fTimeEnd = seconds ();
  
  /* print statistics */
  double fNumFlops = (double) (width-4) * (double) (height-4) * T_MAX * 31.0;
  printf ("FLOPs in stencil code:      %e\n", fNumFlops);
  printf ("Time spent in stencil code: %f\n", fTimeEnd - fTimeStart);
  printf ("Performance in GFlop/s:     %f\n", fNumFlops / (1e9 * (fTimeEnd - fTimeStart)));
  
  return EXIT_SUCCESS;
}

// icc -O3 -fp-model precise heat_1d_np.c -o op-heat-1d-np -lm
// /* @ begin PrimeTile (num_tiling_levels=1; first_depth=1; last_depth=-1; boundary_tiling_level=-1;) @*/
// /* @ begin PrimeRegTile (scalar_replacement=0; T1t5=4; T1t6=4; T1t7=4; T1t8=4; ) @*/
// /* @ end @*/
