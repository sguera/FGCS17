/*
 * Discretized 3D heat equation stencil with non periodic boundary conditions
 * Adapted from Pochoir test bench
 *
 * Irshad Pananilath: irshad@csa.iisc.ernet.in
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
  
  long x_max = X_MAX;
  long y_max = Y_MAX;
  long z_max = Z_MAX;
  
  float ****u = (float ****) malloc(sizeof(float***)*3);
  u[0] = (float ***) malloc(sizeof(float**)*z_max);
  u[1] = (float ***) malloc(sizeof(float**)*z_max);
  u[2] = (float ***) malloc(sizeof(float**)*z_max);
  for(i=0; i<z_max; i++){
    u[0][i] = (float**) malloc(sizeof(float*)*y_max);
    u[1][i] = (float**) malloc(sizeof(float*)*y_max);
    u[2][i] = (float**) malloc(sizeof(float*)*y_max);
    for(j=0;j<y_max;j++){
      u[0][i][j] = (float*) malloc(sizeof(float)*x_max);
      u[1][i][j] = (float*) malloc(sizeof(float)*x_max);
      u[2][i][j] = (float*) malloc(sizeof(float)*x_max);
    }
  }
  
  /* Initialization */
  const float MIN = -1.f;
  const float MAX = 1.f;
  const float DX = (MAX - MIN) / (x_max - 3);
  const float DT = DX / 2.0f;
  const float DT_DX_SQUARE = DT * DT / (DX * DX);
  const int T_MAX = 100;
  /* Define our arrays */
  
  //memset (u, 0, 3 * X_MAX * Y_MAX * Z_MAX * sizeof (float));

  #pragma omp parallel for private(k,j,i)// schedule(static,1)
  for (k = 2; k < Z_MAX - 2; k++)
  {
    for (j = 2; j < Y_MAX - 2; j++)
    {
        for (i = 2; i < X_MAX - 2; i++)
        {
            float x = (i - 1) * DX + MIN;
            float y = (j - 1) * DX + MIN;
            float z = (k - 1) * DX + MIN;

            u[0][k][j][i] = u[1][k][j][i] = u[2][k][j][i] = sinf (2 * M_PI * x) * sinf (2 * M_PI * y) * sinf (2 * M_PI * z);
        }
    }
  }
  
  const float c0 = 2.0f - DT_DX_SQUARE * 7.5f;
  const float c1 = DT_DX_SQUARE * 4.0f / 3.0f;
  const float c2 = DT_DX_SQUARE * (-1.0f / 12.0f);
  
  double fTimeStart = seconds();
  
#pragma scop
  for (t = 0; t < T_MAX; t++)
  {
    for (i = 2; i < z_max - 2; i++)
    {
      for (j = 2; j < y_max - 2; j++)
      {
        for (k = 2; k < x_max - 2; k++)
        {
            u[(t+2)%3][i][j][k] = c0 * u[(t+1)%3][i][j][k] - u[t%3][i][j][k] +
            c1 * (u[(t+1)%3][i+1][j][k] + u[(t+1)%3][i-1][j][k] + u[(t+1)%3][i][j+1][k] + u[(t+1)%3][i][j-1][k] + u[(t+1)%3][i][j][k-1] + u[(t+1)%3][i][j][k+1]) +
            c2 * (u[(t+1)%3][i+2][j][k] + u[(t+1)%3][i-2][j][k] + u[(t+1)%3][i][j+2][k] + u[(t+1)%3][i][j-2][k] + u[(t+1)%3][i][j][k-2] + u[(t+1)%3][i][j][k+2]);
        }
      }
    }
  }
#pragma endscop
  
  double fTimeEnd = seconds ();
  
  /* print statistics */
  double fNumFlops = (double) (x_max-4) * (double) (y_max-4) * (double) (z_max-4) * T_MAX * 16.0;
  printf ("FLOPs in stencil code:      %e\n", fNumFlops);
  printf ("Time spent in stencil code: %f\n", fTimeEnd - fTimeStart);
  printf ("Performance in GFlop/s:     %f\n", fNumFlops / (1e9 * (fTimeEnd - fTimeStart)));
  
  return EXIT_SUCCESS;
}

// icc -O3 -fp-model precise heat_1d_np.c -o op-heat-1d-np -lm
// /* @ begin PrimeTile (num_tiling_levels=1; first_depth=1; last_depth=-1; boundary_tiling_level=-1;) @*/
// /* @ begin PrimeRegTile (scalar_replacement=0; T1t5=4; T1t6=4; T1t7=4; T1t8=4; ) @*/
// /* @ end @*/
