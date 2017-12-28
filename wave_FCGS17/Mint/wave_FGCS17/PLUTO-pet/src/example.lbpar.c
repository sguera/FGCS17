#include <omp.h>
#include <math.h>
#define ceild(n,d)  ceil(((double)(n))/((double)(d)))
#define floord(n,d) floor(((double)(n))/((double)(d)))
#define max(x,y)    ((x) > (y)? (x) : (y))
#define min(x,y)    ((x) < (y)? (x) : (y))

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
  
/* Copyright (C) 1991-2014 Free Software Foundation, Inc.
   This file is part of the GNU C Library.

   The GNU C Library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   The GNU C Library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with the GNU C Library; if not, see
   <http://www.gnu.org/licenses/>.  */
/* This header is separate from features.h so that the compiler can
   include it implicitly at the start of every compilation.  It must
   not itself include <features.h> or any other header that includes
   <features.h> because the implicit include comes before any feature
   test macros that may be defined in a source file before it first
   explicitly includes a system header.  GCC knows the name of this
   header in order to preinclude it.  */
/* glibc's intent is to support the IEC 559 math functionality, real
   and complex.  If the GCC (4.9 and later) predefined macros
   specifying compiler intent are available, use them to determine
   whether the overall intent is to support these features; otherwise,
   presume an older compiler has intent to support these features and
   define these macros by default.  */
/* wchar_t uses ISO/IEC 10646 (2nd ed., published 2011-03-15) /
   Unicode 6.0.  */
/* We do not support C11 <threads.h>.  */
  int t1, t2, t3, t4, t5, t6, t7, t8;
 int lb, ub, lbp, ubp, lb2, ub2;
 register int lbv, ubv;
/* Start of CLooG code */
if ((T_MAX >= 1) && (T_MAX <= 2147483646) && (x_max >= 5) && (x_max <= 2147483648) && (y_max >= 5) && (y_max <= 2147483648) && (z_max >= 5) && (z_max <= 2147483648)) {
  for (t1=-1;t1<=floord(T_MAX-1,4);t1++) {
    lbp=max(ceild(t1,2),ceild(8*t1-T_MAX+2,8));
    ubp=min(floord(2*T_MAX+z_max-5,16),floord(8*t1+z_max+4,16));
#pragma omp parallel for private(lbv,ubv,t3,t4,t5,t6,t7,t8)
    for (t2=lbp;t2<=ubp;t2++) {
      for (t3=max(max(0,ceild(t1-1,2)),ceild(16*t2-z_max-9,16));t3<=min(min(min(floord(2*T_MAX+y_max-5,16),floord(8*t1+y_max+11,16)),floord(16*t2+y_max+9,16)),floord(16*t1-16*t2+z_max+y_max+9,16));t3++) {
        for (t4=max(max(max(0,ceild(t1-124,125)),ceild(16*t2-z_max-993,1000)),ceild(16*t3-y_max-993,1000));t4<=min(min(min(min(floord(2*T_MAX+x_max-5,1000),floord(8*t1+x_max+11,1000)),floord(16*t2+x_max+9,1000)),floord(16*t3+x_max+9,1000)),floord(16*t1-16*t2+z_max+x_max+9,1000));t4++) {
          for (t5=max(max(max(max(max(0,ceild(16*t2-z_max+3,2)),ceild(16*t3-y_max+3,2)),ceild(1000*t4-x_max+3,2)),4*t1),8*t1-8*t2+1);t5<=min(min(min(min(min(floord(16*t1-16*t2+z_max+12,2),T_MAX-1),4*t1+7),8*t2+6),8*t3+6),500*t4+498);t5++) {
            for (t6=max(max(16*t2,2*t5+2),-16*t1+16*t2+4*t5-15);t6<=min(min(16*t2+15,-16*t1+16*t2+4*t5),2*t5+z_max-3);t6++) {
              for (t7=max(16*t3,2*t5+2);t7<=min(16*t3+15,2*t5+y_max-3);t7++) {
                lbv=max(1000*t4,2*t5+2);
                ubv=min(1000*t4+999,2*t5+x_max-3);
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  u[( t5 + 2) % 3][ (-2*t5+t6)][ (-2*t5+t7)][ (-2*t5+t8)] = ((((c0 * u[( t5 + 1) % 3][ (-2*t5+t6)][ (-2*t5+t7)][ (-2*t5+t8)]) - u[ t5 % 3][ (-2*t5+t6)][ (-2*t5+t7)][ (-2*t5+t8)]) + (c1 * (((((u[( t5 + 1) % 3][ (-2*t5+t6) + 1][ (-2*t5+t7)][ (-2*t5+t8)] + u[( t5 + 1) % 3][ (-2*t5+t6) - 1][ (-2*t5+t7)][ (-2*t5+t8)]) + u[( t5 + 1) % 3][ (-2*t5+t6)][ (-2*t5+t7) + 1][ (-2*t5+t8)]) + u[( t5 + 1) % 3][ (-2*t5+t6)][ (-2*t5+t7) - 1][ (-2*t5+t8)]) + u[( t5 + 1) % 3][ (-2*t5+t6)][ (-2*t5+t7)][ (-2*t5+t8) - 1]) + u[( t5 + 1) % 3][ (-2*t5+t6)][ (-2*t5+t7)][ (-2*t5+t8) + 1]))) + (c2 * (((((u[( t5 + 1) % 3][ (-2*t5+t6) + 2][ (-2*t5+t7)][ (-2*t5+t8)] + u[( t5 + 1) % 3][ (-2*t5+t6) - 2][ (-2*t5+t7)][ (-2*t5+t8)]) + u[( t5 + 1) % 3][ (-2*t5+t6)][ (-2*t5+t7) + 2][ (-2*t5+t8)]) + u[( t5 + 1) % 3][ (-2*t5+t6)][ (-2*t5+t7) - 2][ (-2*t5+t8)]) + u[( t5 + 1) % 3][ (-2*t5+t6)][ (-2*t5+t7)][ (-2*t5+t8) - 2]) + u[( t5 + 1) % 3][ (-2*t5+t6)][ (-2*t5+t7)][ (-2*t5+t8) + 2])));;
                }
              }
            }
          }
        }
      }
    }
  }
}
/* End of CLooG code */
  
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
