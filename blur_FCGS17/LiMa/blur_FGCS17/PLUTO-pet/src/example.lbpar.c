#include <omp.h>
#include <math.h>
#define ceild(n,d)  ceil(((double)(n))/((double)(d)))
#define floord(n,d) floor(((double)(n))/((double)(d)))
#define max(x,y)    ((x) > (y)? (x) : (y))
#define min(x,y)    ((x) < (y)? (x) : (y))

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
  
/* Copyright (C) 1991-2012 Free Software Foundation, Inc.
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
/* We do support the IEC 559 math functionality, real and complex.  */
/* wchar_t uses ISO/IEC 10646 (2nd ed., published 2011-03-15) /
   Unicode 6.0.  */
/* We do not support C11 <threads.h>.  */
  int t1, t2, t3, t4, t5, t6;
 int lb, ub, lbp, ubp, lb2, ub2;
 register int lbv, ubv;
/* Start of CLooG code */
if ((T_MAX >= 1) && (height >= 5) && (height <= 2147483648) && (width >= 5) && (width <= 2147483648)) {
  for (t1=-1;t1<=floord(T_MAX-1,4);t1++) {
    lbp=max(ceild(t1,2),ceild(8*t1-T_MAX+2,8));
    ubp=min(floord(2*T_MAX+height-5,16),floord(8*t1+height+4,16));
#pragma omp parallel for private(lbv,ubv,t3,t4,t5,t6)
    for (t2=lbp;t2<=ubp;t2++) {
      for (t3=max(max(0,ceild(t1-1,2)),ceild(16*t2-height-9,16));t3<=min(min(min(floord(2*T_MAX+width-5,16),floord(8*t1+width+11,16)),floord(16*t2+width+9,16)),floord(16*t1-16*t2+height+width+9,16));t3++) {
        for (t4=max(max(max(max(0,ceild(16*t2-height+3,2)),ceild(16*t3-width+3,2)),4*t1),8*t1-8*t2+1);t4<=min(min(min(min(floord(16*t1-16*t2+height+12,2),T_MAX-1),4*t1+7),8*t2+6),8*t3+6);t4++) {
          for (t5=max(max(16*t2,2*t4+2),-16*t1+16*t2+4*t4-15);t5<=min(min(16*t2+15,-16*t1+16*t2+4*t4),2*t4+height-3);t5++) {
            lbv=max(16*t3,2*t4+2);
            ubv=min(16*t3+15,2*t4+width-3);
#pragma ivdep
#pragma vector always
            for (t6=lbv;t6<=ubv;t6++) {
              u[( t4 + 1) % 2][ (-2*t4+t5)][ (-2*t4+t6)] = (f * ((((((s0 * u[ t4 % 2][ (-2*t4+t5)][ (-2*t4+t6)]) + (s1 * (((u[ t4 % 2][ (-2*t4+t5)][ (-2*t4+t6) - 1] + u[ t4 % 2][ (-2*t4+t5)][ (-2*t4+t6) + 1]) + u[ t4 % 2][ (-2*t4+t5) - 1][ (-2*t4+t6)]) + u[ t4 % 2][ (-2*t4+t5) + 1][ (-2*t4+t6)]))) + (s2 * (((u[ t4 % 2][ (-2*t4+t5) - 1][ (-2*t4+t6) - 1] + u[ t4 % 2][ (-2*t4+t5) - 1][ (-2*t4+t6) + 1]) + u[ t4 % 2][ (-2*t4+t5) + 1][ (-2*t4+t6) - 1]) + u[ t4 % 2][ (-2*t4+t5) + 1][ (-2*t4+t6) + 1]))) + (s4 * (((u[ t4 % 2][ (-2*t4+t5)][ (-2*t4+t6) - 2] + u[ t4 % 2][ (-2*t4+t5)][ (-2*t4+t6) + 2]) + u[ t4 % 2][ (-2*t4+t5) - 2][ (-2*t4+t6)]) + u[ t4 % 2][ (-2*t4+t5) + 2][ (-2*t4+t6)]))) + (s5 * (((((((u[ t4 % 2][ (-2*t4+t5) - 1][ (-2*t4+t6) - 2] + u[ t4 % 2][ (-2*t4+t5) - 2][ (-2*t4+t6) - 1]) + u[ t4 % 2][ (-2*t4+t5) - 2][ (-2*t4+t6) + 1]) + u[ t4 % 2][ (-2*t4+t5) - 1][ (-2*t4+t6) + 2]) + u[ t4 % 2][ (-2*t4+t5) + 1][ (-2*t4+t6) - 2]) + u[ t4 % 2][ (-2*t4+t5) + 2][ (-2*t4+t6) - 1]) + u[ t4 % 2][ (-2*t4+t5) + 2][ (-2*t4+t6) + 1]) + u[ t4 % 2][ (-2*t4+t5) + 1][ (-2*t4+t6) + 2]))) + (s8 * (((u[ t4 % 2][ (-2*t4+t5) - 2][ (-2*t4+t6) - 2] + u[ t4 % 2][ (-2*t4+t5) - 2][ (-2*t4+t6) + 2]) + u[ t4 % 2][ (-2*t4+t5) + 2][ (-2*t4+t6) - 2]) + u[ t4 % 2][ (-2*t4+t5) + 2][ (-2*t4+t6) + 2]))));;
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
