#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>


#include <omp.h>
#include <xmmintrin.h>
#include <emmintrin.h>
#include <ssemathlib.h>
#include <math.h>
#include "patusrt.h"

// forward_decls -->
void initialize_wave_parm(float *  U_0_m1, float *  U_0_0, float *  U_0_1, float fMin, float fDX, float fDT_DX_sq, long X_MAX, long Y_MAX, long Z_MAX, long cb_x, long cb_y, long cb_z, long chunk);
void wave_parm(float *  *  U_0_1_out, float *  U_0_m1, float *  U_0_0, float *  U_0_1, float fMin, float fDX, float fDT_DX_sq, long X_MAX, long Y_MAX, long Z_MAX, long cb_x, long cb_y, long cb_z, long chunk, int _unroll_p3);

// <--


int main (int argc, char** argv)
{
	int i;
	
	// prepare grids
	// declare_grids -->
	float *  U_0_1_out;
	float *  U_0_1_out_ref;
	float *  U_0_m1;
	float *  U_0_m1_ref;
	float *  U_0_0;
	float *  U_0_0_ref;
	float *  U_0_1;
	float *  U_0_1_ref;
	int nParamArgsCount = 0;
	char *  *  rgArgs = ((char *  * )malloc((argc*sizeof (char * ))));
	int ii;
	for (ii=1; ii<argc;  ++ ii)
	{
		if ((( * argv[ii])!='-'))
		{
			rgArgs[nParamArgsCount]=argv[ii];
			 ++ nParamArgsCount;
		}
	}
	if ((nParamArgsCount!=8))
	{
		printf("Wrong number of parameters. Syntax:\n%s <X_MAX> <Y_MAX> <Z_MAX> <cb_x> <cb_y> <cb_z> <chunk> <_unroll_p3>\n", argv[0]);
		exit(-1);
	}
	long X_MAX = atoi(rgArgs[0]);
	long Y_MAX = atoi(rgArgs[1]);
	long Z_MAX = atoi(rgArgs[2]);
	long cb_x = atoi(rgArgs[3]);
	long cb_y = atoi(rgArgs[4]);
	long cb_z = atoi(rgArgs[5]);
	long chunk = atoi(rgArgs[6]);
	long _unroll_p3 = atoi(rgArgs[7]);
	free(rgArgs);
	// <--
	
	// allocate_grids -->
	if (((X_MAX%4)!=0))
	{
		printf("Non-native, aligned SIMD type mode requires that X_MAX is divisible by 4 [X_MAX = %ld].\n", X_MAX);
		return -1;
	}
	U_0_m1=((float * )malloc(((((X_MAX*Y_MAX)*Z_MAX)*sizeof (float))+15)));
	U_0_m1_ref=((float * )malloc(((((X_MAX*Y_MAX)*Z_MAX)*sizeof (float))+15)));
	if (((X_MAX%4)!=0))
	{
		printf("Non-native, aligned SIMD type mode requires that X_MAX is divisible by 4 [X_MAX = %ld].\n", X_MAX);
		return -1;
	}
	U_0_0=((float * )malloc(((((X_MAX*Y_MAX)*Z_MAX)*sizeof (float))+15)));
	U_0_0_ref=((float * )malloc(((((X_MAX*Y_MAX)*Z_MAX)*sizeof (float))+15)));
	if (((X_MAX%4)!=0))
	{
		printf("Non-native, aligned SIMD type mode requires that X_MAX is divisible by 4 [X_MAX = %ld].\n", X_MAX);
		return -1;
	}
	U_0_1=((float * )malloc(((((X_MAX*Y_MAX)*Z_MAX)*sizeof (float))+15)));
	U_0_1_ref=((float * )malloc(((((X_MAX*Y_MAX)*Z_MAX)*sizeof (float))+15)));
	// <--
	
	
	// initialize
#pragma omp parallel
	{
		// initialize_grids -->
		initialize_wave_parm(((float * )((((uintptr_t)U_0_m1)+15)&( ~ ((uintptr_t)15)))), ((float * )((((uintptr_t)U_0_0)+15)&( ~ ((uintptr_t)15)))), ((float * )((((uintptr_t)U_0_1)+15)&( ~ ((uintptr_t)15)))), -1.0, (2.0/(X_MAX-3.0)), 0.25, X_MAX, Y_MAX, Z_MAX, cb_x, cb_y, cb_z, chunk);
		initialize_wave_parm(U_0_m1_ref, U_0_0_ref, U_0_1_ref, -1.0, (2.0/(X_MAX-3.0)), 0.25, X_MAX, Y_MAX, Z_MAX, cb_x, cb_y, cb_z, chunk);
		// <--
		
	}
	
	// write output
	if (has_arg ("-o", argc, argv))
	{
		// write_grids -->
		write_data_f("U_0_m1.0.data", 3, ((float * )((((uintptr_t)U_0_m1)+15)&( ~ ((uintptr_t)15)))), X_MAX, Y_MAX, Z_MAX);
		write_data_f("U_0_0.0.data", 3, ((float * )((((uintptr_t)U_0_0)+15)&( ~ ((uintptr_t)15)))), X_MAX, Y_MAX, Z_MAX);
		write_data_f("U_0_1.0.data", 3, ((float * )((((uintptr_t)U_0_1)+15)&( ~ ((uintptr_t)15)))), X_MAX, Y_MAX, Z_MAX);
		// <--
		
	}
	
	long nFlopsPerStencil = 20;
	long nGridPointsCount = 5 * (((100*(X_MAX-4))*(Y_MAX-4))*(Z_MAX-4));
	long nBytesTransferred = 5 * (100*(((((X_MAX*Y_MAX)*Z_MAX)*sizeof (float))+(((X_MAX*Y_MAX)*Z_MAX)*sizeof (float)))+(((X_MAX*Y_MAX)*Z_MAX)*sizeof (float))));
	
	// warm up
#pragma omp parallel
	{
		// compute_stencil -->
		wave_parm(( & U_0_1_out), ((float * )((((uintptr_t)U_0_m1)+15)&( ~ ((uintptr_t)15)))), ((float * )((((uintptr_t)U_0_0)+15)&( ~ ((uintptr_t)15)))), ((float * )((((uintptr_t)U_0_1)+15)&( ~ ((uintptr_t)15)))), -1.0, (2.0/(X_MAX-3.0)), 0.25, X_MAX, Y_MAX, Z_MAX, cb_x, cb_y, cb_z, chunk, _unroll_p3);
		// <--
		
	}
	
	// run the benchmark
	tic ();
#pragma omp parallel private(i)
	for (i = 0; i < 5; i++)
	{
		// compute_stencil -->
		wave_parm(( & U_0_1_out), ((float * )((((uintptr_t)U_0_m1)+15)&( ~ ((uintptr_t)15)))), ((float * )((((uintptr_t)U_0_0)+15)&( ~ ((uintptr_t)15)))), ((float * )((((uintptr_t)U_0_1)+15)&( ~ ((uintptr_t)15)))), -1.0, (2.0/(X_MAX-3.0)), 0.25, X_MAX, Y_MAX, Z_MAX, cb_x, cb_y, cb_z, chunk, _unroll_p3);
		// <--
		
#pragma omp barrier
	}
	toc (nFlopsPerStencil, nGridPointsCount, nBytesTransferred);
	
	// write output
	if (has_arg ("-o", argc, argv))
	{
#pragma omp parallel
		{
			// initialize_grids -->
			initialize_wave_parm(((float * )((((uintptr_t)U_0_m1)+15)&( ~ ((uintptr_t)15)))), ((float * )((((uintptr_t)U_0_0)+15)&( ~ ((uintptr_t)15)))), ((float * )((((uintptr_t)U_0_1)+15)&( ~ ((uintptr_t)15)))), -1.0, (2.0/(X_MAX-3.0)), 0.25, X_MAX, Y_MAX, Z_MAX, cb_x, cb_y, cb_z, chunk);
			initialize_wave_parm(U_0_m1_ref, U_0_0_ref, U_0_1_ref, -1.0, (2.0/(X_MAX-3.0)), 0.25, X_MAX, Y_MAX, Z_MAX, cb_x, cb_y, cb_z, chunk);
			// <--
			
#pragma omp barrier
			// compute_stencil -->
			wave_parm(( & U_0_1_out), ((float * )((((uintptr_t)U_0_m1)+15)&( ~ ((uintptr_t)15)))), ((float * )((((uintptr_t)U_0_0)+15)&( ~ ((uintptr_t)15)))), ((float * )((((uintptr_t)U_0_1)+15)&( ~ ((uintptr_t)15)))), -1.0, (2.0/(X_MAX-3.0)), 0.25, X_MAX, Y_MAX, Z_MAX, cb_x, cb_y, cb_z, chunk, _unroll_p3);
			// <--
			
		}
		// write_grids -->
		write_data_f("U_0_1_out.1.data", 3, U_0_1_out, X_MAX, Y_MAX, Z_MAX);
		// <--
		
	}
	
	// validate
	if (1)
	{
#pragma omp parallel
		{
			// initialize_grids -->
			initialize_wave_parm(((float * )((((uintptr_t)U_0_m1)+15)&( ~ ((uintptr_t)15)))), ((float * )((((uintptr_t)U_0_0)+15)&( ~ ((uintptr_t)15)))), ((float * )((((uintptr_t)U_0_1)+15)&( ~ ((uintptr_t)15)))), -1.0, (2.0/(X_MAX-3.0)), 0.25, X_MAX, Y_MAX, Z_MAX, cb_x, cb_y, cb_z, chunk);
			initialize_wave_parm(U_0_m1_ref, U_0_0_ref, U_0_1_ref, -1.0, (2.0/(X_MAX-3.0)), 0.25, X_MAX, Y_MAX, Z_MAX, cb_x, cb_y, cb_z, chunk);
			// <--
			
#pragma omp barrier
			// compute_stencil -->
			wave_parm(( & U_0_1_out), ((float * )((((uintptr_t)U_0_m1)+15)&( ~ ((uintptr_t)15)))), ((float * )((((uintptr_t)U_0_0)+15)&( ~ ((uintptr_t)15)))), ((float * )((((uintptr_t)U_0_1)+15)&( ~ ((uintptr_t)15)))), -1.0, (2.0/(X_MAX-3.0)), 0.25, X_MAX, Y_MAX, Z_MAX, cb_x, cb_y, cb_z, chunk, _unroll_p3);
			// <--
			
		}
		// validate_computation -->
		int bHasErrors = 0;
		int _idx0;
		int _idx10;
		int _idx11;
		int _idx12;
		int _idx13;
		int _idx2;
		int _idx3;
		int _idx4;
		int _idx5;
		int _idx6;
		int _idx7;
		int _idx8;
		int _idx9;
		float c1;
		float c2;
		float c3;
		int p3_loc_idx_x;
		int p3_loc_idx_y;
		int p3_loc_idx_z;
		int t_ref;
		float *  tmp_swap_2;
		{
			/*
			for t_ref = 1..100 by 1 parallel 1 <level 0> schedule 1 { ... }
			*/
			for (t_ref=1; t_ref<=100; t_ref+=1)
			{
				/*
				for POINT pt_ref[t=t][0] of size [1, 1, 1] in u0[t=t][0] + [ min=[0, 0, 0], max=[0, 0, 0] ] parallel 1 <level 0> schedule default { ... }
				*/
				{
					/* Index bounds calculations for iterators in pt_ref[t=t][0] */
					for (p3_loc_idx_z=2; p3_loc_idx_z<((Z_MAX-3)+1); p3_loc_idx_z+=1)
					{
						for (p3_loc_idx_y=2; p3_loc_idx_y<((Y_MAX-3)+1); p3_loc_idx_y+=1)
						{
							for (p3_loc_idx_x=2; p3_loc_idx_x<((X_MAX-3)+1); p3_loc_idx_x+=1)
							{
								/*
								pt_ref[t=t][0]=stencil(pt_ref[t=t][0])
								*/
								/* _idx0 = ((X_MAX*((Y_MAX*p3_loc_idx_z)+p3_loc_idx_y))+p3_loc_idx_x) */
								_idx0=((X_MAX*((Y_MAX*p3_loc_idx_z)+p3_loc_idx_y))+p3_loc_idx_x);
								/* _idx2 = ((X_MAX*((Y_MAX*p3_loc_idx_z)+p3_loc_idx_y))+(p3_loc_idx_x+1)) */
								_idx2=(_idx0+1);
								/* _idx3 = ((X_MAX*((Y_MAX*p3_loc_idx_z)+p3_loc_idx_y))+(p3_loc_idx_x-1)) */
								_idx3=(_idx0-1);
								/* _idx4 = ((X_MAX*((Y_MAX*p3_loc_idx_z)+(p3_loc_idx_y+1)))+p3_loc_idx_x) */
								_idx4=(X_MAX+_idx0);
								/* _idx5 = ((X_MAX*((Y_MAX*p3_loc_idx_z)+(p3_loc_idx_y-1)))+p3_loc_idx_x) */
								_idx5=(_idx0-X_MAX);
								/* _idx6 = ((X_MAX*((Y_MAX*(p3_loc_idx_z+1))+p3_loc_idx_y))+p3_loc_idx_x) */
								_idx6=((X_MAX*Y_MAX)+_idx0);
								/* _idx7 = ((X_MAX*((Y_MAX*(p3_loc_idx_z-1))+p3_loc_idx_y))+p3_loc_idx_x) */
								_idx7=(_idx0-(X_MAX*Y_MAX));
								/* _idx8 = ((X_MAX*((Y_MAX*p3_loc_idx_z)+p3_loc_idx_y))+(p3_loc_idx_x+2)) */
								_idx8=(_idx0+2);
								/* _idx9 = ((X_MAX*((Y_MAX*p3_loc_idx_z)+p3_loc_idx_y))+(p3_loc_idx_x-2)) */
								_idx9=(_idx0-2);
								/* _idx10 = ((X_MAX*((Y_MAX*p3_loc_idx_z)+(p3_loc_idx_y+2)))+p3_loc_idx_x) */
								_idx10=((2*X_MAX)+_idx0);
								/* _idx11 = ((X_MAX*((Y_MAX*p3_loc_idx_z)+(p3_loc_idx_y-2)))+p3_loc_idx_x) */
								_idx11=(_idx5-X_MAX);
								/* _idx12 = ((X_MAX*((Y_MAX*(p3_loc_idx_z+2))+p3_loc_idx_y))+p3_loc_idx_x) */
								_idx12=((X_MAX*Y_MAX)+_idx6);
								/* _idx13 = ((X_MAX*((Y_MAX*(p3_loc_idx_z-2))+p3_loc_idx_y))+p3_loc_idx_x) */
								_idx13=(_idx7-(X_MAX*Y_MAX));
								c1=(2.0f-(0.25*7.5f));
								c2=(0.25*1.3333333333333333f);
								c3=(0.25*-0.08333333333333333f);
								U_0_1_ref[_idx0]=(((c1*U_0_0_ref[_idx0])-U_0_m1_ref[_idx0])+((c2*((U_0_0_ref[_idx2]+(U_0_0_ref[_idx3]+U_0_0_ref[_idx4]))+(U_0_0_ref[_idx5]+(U_0_0_ref[_idx6]+U_0_0_ref[_idx7]))))+(c3*((U_0_0_ref[_idx8]+(U_0_0_ref[_idx9]+U_0_0_ref[_idx10]))+(U_0_0_ref[_idx11]+(U_0_0_ref[_idx12]+U_0_0_ref[_idx13]))))));
							}
						}
					}
				}
				U_0_1_out_ref=U_0_1_ref;
				tmp_swap_2=U_0_m1_ref;
				U_0_m1_ref=U_0_0_ref;
				U_0_0_ref=U_0_1_ref;
				U_0_1_ref=tmp_swap_2;
#pragma omp barrier 
			}
		}
		/*
		for POINT pt_ref[t=t][0] of size [1, 1, 1] in u0[t=t][0] + [ min=[0, 0, 0], max=[0, 0, 0] ] parallel 1 <level 0> schedule default { ... }
		*/
		{
			/* Index bounds calculations for iterators in pt_ref[t=t][0] */
			for (p3_loc_idx_z=2; p3_loc_idx_z<((Z_MAX-3)+1); p3_loc_idx_z+=1)
			{
				for (p3_loc_idx_y=2; p3_loc_idx_y<((Y_MAX-3)+1); p3_loc_idx_y+=1)
				{
					for (p3_loc_idx_x=2; p3_loc_idx_x<((X_MAX-3)+1); p3_loc_idx_x+=1)
					{
						/* _idx0 = ((X_MAX*((Y_MAX*p3_loc_idx_z)+p3_loc_idx_y))+p3_loc_idx_x) */
						_idx0=((X_MAX*((Y_MAX*p3_loc_idx_z)+p3_loc_idx_y))+p3_loc_idx_x);
						if ((fabs(((U_0_1_out[_idx0]-U_0_1_out_ref[_idx0])/U_0_1_out_ref[_idx0]))>1.0E-4))
						{
							bHasErrors=1;
							break;
						}
					}
				}
			}
		}
		// <--
		
		if (( ! bHasErrors))
		puts ("Validation OK.");
		else
		{
			// deallocate_grids -->
			free(U_0_m1);
			free(U_0_m1_ref);
			free(U_0_0);
			free(U_0_0_ref);
			free(U_0_1);
			free(U_0_1_ref);
			// <--
			
			puts ("Validation failed.");
			return -1;
		}
	}
	
	// free memory
	// deallocate_grids -->
	free(U_0_m1);
	free(U_0_m1_ref);
	free(U_0_0);
	free(U_0_0_ref);
	free(U_0_1);
	free(U_0_1_ref);
	// <--
	
	
	return 0;
}
