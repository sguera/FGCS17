#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>


#include <omp.h>
#include <immintrin.h>
#include <x86intrin.h>
#include <avxmathlib.h>
#include <math.h>
#include "patusrt.h"

// forward_decls -->
void initialize_blur_parm(float *  U_0_0, float *  U_0_1, float sigma, float f0, float s0, float s1, float s2, float s4, float s5, float s8, float f, long WIDTH, long HEIGHT, long cb_x, long cb_y, long chunk);
void blur_parm(float *  *  U_0_1_out, float *  U_0_0, float *  U_0_1, float sigma, float f0, float s0, float s1, float s2, float s4, float s5, float s8, float f, long WIDTH, long HEIGHT, long cb_x, long cb_y, long chunk, int _unroll_p3);

// <--


int main (int argc, char** argv)
{
	int i;
	
	// prepare grids
	// declare_grids -->
	float *  U_0_1_out;
	float *  U_0_1_out_ref;
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
	if ((nParamArgsCount!=6))
	{
		printf("Wrong number of parameters. Syntax:\n%s <WIDTH> <HEIGHT> <cb_x> <cb_y> <chunk> <_unroll_p3>\n", argv[0]);
		exit(-1);
	}
	long WIDTH = atoi(rgArgs[0]);
	long HEIGHT = atoi(rgArgs[1]);
	long cb_x = atoi(rgArgs[2]);
	long cb_y = atoi(rgArgs[3]);
	long chunk = atoi(rgArgs[4]);
	long _unroll_p3 = atoi(rgArgs[5]);
	free(rgArgs);
	// <--
	
	// allocate_grids -->
	if (((WIDTH%8)!=0))
	{
		printf("Non-native, aligned SIMD type mode requires that WIDTH is divisible by 8 [WIDTH = %ld].\n", WIDTH);
		return -1;
	}
	U_0_0=((float * )malloc((((WIDTH*HEIGHT)*sizeof (float))+31)));
	U_0_0_ref=((float * )malloc((((WIDTH*HEIGHT)*sizeof (float))+31)));
	if (((WIDTH%8)!=0))
	{
		printf("Non-native, aligned SIMD type mode requires that WIDTH is divisible by 8 [WIDTH = %ld].\n", WIDTH);
		return -1;
	}
	U_0_1=((float * )malloc((((WIDTH*HEIGHT)*sizeof (float))+31)));
	U_0_1_ref=((float * )malloc((((WIDTH*HEIGHT)*sizeof (float))+31)));
	// <--
	
	
	// initialize
#pragma omp parallel
	{
		// initialize_grids -->
		initialize_blur_parm(((float * )((((uintptr_t)U_0_0)+31)&( ~ ((uintptr_t)31)))), ((float * )((((uintptr_t)U_0_1)+31)&( ~ ((uintptr_t)31)))), 3.0, 0.05555555555555555, exp(1.0), exp(-0.05555555555555555), exp(-0.1111111111111111), exp(-0.2222222222222222), exp(-0.2777777777777778), exp(-0.4444444444444444), (1.0/(exp(1.0)+((((exp(-0.05555555555555555)+exp(-0.1111111111111111))+(exp(-0.2222222222222222)+exp(-0.4444444444444444)))*4.0)+(exp(-0.2777777777777778)*8.0)))), WIDTH, HEIGHT, cb_x, cb_y, chunk);
		initialize_blur_parm(U_0_0_ref, U_0_1_ref, 3.0, 0.05555555555555555, exp(1.0), exp(-0.05555555555555555), exp(-0.1111111111111111), exp(-0.2222222222222222), exp(-0.2777777777777778), exp(-0.4444444444444444), (1.0/(exp(1.0)+((((exp(-0.05555555555555555)+exp(-0.1111111111111111))+(exp(-0.2222222222222222)+exp(-0.4444444444444444)))*4.0)+(exp(-0.2777777777777778)*8.0)))), WIDTH, HEIGHT, cb_x, cb_y, chunk);
		// <--
		
	}
	
	// write output
	if (has_arg ("-o", argc, argv))
	{
		// write_grids -->
		write_data_f("U_0_0.0.data", 2, ((float * )((((uintptr_t)U_0_0)+31)&( ~ ((uintptr_t)31)))), WIDTH, HEIGHT);
		write_data_f("U_0_1.0.data", 2, ((float * )((((uintptr_t)U_0_1)+31)&( ~ ((uintptr_t)31)))), WIDTH, HEIGHT);
		// <--
		
	}
	
	long nFlopsPerStencil = 31;
	long nGridPointsCount = 5 * ((50*(HEIGHT-4))*(WIDTH-4));
	long nBytesTransferred = 5 * (50*(((WIDTH*HEIGHT)*sizeof (float))+((WIDTH*HEIGHT)*sizeof (float))));
	
	// warm up
#pragma omp parallel
	{
		// compute_stencil -->
		blur_parm(( & U_0_1_out), ((float * )((((uintptr_t)U_0_0)+31)&( ~ ((uintptr_t)31)))), ((float * )((((uintptr_t)U_0_1)+31)&( ~ ((uintptr_t)31)))), 3.0, 0.05555555555555555, exp(1.0), exp(-0.05555555555555555), exp(-0.1111111111111111), exp(-0.2222222222222222), exp(-0.2777777777777778), exp(-0.4444444444444444), (1.0/(exp(1.0)+((((exp(-0.05555555555555555)+exp(-0.1111111111111111))+(exp(-0.2222222222222222)+exp(-0.4444444444444444)))*4.0)+(exp(-0.2777777777777778)*8.0)))), WIDTH, HEIGHT, cb_x, cb_y, chunk, _unroll_p3);
		// <--
		
	}
	
	// run the benchmark
	tic ();
#pragma omp parallel private(i)
	for (i = 0; i < 5; i++)
	{
		// compute_stencil -->
		blur_parm(( & U_0_1_out), ((float * )((((uintptr_t)U_0_0)+31)&( ~ ((uintptr_t)31)))), ((float * )((((uintptr_t)U_0_1)+31)&( ~ ((uintptr_t)31)))), 3.0, 0.05555555555555555, exp(1.0), exp(-0.05555555555555555), exp(-0.1111111111111111), exp(-0.2222222222222222), exp(-0.2777777777777778), exp(-0.4444444444444444), (1.0/(exp(1.0)+((((exp(-0.05555555555555555)+exp(-0.1111111111111111))+(exp(-0.2222222222222222)+exp(-0.4444444444444444)))*4.0)+(exp(-0.2777777777777778)*8.0)))), WIDTH, HEIGHT, cb_x, cb_y, chunk, _unroll_p3);
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
			initialize_blur_parm(((float * )((((uintptr_t)U_0_0)+31)&( ~ ((uintptr_t)31)))), ((float * )((((uintptr_t)U_0_1)+31)&( ~ ((uintptr_t)31)))), 3.0, 0.05555555555555555, exp(1.0), exp(-0.05555555555555555), exp(-0.1111111111111111), exp(-0.2222222222222222), exp(-0.2777777777777778), exp(-0.4444444444444444), (1.0/(exp(1.0)+((((exp(-0.05555555555555555)+exp(-0.1111111111111111))+(exp(-0.2222222222222222)+exp(-0.4444444444444444)))*4.0)+(exp(-0.2777777777777778)*8.0)))), WIDTH, HEIGHT, cb_x, cb_y, chunk);
			initialize_blur_parm(U_0_0_ref, U_0_1_ref, 3.0, 0.05555555555555555, exp(1.0), exp(-0.05555555555555555), exp(-0.1111111111111111), exp(-0.2222222222222222), exp(-0.2777777777777778), exp(-0.4444444444444444), (1.0/(exp(1.0)+((((exp(-0.05555555555555555)+exp(-0.1111111111111111))+(exp(-0.2222222222222222)+exp(-0.4444444444444444)))*4.0)+(exp(-0.2777777777777778)*8.0)))), WIDTH, HEIGHT, cb_x, cb_y, chunk);
			// <--
			
#pragma omp barrier
			// compute_stencil -->
			blur_parm(( & U_0_1_out), ((float * )((((uintptr_t)U_0_0)+31)&( ~ ((uintptr_t)31)))), ((float * )((((uintptr_t)U_0_1)+31)&( ~ ((uintptr_t)31)))), 3.0, 0.05555555555555555, exp(1.0), exp(-0.05555555555555555), exp(-0.1111111111111111), exp(-0.2222222222222222), exp(-0.2777777777777778), exp(-0.4444444444444444), (1.0/(exp(1.0)+((((exp(-0.05555555555555555)+exp(-0.1111111111111111))+(exp(-0.2222222222222222)+exp(-0.4444444444444444)))*4.0)+(exp(-0.2777777777777778)*8.0)))), WIDTH, HEIGHT, cb_x, cb_y, chunk, _unroll_p3);
			// <--
			
		}
		// write_grids -->
		write_data_f("U_0_1_out.1.data", 2, U_0_1_out, WIDTH, HEIGHT);
		// <--
		
	}
	
	// validate
	if (1)
	{
#pragma omp parallel
		{
			// initialize_grids -->
			initialize_blur_parm(((float * )((((uintptr_t)U_0_0)+31)&( ~ ((uintptr_t)31)))), ((float * )((((uintptr_t)U_0_1)+31)&( ~ ((uintptr_t)31)))), 3.0, 0.05555555555555555, exp(1.0), exp(-0.05555555555555555), exp(-0.1111111111111111), exp(-0.2222222222222222), exp(-0.2777777777777778), exp(-0.4444444444444444), (1.0/(exp(1.0)+((((exp(-0.05555555555555555)+exp(-0.1111111111111111))+(exp(-0.2222222222222222)+exp(-0.4444444444444444)))*4.0)+(exp(-0.2777777777777778)*8.0)))), WIDTH, HEIGHT, cb_x, cb_y, chunk);
			initialize_blur_parm(U_0_0_ref, U_0_1_ref, 3.0, 0.05555555555555555, exp(1.0), exp(-0.05555555555555555), exp(-0.1111111111111111), exp(-0.2222222222222222), exp(-0.2777777777777778), exp(-0.4444444444444444), (1.0/(exp(1.0)+((((exp(-0.05555555555555555)+exp(-0.1111111111111111))+(exp(-0.2222222222222222)+exp(-0.4444444444444444)))*4.0)+(exp(-0.2777777777777778)*8.0)))), WIDTH, HEIGHT, cb_x, cb_y, chunk);
			// <--
			
#pragma omp barrier
			// compute_stencil -->
			blur_parm(( & U_0_1_out), ((float * )((((uintptr_t)U_0_0)+31)&( ~ ((uintptr_t)31)))), ((float * )((((uintptr_t)U_0_1)+31)&( ~ ((uintptr_t)31)))), 3.0, 0.05555555555555555, exp(1.0), exp(-0.05555555555555555), exp(-0.1111111111111111), exp(-0.2222222222222222), exp(-0.2777777777777778), exp(-0.4444444444444444), (1.0/(exp(1.0)+((((exp(-0.05555555555555555)+exp(-0.1111111111111111))+(exp(-0.2222222222222222)+exp(-0.4444444444444444)))*4.0)+(exp(-0.2777777777777778)*8.0)))), WIDTH, HEIGHT, cb_x, cb_y, chunk, _unroll_p3);
			// <--
			
		}
		// validate_computation -->
		int bHasErrors = 0;
		int _idx0;
		int _idx10;
		int _idx11;
		int _idx12;
		int _idx13;
		int _idx14;
		int _idx15;
		int _idx16;
		int _idx17;
		int _idx18;
		int _idx19;
		int _idx2;
		int _idx20;
		int _idx21;
		int _idx22;
		int _idx23;
		int _idx24;
		int _idx25;
		int _idx3;
		int _idx4;
		int _idx5;
		int _idx6;
		int _idx7;
		int _idx8;
		int _idx9;
		int p3_loc_idx_x;
		int p3_loc_idx_y;
		int t_ref;
		float *  tmp_swap_2;
		{
			/*
			for t_ref = 1..50 by 1 parallel 1 <level 0> schedule 1 { ... }
			*/
			for (t_ref=1; t_ref<=50; t_ref+=1)
			{
				/*
				for POINT pt_ref[t=t][0] of size [1, 1] in u0[t=t][0] + [ min=[0, 0], max=[0, 0] ] parallel 1 <level 0> schedule default { ... }
				*/
				{
					/* Index bounds calculations for iterators in pt_ref[t=t][0] */
					for (p3_loc_idx_y=2; p3_loc_idx_y<((HEIGHT-3)+1); p3_loc_idx_y+=1)
					{
						for (p3_loc_idx_x=2; p3_loc_idx_x<((WIDTH-3)+1); p3_loc_idx_x+=1)
						{
							/*
							pt_ref[t=t][0]=stencil(pt_ref[t=t][0])
							*/
							/* _idx0 = ((WIDTH*p3_loc_idx_y)+p3_loc_idx_x) */
							_idx0=((WIDTH*p3_loc_idx_y)+p3_loc_idx_x);
							/* _idx2 = ((WIDTH*p3_loc_idx_y)+(p3_loc_idx_x-1)) */
							_idx2=(_idx0-1);
							/* _idx3 = ((WIDTH*p3_loc_idx_y)+(p3_loc_idx_x+1)) */
							_idx3=(_idx2+2);
							/* _idx4 = ((WIDTH*(p3_loc_idx_y-1))+p3_loc_idx_x) */
							_idx4=((( - WIDTH)+_idx2)+1);
							/* _idx5 = ((WIDTH*(p3_loc_idx_y+1))+p3_loc_idx_x) */
							_idx5=((2*WIDTH)+_idx4);
							/* _idx6 = ((WIDTH*(p3_loc_idx_y-1))+(p3_loc_idx_x-1)) */
							_idx6=(_idx4-1);
							/* _idx7 = ((WIDTH*(p3_loc_idx_y-1))+(p3_loc_idx_x+1)) */
							_idx7=(_idx4+1);
							/* _idx8 = ((WIDTH*(p3_loc_idx_y+1))+(p3_loc_idx_x-1)) */
							_idx8=(_idx5-1);
							/* _idx9 = ((WIDTH*(p3_loc_idx_y+1))+(p3_loc_idx_x+1)) */
							_idx9=(_idx5+1);
							/* _idx10 = ((WIDTH*p3_loc_idx_y)+(p3_loc_idx_x-2)) */
							_idx10=((WIDTH+_idx4)-2);
							/* _idx11 = ((WIDTH*p3_loc_idx_y)+(p3_loc_idx_x+2)) */
							_idx11=((( - WIDTH)+_idx5)+2);
							/* _idx12 = ((WIDTH*(p3_loc_idx_y-2))+p3_loc_idx_x) */
							_idx12=(_idx4-WIDTH);
							/* _idx13 = ((WIDTH*(p3_loc_idx_y+2))+p3_loc_idx_x) */
							_idx13=(WIDTH+_idx5);
							/* _idx14 = ((WIDTH*(p3_loc_idx_y-1))+(p3_loc_idx_x-2)) */
							_idx14=(_idx4-2);
							/* _idx15 = ((WIDTH*(p3_loc_idx_y-2))+(p3_loc_idx_x-1)) */
							_idx15=(_idx12-1);
							/* _idx16 = ((WIDTH*(p3_loc_idx_y-2))+(p3_loc_idx_x+1)) */
							_idx16=((( - WIDTH)+_idx4)+1);
							/* _idx17 = ((WIDTH*(p3_loc_idx_y-1))+(p3_loc_idx_x+2)) */
							_idx17=(_idx4+2);
							/* _idx18 = ((WIDTH*(p3_loc_idx_y+1))+(p3_loc_idx_x-2)) */
							_idx18=(_idx5-2);
							/* _idx19 = ((WIDTH*(p3_loc_idx_y+2))+(p3_loc_idx_x-1)) */
							_idx19=((WIDTH+_idx5)-1);
							/* _idx20 = ((WIDTH*(p3_loc_idx_y+2))+(p3_loc_idx_x+1)) */
							_idx20=((WIDTH+_idx5)+1);
							/* _idx21 = ((WIDTH*(p3_loc_idx_y+1))+(p3_loc_idx_x+2)) */
							_idx21=(_idx5+2);
							/* _idx22 = ((WIDTH*(p3_loc_idx_y-2))+(p3_loc_idx_x-2)) */
							_idx22=(_idx12-2);
							/* _idx23 = ((WIDTH*(p3_loc_idx_y-2))+(p3_loc_idx_x+2)) */
							_idx23=((( - WIDTH)+_idx4)+2);
							/* _idx24 = ((WIDTH*(p3_loc_idx_y+2))+(p3_loc_idx_x-2)) */
							_idx24=((WIDTH+_idx5)-2);
							/* _idx25 = ((WIDTH*(p3_loc_idx_y+2))+(p3_loc_idx_x+2)) */
							_idx25=((WIDTH+_idx5)+2);
							U_0_1_ref[_idx0]=((1.0/(exp(1.0)+((((exp(-0.05555555555555555)+exp(-0.1111111111111111))+(exp(-0.2222222222222222)+exp(-0.4444444444444444)))*4.0)+(exp(-0.2777777777777778)*8.0))))*(((exp(1.0)*U_0_0_ref[_idx0])+((exp(-0.05555555555555555)*((U_0_0_ref[_idx2]+U_0_0_ref[_idx3])+(U_0_0_ref[_idx4]+U_0_0_ref[_idx5])))+(exp(-0.1111111111111111)*((U_0_0_ref[_idx6]+U_0_0_ref[_idx7])+(U_0_0_ref[_idx8]+U_0_0_ref[_idx9])))))+((exp(-0.2222222222222222)*((U_0_0_ref[_idx10]+U_0_0_ref[_idx11])+(U_0_0_ref[_idx12]+U_0_0_ref[_idx13])))+((exp(-0.2777777777777778)*(((U_0_0_ref[_idx14]+U_0_0_ref[_idx15])+(U_0_0_ref[_idx16]+U_0_0_ref[_idx17]))+((U_0_0_ref[_idx18]+U_0_0_ref[_idx19])+(U_0_0_ref[_idx20]+U_0_0_ref[_idx21]))))+(exp(-0.4444444444444444)*((U_0_0_ref[_idx22]+U_0_0_ref[_idx23])+(U_0_0_ref[_idx24]+U_0_0_ref[_idx25])))))));
						}
					}
				}
				U_0_1_out_ref=U_0_1_ref;
				tmp_swap_2=U_0_0_ref;
				U_0_0_ref=U_0_1_ref;
				U_0_1_ref=tmp_swap_2;
#pragma omp barrier 
			}
		}
		/*
		for POINT pt_ref[t=t][0] of size [1, 1] in u0[t=t][0] + [ min=[0, 0], max=[0, 0] ] parallel 1 <level 0> schedule default { ... }
		*/
		{
			/* Index bounds calculations for iterators in pt_ref[t=t][0] */
			for (p3_loc_idx_y=2; p3_loc_idx_y<((HEIGHT-3)+1); p3_loc_idx_y+=1)
			{
				for (p3_loc_idx_x=2; p3_loc_idx_x<((WIDTH-3)+1); p3_loc_idx_x+=1)
				{
					/* _idx0 = ((WIDTH*p3_loc_idx_y)+p3_loc_idx_x) */
					_idx0=((WIDTH*p3_loc_idx_y)+p3_loc_idx_x);
					if ((fabs(((U_0_1_out[_idx0]-U_0_1_out_ref[_idx0])/U_0_1_out_ref[_idx0]))>1.0E-4))
					{
						bHasErrors=1;
						break;
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
	free(U_0_0);
	free(U_0_0_ref);
	free(U_0_1);
	free(U_0_1_ref);
	// <--
	
	
	return 0;
}
