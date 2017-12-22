#include <iostream>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "definitions.cuh"
#include <time.h>
#include <stdio.h>
#include "util.h"
#include "dim_input.h"

//Number of elements on which to perform CFD
unsigned int Ni = PAR_I; //512; // Y elements
unsigned int Nj = PAR_J; //512; // X elements
unsigned int nIterations = PAR_ITER; //10000; // No Of Iterations

int main( int argc, char** argv ) {

    //Variables for graph
	char metric[] = "ms";
	
	//Variables for Timing
	float gpuTime;

	// CPU and GPU Pointers ( d_XX : refers to pointer pointing to GPU memory. This is just a convention)
	float *t = NULL, *t_prev = NULL;
	float *d_t = NULL,*d_t_prev= NULL;

    if ( ( Ni % THREADS_PER_BLOCK_Y != 0 )  || ( Nj % THREADS_PER_BLOCK_X != 0 ) ) {
        fprintf( stderr, "Please specify Ni & Nj as multiple of 16 !!!!" );
        exit( 0 );
    }

	printf("\n Ni= %d, Nj=%d nIteration=%d",Ni,Nj,nIterations);
	
	//unsigned int size = Ni * Nj * sizeof(float);

	if ( !initializeCPU( &t, &t_prev ) ) {
		printf( "\n Error in allocating memory on CPU!!!" );
		unInitializeCPU( &t, &t_prev );
		return 0;
	}

	if ( !initializeGPU( &d_t, &d_t_prev ) ) {
		printf( "\n Error in allocating memory on GPU!!!" );
		unInitializeCPU( &t, &t_prev );
		unInitializeGPU( &d_t, &d_t_prev );
		return 0;
	}

	//Perform CFD on GPU
	if ( !performGPUCFD( d_t,d_t_prev, t, t_prev, &gpuTime ) ) {
		printf( "\n GPU Kernel failed !!!" );
		unInitializeCPU( &t, &t_prev );
		unInitializeGPU( &d_t, &d_t_prev );
		return 0;
	}

    /* information read by the tool for printing graphs*/
    printMetric( metric, gpuTime );
    
	unInitializeCPU( &t, &t_prev );
	unInitializeGPU( &d_t, &d_t_prev );
	
	return EXIT_SUCCESS;

}

int initializeCPU( float **t, float **t_prev ) {
	*t = (float*) calloc( Ni*Nj, sizeof( float ) );
	*t_prev = (float*) calloc( Ni*Nj, sizeof( float ) );

	if ( (*t) == NULL || (*t_prev) == NULL )
		return 0;
	else
		return 1;
}

void unInitializeCPU( float **t, float **t_prev ) {
	if ( (*t) != NULL )
		free( *t );
	if ( (*t_prev) != NULL )
		free( *t_prev );
}

int initializeGPU( float **d_t, float **d_t_prev ) {

	unsigned int size = Ni * Nj * sizeof(float);

	// Choose which GPU to run on, change this on a multi-GPU system.
    cudaError_t cudaStatus = cudaSetDevice(0);
    if ( cudaStatus != cudaSuccess ) {
        fprintf( stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?" );
        return 0;
    }
	 // Allocate GPU buffers.
    cudaStatus = cudaMalloc( (void**) &(*d_t), size );
    if ( cudaStatus != cudaSuccess ) {
        fprintf( stderr, "cudaMalloc failed!" );
        return 0;
    }

	 // Allocate GPU buffers   .
    cudaStatus = cudaMalloc( (void**) &(*d_t_prev), size );
    if ( cudaStatus != cudaSuccess ) {
        fprintf( stderr, "cudaMalloc failed!" );
        return 0;
    }

	 // Memset GPU buffers
    cudaStatus = cudaMemset( (*d_t), 0, size );
    if ( cudaStatus != cudaSuccess ) {
        fprintf( stderr, "cudaMemset failed!" );
        return 0;
    }

	// Memset GPU buffers
    cudaStatus = cudaMemset( (*d_t_prev), 0, size );
    if ( cudaStatus != cudaSuccess ) {
        fprintf( stderr, "cudaMemset failed!" );
        return 0;
    }

	return 1;
}

void unInitializeGPU( float **d_t, float **d_t_prev ) {
	cudaError_t cudaStatus;

	if ( (*d_t) != NULL )
        cudaStatus = cudaFree( (*d_t) );
    if ( cudaStatus != cudaSuccess ) {
        fprintf( stderr, "cudaFree failed!" );
        return;
    }

	if ( (*d_t_prev) != NULL )
        cudaStatus = cudaFree( (*d_t_prev) );
    if ( cudaStatus != cudaSuccess ) {
        fprintf( stderr, "cudaFree failed!" );
        return;
    }

	cudaStatus = cudaDeviceReset();
    if ( cudaStatus != cudaSuccess ) {
        fprintf( stderr, "cudaDeviceReset failed!" );
        return;
    }
}

int performGPUCFD( float *d_t, float *d_t_prev, float *t, float *t_prev, float *gpuTime ) {

	float h,x,y;
	const char *str = (char*) malloc( 1024 ); // To store error string
	
	//Decide how many blocks per thread and how many blocks per grid
	dim3 dimBlock( THREADS_PER_BLOCK_X, THREADS_PER_BLOCK_Y );
	dim3 dimGrid( Nj/dimBlock.x, Ni/dimBlock.y );
	
	h = 1.0f/(Ni-1);
	memset( t_prev, 0, sizeof(float) * Ni * Nj );
  
	for ( unsigned int i=0; i<Ni; i++ ) {
		x = i*h;
		t_prev[i*Nj+0] = x*x;
		t_prev[i*Nj+(Nj-1)] = x*x + 1.0f;
	}

	for ( unsigned int j=0; j<Nj; j++ ) {
		y = j*h;
		t_prev[0*Nj+j] = y*y;
		t_prev[((Ni-1) * Nj) + j] = 1.0f + y*y;
	}

	//Copy data to device
	cudaMemcpy( d_t_prev, t_prev, sizeof(float) * Ni * Nj , cudaMemcpyHostToDevice );

	//Insert event to calculate time
	cudaEvent_t start, stop;
	cudaEventCreate( &start );
	cudaEventCreate( &stop );

	//This calls Version 1 of kernel which uses Global memory
	cudaEventRecord( start, 0 );
	
	for ( unsigned int k=0; k<nIterations; k++)	{
		// Launch a kernel on the GPU with one thread for each element.
		calculateCFD_V1<<<dimGrid,dimBlock>>>( d_t_prev,d_t, Ni, Nj, h );
		
		float* pingPong = d_t_prev;
		d_t_prev = d_t;
		d_t = pingPong;
	}
	cudaEventRecord( stop, 0 );
	cudaEventSynchronize( stop );

	float elapsedTime;
	cudaEventElapsedTime( &elapsedTime, start, stop );
	//printf( "\n GPU Time:: %f ms", elapsedTime );

	*gpuTime = elapsedTime;

	cudaError_t cudaStatus = cudaMemcpy( t, d_t_prev, sizeof(float) * Ni * Nj , cudaMemcpyDeviceToHost );
	if ( cudaStatus != cudaSuccess ) {
		fprintf( stderr, "cudaMemcpy failed!" );
		str = cudaGetErrorString( cudaStatus );
		fprintf( stderr, "CUDA Error!:: %s\n", str );
		return 0;
	}
	
	return 1;
}

