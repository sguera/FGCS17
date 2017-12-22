#include "definitions.cuh"

//Performs CFD calculation on global memory. This code does not use any advance optimization technique on GPU
// But still acheives many fold performance gain
__global__ void calculateCFD_V1( float* input,  float* output, unsigned int Ni, unsigned int Nj, 
								   float h)
{
	unsigned int i = blockDim.x * blockIdx.x + threadIdx.x; // Y - ID
	unsigned int j = blockDim.y * blockIdx.y + threadIdx.y; // X - ID

	unsigned int iPrev = i-1; // Previous Y element
	unsigned int iNext = i+1; // Next Y element

	unsigned int jPrev = j-1; //Previous X element
	unsigned int jNext = j+1; // Next X element


	unsigned int index = i * Nj + j;

	if( i > 0 && j > 0 && i < (Ni-1) && j <(Nj-1))
		output[index] = 0.25f * (input[iPrev * Nj + j] + input[iNext* Nj + j] + input[i * Nj+ jPrev] 
			+ input[i* Nj + jNext] - 4*h*h);
}
