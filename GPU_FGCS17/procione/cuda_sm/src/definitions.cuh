//Number of Threads per bolock
#define THREADS_PER_BLOCK_X 32
#define THREADS_PER_BLOCK_Y 16
#define RADIUS 1

int initializeCPU( float **t, float **t_prev );
int initializeGPU( float **d_t, float **d_t_prev );

void unInitializeCPU( float **t, float **t_prev );
void unInitializeGPU( float **d_t, float **d_t_prev );

int performGPUCFD( float *d_t, float *d_t_prev, float *t, float *t_prev, float *gpuTime );

//Kernel Declaration
__global__ void calculateCFD_V2( float* input,  float* output, unsigned int Ni, unsigned int Nj, float h );


