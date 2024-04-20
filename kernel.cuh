#ifndef KERNEL_CUH
#define KERNEL_CUH

#include <cuda_runtime.h>

__global__ void localSearch(int *targetBuffer_d, int *solutionBuffer_d, float *solutionBufferValue_d, int *permutationPool_d, float *W_d, int numBits, int numPermutations, int numSegments);

#endif