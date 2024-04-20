#include <float.h>
#include <stdio.h>
#include <time.h>
#include "kernel.cuh"

#define true 1
#define false 0

__global__ void localSearch(int *targetBuffer_d, int *solutionBuffer_d, float *solutionBufferValue_d, int *permutationPool_d, float *W_d, int numBits, int numPermutations, int numSegments)
{
    // Energy difference delta_i declared in shared memory (see definition of delta_i in paper)
    // Array size is static, either set an upper bound, or use extern keyword to dynamically specify size (where there should be modifications on calling the kernel function)
    // We use static array size for now, can change later

    // Record the current solution, best solution, and their corresponding values in shared memory, so that all the threads can see them
    __shared__ float delta[512];
    __shared__ int currentSolution[512]; // Initialize the current solution to be a vector of zeros
    __shared__ int bestSolution[512];
    __shared__ float min_delta;
    __shared__ float currentValue;
    __shared__ float bestValue;

    // NOTE: shared variable is not allowed to initialize upon declaration
    //  newly added for straight search and cyclic-min
    __shared__ int k;
    __shared__ int flag;
    __shared__ int numSegBits;
    // Record target solution we want to approach
    // And do bit permutation, the permutation used is blockIdx.x mod numPermutations
    __shared__ int targetSolution[512];

    // Auxiliary variables in the reduction algorithm
    __shared__ float aux_value[512];
    __shared__ unsigned int aux_index[512];

    int pIndex = blockIdx.x % numPermutations;

    // Permutation & initilization
    if (threadIdx.x < numBits)
    {
        targetSolution[permutationPool_d[pIndex * numBits + threadIdx.x]] = targetBuffer_d[blockIdx.x * numBits + threadIdx.x];
        delta[threadIdx.x] = W_d[pIndex * numBits * numBits + threadIdx.x * numBits + threadIdx.x];

        // initially starting at 0, the flip will change the i-th element by W[i][i]
        currentSolution[threadIdx.x] = 0;
        bestSolution[threadIdx.x] = 0;
    }
    __syncthreads();

    if (threadIdx.x == 0)
    {
        currentValue = 0;
        min_delta = FLT_MAX;
        bestValue = 0;
        k = 0;
        flag = numBits;
        numSegBits = numBits - 1 / numSegments + 1;
    }
    __syncthreads();

    // REMARK: the weight matrix used here is W_d[pIndex] (it is a pointer to pointer of float)
    // One question needs discussion: should we copy the permutated weight matrix on shared memory of this block? The matrix is large

    // Implement straight search (Exactly the same as Algorithm 3 in paper)
    while (flag != 0) // flag is set to be numBits, when every elem of array does not change, we get T from 0.
    {
        flag = numBits; // reset flag
        min_delta = FLT_MAX;

        __syncthreads();

        // This step is to find the smallest among energy difference "delta"
        // Use Reduction algorithm
        if (currentSolution[threadIdx.x] != targetSolution[threadIdx.x])
        {
            aux_index[threadIdx.x] = threadIdx.x;
            aux_value[threadIdx.x] = delta[threadIdx.x];
        }
        else
        {
            aux_index[threadIdx.x] = threadIdx.x;
            aux_value[threadIdx.x] = FLT_MAX;
        }
        __syncthreads();

        for (int s = numBits / 2; s > 0; s >>= 1)
        {
            if (threadIdx.x < s)
            {
                if (aux_value[threadIdx.x] > aux_value[threadIdx.x + s])
                {
                    aux_value[threadIdx.x] = aux_value[threadIdx.x + s];
                    aux_index[threadIdx.x] = aux_index[threadIdx.x + s];
                }
            }
            __syncthreads();
        }

        if (threadIdx.x == 0)
        {
            k = aux_index[0];
            min_delta = aux_value[0];
        }
        __syncthreads();

        // update each delta i
        if (threadIdx.x == k)
            delta[threadIdx.x] *= -1;
        else
            delta[threadIdx.x] += 2 * W_d[pIndex * numBits * numBits + threadIdx.x * numBits + k] * (-2 * currentSolution[threadIdx.x] + 1) * (-2 * currentSolution[k] + 1);

        // Update the best solution & value
        if (threadIdx.x == k)
        {
            currentSolution[k] = 1 - currentSolution[k]; // flip
            currentValue += min_delta;                   //// current energy function
        }
        __syncthreads();

        if (currentValue < bestValue)
        {
            bestSolution[threadIdx.x] = currentSolution[threadIdx.x];
        }
        __syncthreads();

        // Do not change this. Handle condition related to shared variable with care
        if (currentValue < bestValue && threadIdx.x == k)
            bestValue = currentValue;
        __syncthreads();

        // judge whether to get out of while-loop
        if (currentSolution[threadIdx.x] == targetSolution[threadIdx.x])
            atomicAdd(&flag, -1); // need atomic operations
        __syncthreads();
    }

    // Now that we know energy value of the target solution,
    // Implement cyclic-Min
    for (int i = 0; i < numBits; i += numSegBits)
    {
        if (threadIdx.x >= i && threadIdx.x < i + numSegBits)
        {
            aux_index[threadIdx.x] = threadIdx.x;
            aux_value[threadIdx.x] = delta[threadIdx.x];
        }
        else
        {
            aux_index[threadIdx.x] = threadIdx.x;
            aux_value[threadIdx.x] = FLT_MAX;
        }
        __syncthreads();

        for (int s = numBits / 2; s > 0; s >>= 1)
        {
            if (threadIdx.x < s)
            {
                if (aux_value[threadIdx.x] > aux_value[threadIdx.x + s])
                {
                    aux_value[threadIdx.x] = aux_value[threadIdx.x + s];
                    aux_index[threadIdx.x] = aux_index[threadIdx.x + s];
                }
            }
            __syncthreads();
        }

        if (threadIdx.x == 0)
        {
            k = aux_index[0];
            min_delta = aux_value[0];
        }
        __syncthreads();

        if (threadIdx.x == k)
            delta[threadIdx.x] *= -1;
        else
            delta[threadIdx.x] += 2 * W_d[pIndex * numBits * numBits + numBits * threadIdx.x + k] * (-2 * currentSolution[threadIdx.x] + 1) * (-2 * currentSolution[k] + 1);

        __syncthreads(); //// wait for all the delta_i been calculated and find a minimum delta.

        // Update the best solution & value
        if (threadIdx.x == k)
        {
            currentSolution[k] = 1 - currentSolution[k]; // flip
            currentValue += min_delta;                   //// current energy function
        }
        __syncthreads();

        if (currentValue < bestValue)
            bestSolution[threadIdx.x] = currentSolution[threadIdx.x];
        __syncthreads();

        if (currentValue < bestValue && threadIdx.x == k)
            bestValue = currentValue;

        __syncthreads();
    }

    // Reverse back the bestSolution, and write it on solution buffer
    if (threadIdx.x < numBits)
        solutionBuffer_d[blockIdx.x * numBits + threadIdx.x] = bestSolution[permutationPool_d[pIndex * numBits + threadIdx.x]];

    solutionBufferValue_d[blockIdx.x] = bestValue;
}