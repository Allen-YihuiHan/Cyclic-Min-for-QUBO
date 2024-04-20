#include <stdio.h>
#include <stdlib.h>
#include <time.h>
extern "C"
{
#include "binary_search_test.h"
#include "genetic_algorithm.h"
}
#include "kernel.cuh"

#define MAX_LOOP 10000

void freeMatrixFloat(float **matrix, int rows)
{
    for (int i = 0; i < rows; i++)
        free(matrix[i]);
    free(matrix);
}

void freeMatrixInt(int **matrix, int rows)
{
    for (int i = 0; i < rows; i++)
        free(matrix[i]);
    free(matrix);
}

// Use Fisher-Yates Algorithm to shuffle an array
void shuffle(int *list, int size)
{
    int pos, temp;

    for (int i = 0; i < size; i++)
        list[i] = i;

    for (int i = 0; i < size - 1; i++)
    {
        pos = i + rand() % (size - i);
        temp = list[i];
        list[i] = list[pos];
        list[pos] = temp;
    }
}

// Sort a matrix according to the provided order
// Specifically, turn W to P^(-1)WP
void sort(float **matrix, int *permutation, int size)
{
    // first get a copy of the original matrix
    float **copy = (float **)malloc(size * sizeof(float *));
    for (int i = 0; i < size; i++)
        copy[i] = (float *)malloc(size * sizeof(float));

    for (int i = 0; i < size; i++)
        for (int j = 0; j < size; j++)
            copy[i][j] = matrix[i][j];

    // do right multiplication
    for (int i = 0; i < size; i++)
        for (int j = 0; j < size; j++)
            matrix[permutation[i]][j] = copy[i][j];

    for (int i = 0; i < size; i++)
        for (int j = 0; j < size; j++)
            copy[i][j] = matrix[i][j];

    // do left multiplication
    for (int i = 0; i < size; i++)
        for (int j = 0; j < size; j++)
            matrix[i][permutation[j]] = copy[i][j];

    for (int i = 0; i < size; i++)
        free(copy[i]);

    free(copy);
}

float calculate(int *solution, float **W, int numBits)
{
    float result = 0;

    for (int i = 0; i < numBits; i++)
        for (int j = 0; j < numBits; j++)
            result = result + solution[i] * solution[j] * W[i][j];

    return result;
}

// insert the good solutions from solution buffer to solution pool, while maintaining the energy order
// If nothing can be inserted to solution buffer, and the energy difference within solution pool is small, establish convergence by returning true
int updateSolutionPool(int **solutionPool, int **solutionBuffer_h, float *solutionPoolValue, float *solutionBufferValue_h, int numBlocks, int solutionPoolSize, int numBits, float epsilon)
{
    int ifConvergence = false;
    int ifInsert = false;
    int pos = 0;

    for (int i = 0; i < numBlocks; i++)
    {
        if (solutionBufferValue_h[i] < solutionPoolValue[solutionPoolSize - 1])
        {
            pos = 0;
            ifInsert = true;
            while (solutionBufferValue_h[i] > solutionPoolValue[pos])
                pos++;

            for (int j = solutionPoolSize - 1; j > pos; j--)
            {
                solutionPoolValue[j] = solutionPoolValue[j - 1];
                for (int k = 0; k < numBits; k++)
                    solutionPool[j][k] = solutionPool[j - 1][k];
            }

            solutionPoolValue[pos] = solutionBufferValue_h[i];

            for (int k = 0; k < numBits; k++)
                solutionPool[pos][k] = solutionBuffer_h[i][k];
        }
    }

    if (ifInsert == false && solutionPoolValue[solutionPoolSize - 1] - solutionPoolValue[0] <= epsilon)
        ifConvergence = true;

    return ifConvergence;
}

int main(int argc, char **argv)
{
    clock_t begin, end;
    const unsigned int numBlocks = 20;
    int numPermutations = 5;   // This is the number of permutations used in device algorithms
    int solutionPoolSize = 20; // This is the number of solutions stored in host
    int numSegments = 8;       // The number of segments to break the solution into, in cyclic-min algorithm
    float epsilon = 1.0;
    int *permutationPool_d; // Flatten out the permutation arrays by array1 + array2+ ...
    float *W_d;             // Flatten out the weight matrices by W1row1 + W1row2 + ... + W2row1 + W2row2 + ...
    float optimalValue;
    int numBits;

    if (argc < 2)
    {
        fprintf(stderr, "Error: Too few arguments. Please provide at least one file.\n");
        exit(1);
    }

    for (int n_files = 1; n_files < argc; n_files++)
    {
        FILE *fp = fopen(argv[n_files], "r");

        if (fp == NULL)
        {
            fprintf(stderr, "Error opening file\n");
            return 1;
        }

        fscanf(file, "%f\n", &optimalValue);
        fscanf(file, "%d\n", &numBits);

        float **W = (float **)malloc(numBits * sizeof(float *));
        for (int i = 0; i < numBits; i++)
            W[i] = (float *)malloc(numBits * sizeof(float));

        for (i = 0; i < numBits; i++)
            for (j = 0; j < numBits; j++)
                fscanf(file, "%f\n", &W[i][j]);

        fclose(fp);

        // Define the Objects and allocate host&device space

        int **solutionPool = (int **)malloc(solutionPoolSize * sizeof(int *));
        for (int i = 0; i < solutionPoolSize; i++)
            solutionPool[i] = (int *)calloc(numBits, sizeof(int)); // initialize the solution pool with all zero using calloc.

        float *solutionPoolValue = (float *)calloc(solutionPoolSize, sizeof(float *));

        // Allocate memory for permutation information
        int **permutationPool_h = (int **)malloc(numPermutations * sizeof(int *));

        // Here I use triple pointers for weight matrices corresponding to different permutations, where it can be simplified
        float ***W_h = (float ***)malloc(numPermutations * sizeof(float **));

        //// allocate memory
        for (int i = 0; i < numPermutations; i++)
        {
            permutationPool_h[i] = (int *)malloc(numBits * sizeof(int));

            W_h[i] = (float **)malloc(numBits * sizeof(float *));

            for (int j = 0; j < numBits; j++)
                W_h[i][j] = (float *)malloc(numBits * sizeof(float));
        }

        for (int i = 0; i < numPermutations; i++)
            for (int j = 0; j < numBits; j++)
                for (int k = 0; k < numBits; k++)
                    W_h[i][j][k] = W[j][k];

        ////first dimension in cpu, second dimension in gpu???
        cudaMalloc((void **)&(permutationPool_d), numPermutations * numBits * sizeof(int));

        cudaMalloc((void **)&(W_d), numPermutations * numBits * numBits * sizeof(float *));

        cudaDeviceSynchronize();

        // Copy the initialized host permutation information to their correspondence in global memory
        for (int i = 0; i < numPermutations; i++)
        {
            shuffle(permutationPool_h[i], numBits); // fisher-yates shuffle the array

            sort(W_h[i], permutationPool_h[i], numBits);

            cudaMemcpy(permutationPool_d + i * numBits, permutationPool_h[i], sizeof(int) * numBits, cudaMemcpyHostToDevice);

            for (int j = 0; j < numBits; j++)
                cudaMemcpy(W_d + (i * numBits * numBits) + (j * numBits), W_h[i][j], sizeof(float) * numBits, cudaMemcpyHostToDevice);
        }

        cudaDeviceSynchronize();

        // Define the number of blocks to be used, assume the maximum number of threads in each block is 512
        // if number of bits of a solution exceeds 512, it implies that each thread will handle more than one bit.
        // --------------------------------------------------------------------------------------------------------
        // WARNING: The current implementation does not cosider numBits > 512, waiting for future implementation
        // --------------------------------------------------------------------------------------------------------

        const unsigned int THREADS_PER_BLOCK = numBits;

        dim3 gridDim(numBlocks, 1, 1), blockDim(THREADS_PER_BLOCK, 1, 1);

        int **targetBuffer_h = (int **)malloc(solutionPoolSize * sizeof(int *));
        int **solutionBuffer_h = (int **)malloc(numBlocks * sizeof(int *));
        float *solutionBufferValue_h = (float *)malloc(numBlocks * sizeof(float));

        int *targetBuffer_d;
        int *solutionBuffer_d;
        float *solutionBufferValue_d;

        cudaDeviceSynchronize();

        for (int i = 0; i < numBlocks; i++)
        {
            targetBuffer_h[i] = (int *)malloc(numBits * sizeof(int));
            solutionBuffer_h[i] = (int *)malloc(numBits * sizeof(int));
        }

        //----------------------------------------------------------------------------------------------------------------
        // GPU kernels are asynchronous with host by default
        // to allow more concurrency, use different streams, where the kernel methods run by different streams may overlap
        // The paper says using cudaMemcpyAsync, and "repreatedly reads" a counter
        // The problem lies in how to create a task buffer, where the cuda blocks continues another task automatically after finishing one
        // For now, we will implement the version without asychronous kernels
        // Initialize solution pool with genetic algorithm, while sorting them according to energy value
        // initialize target buffer with genetic algorithm
        // Use tempPool, first copy the value of solutionPool and perform GA, then choose the first numBlocks chromsomes give to TargetBuffer_h, then TargetBuffer_d, and perfrom the next iteration's local search.
        // This also means that the GA will take a tempPool as input and output a same size tempPool.

        int convergence = false;
        int iteration = MAX_LOOP;
        // The MAIN LOOP:

        // Create a population called tempPool_chrom here to avoid repeatly creation.
        Chromosome **tempPool_chrom = (Chromosome **)malloc(solutionPoolSize * sizeof(Chromosome *));
        for (int i = 0; i < solutionPoolSize; i++)
        {
            tempPool_chrom[i] = (Chromosome *)malloc(sizeof(Chromosome));
            tempPool_chrom[i]->genes = (int *)malloc(numBits * sizeof(int));
            tempPool_chrom[i]->fitness = 0.0;
        }

        begin = clock();
        while (iteration--)
        {
            // Copy the solutionPool to a tempPool and perfrom GA on tempPool.
            for (int i = 0; i < solutionPoolSize; i++)
                for (int j = 0; j < numBits; j++)
                    tempPool_chrom[i]->genes[j] = solutionPool[i][j];

            //--------------------------Perfrom GA!-----------------------------------------------------
            runGeneticAlgor(tempPool_chrom, numBits, solutionPoolSize, W);
            // select the first numBlocks elements to targetBuffer_h:

            for (int i = 0; i < numBlocks && i < solutionPoolSize; i++)
                for (int j = 0; j < numBits; j++)
                    targetBuffer_h[i][j] = tempPool_chrom[i]->genes[j];

            // copy the value from the targetBuffer_h to the targetBuffer_d:
            for (int i = 0; i < numBlocks; i++)
            {
                err = cudaMemcpy(targetBuffer_d + i * numBits, targetBuffer_h[i], numBits * sizeof(int), cudaMemcpyHostToDevice);
                if (err != cudaSuccess)
                    fprintf(stderr, "CUDA Error: %s\n", cudaGetErrorString(err));
            }

            // call the local search, if converge, then we break out of the loop. if not, we reuse the genetic algorithm.
            localSearch<<<gridDim, blockDim>>>(targetBuffer_d, solutionBuffer_d, solutionBufferValue_d, permutationPool_d, W_d, numBits, numPermutations, numSegments);

            err = cudaGetLastError();
            if (err != cudaSuccess)
                fprintf(stderr, "CUDA Kernel Launch Error: %s\n", cudaGetErrorString(err));

            cudaDeviceSynchronize();
            err = cudaMemcpy(solutionBufferValue_h, solutionBufferValue_d, sizeof(float) * numBlocks, cudaMemcpyDeviceToHost);

            if (err != cudaSuccess)
                fprintf(stderr, "CUDA Error: %s\n", cudaGetErrorString(err));

            for (int i = 0; i < numBlocks; i++)
            {
                cudaMemcpy(solutionBuffer_h[i], solutionBuffer_d + i * numBits, sizeof(int) * numBits, cudaMemcpyDeviceToHost);
                cudaDeviceSynchronize();
            }

            cudaDeviceSynchronize();
            // update solution pool by inserting the better solutions from solutionBuffer to solutionPool
            convergence = updateSolutionPool(solutionPool, solutionBuffer_h, solutionPoolValue, solutionBufferValue_h, numBlocks, solutionPoolSize, numBits, epsilon);

            // Use genetic algorithm to generate another round of new target solutions
            // And then update targetBuffer_h, targetBuffer_d

            if (iteration % 2 == 0)
            {
                printf("\nIteration %d\n", iteration);
                printf("Optimal Value is %f \n", optimalValue);

                for (int i = 0; i < 1; i++)
                    printf("bestValue: %f\n", solutionPoolValue[i]);

                if (solutionPoolValue[0] < optimalValue + 2)
                {
                    printf("Convergence Criterion Satisfied.\n");
                    break;
                }

                printf("--------------------------------------------------------------------");
            }
        }
        end = clock();

        free(solutionPoolValue);
        free(solutionBufferValue_h);
        freeMatrixInt(targetBuffer_h, numBlocks);
        freeMatrixInt(solutionBuffer_h, numBlocks);
        freeMatrixInt(permutationPool_h, numPermutations);
        freeMatrixInt(solutionPool, solutionPoolSize);
        freeMatrixFloat(W, numBits);
        for (int i = 0; i < numPermutations; i++)
            freeMatrixFloat(W_h[i], numBits);
        cudaFree(W_d);
        cudaFree(permutationPool_d);
        cudaFree(targetBuffer_d);
        cudaFree(solutionBuffer_d);
        cudaFree(solutionBufferValue_d);
        free_population(tempPool_chrom, solutionPoolSize);

        printf("\nTotal Time Spent is %f seconds", (float)(end - begin) / CLOCKS_PER_SEC);
    }

    return 0;
}