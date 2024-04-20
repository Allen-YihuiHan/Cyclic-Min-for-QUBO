// binary_search_test.h
#ifndef BINARY_SEARCH_TEST_H
#define BINARY_SEARCH_TEST_H
#define false 0
#define true 1

#include "quicksort.h" // Include the quicksort header if binarySearchInsert uses functions from it

// Function prototypes for the functions defined in binary_search_test.c
int binarySearch(float *arr, int l, int r, float x,int solutionPoolSize);
int binarySearchInsert(int **solutionPool, int **solutionBuffer_h, float *solutionPoolValue, float *solutionBufferValue_h, int solutionPoolSize,int numBlocks);

#endif // BINARY_SEARCH_TEST_H
