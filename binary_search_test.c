#include <stdio.h>
#include <stdlib.h>
#include "quicksort.h"
#include "binary_search_test.h"
#include <stdbool.h>
#define false 0
#define true 1



int binarySearch(float *arr, int l, int r, float x,int solutionPoolSize) {
    ///m here
    int m;
    while (l <= r) {
        m = l + (r - l) / 2;
        if (arr[m] < x) {
            l = m + 1;
        } else if (arr[m] > x) {
            r = m - 1;
        } else {
            // if the new value insert already appeared in the solutionpool_value, we will not insert it.    //still insert! need modication!
            return m+1;
        }
    }
    // if the newly inserted solution is larger than the worst value in the solutionpool, we will not insert it.
    if (l >= solutionPoolSize)
    {
        return -1;
    }
    else
    {
        return l; // In other cases, we return the insert position, the position range should in [0, solutionPoolSize -1]

    }
    
}


////The input of solutionpool and solutionPoolValue should already been sorted.
int binarySearchInsert(int **solutionPool, int **solutionBuffer_h, float *solutionPoolValue, float *solutionBufferValue_h, int solutionPoolSize, int numBlocks) {
    int ifInsert = false;
    // we iterate all element in the solution buffer and decide if we insert it.
    for (int i = 0; i < numBlocks; i++) 
    {
        float bufferValue = solutionBufferValue_h[i];
        int *bufferSolution = solutionBuffer_h[i];

        // Find the position to insert the new solution
        int position = binarySearch(solutionPoolValue, 0, solutionPoolSize - 1, bufferValue,solutionPoolSize);
        // If same value already appeared or the value is worse than all value in the solutionpool_value, we will not insert it. 


        // If there's a better solution to insert
        if (position != -1) 
        {
            ifInsert = true;
            if (position == solutionPoolSize - 1)
            {
                solutionPoolValue[position] = bufferValue;
                solutionPool[position] = bufferSolution;
            }
            else
            {
                // Shift worse solutions down by one to make room for the new one
                for (int j = solutionPoolSize - 1; j > position; j--) 
                {
                    solutionPoolValue[j] = solutionPoolValue[j - 1];
                    solutionPool[j] = solutionPool[j - 1];
                }

                // Insert the new solution
                solutionPoolValue[position] = bufferValue;
                solutionPool[position] = bufferSolution;
            }
        }
    }
    return ifInsert;
}



////When compling, running this command: "gcc -o binary_search_test binary_search_test.c quicksort.c"

//////-----------------------------------TEST CASE FOR BINARY SEARCH INSERT------------------------------------------------------------------------------ 
// int main(){
//     #define NUM_SOLUTIONS 8
//     #define NUM_GENES 10 
//     int **solutionpool = (int **)malloc(NUM_SOLUTIONS * sizeof(int *));
//     int **solutionbuffer = (int **)malloc(NUM_SOLUTIONS * sizeof(int *));

//     float *solutionpool_value  = (float *)malloc(NUM_SOLUTIONS * sizeof(float));
//     float *solutionbuffer_value  = (float *)malloc(NUM_SOLUTIONS * sizeof(float));

//     ////randomly generate solutionpool_value      
//     for (int i = 0; i<NUM_SOLUTIONS; i++ )
//     {
//         solutionpool_value[i] = rand() % 100;
//         solutionbuffer_value[i] = rand() % 100;

//     }


//     ////randomly generate solutionpool
//     for(int i = 0; i< NUM_SOLUTIONS; i++)
//     {
//         solutionpool[i] = (int *)malloc(NUM_GENES * sizeof(int));
//         solutionbuffer[i] = (int *)malloc(NUM_GENES * sizeof(int));
//         for (int j = 0; j < NUM_GENES; j++) {
//             solutionpool[i][j] = rand() % 2;  
//             solutionbuffer[i][j] = rand() % 2;  

//         }
//     }

//     ////------------------------------------------------------------------------------
//     quickSort(solutionpool, solutionpool_value, 0, NUM_SOLUTIONS-1);
//     binarySearchInsert(solutionpool, solutionbuffer, solutionpool_value, solutionbuffer_value, NUM_SOLUTIONS);


//  //Test if the solutionpool_value generated correctly
//     for (int i = 0; i < NUM_SOLUTIONS; i++) 
//     {
//         printf("the %dth solution pool value is: %f\n", i+1, solutionpool_value[i]);
//     }

//     for (int i = 0; i < NUM_SOLUTIONS; i++) 
//     {
//         printf("the %dth solution buffer value is: %f\n", i, solutionbuffer_value[i]);
//     }

//     ////Test if the solutionpool generated correctly
//     for (int i = 0; i < NUM_SOLUTIONS; i++)
//     {   
//         printf("the %dth chromosome in solution pool is: ", i+1);
//         for (int j = 0; j < NUM_GENES; j++)
//         {
//             printf("%d",solutionpool[i][j]);
//         }
//         printf("\n");
//     }

//     for (int i = 0; i < NUM_SOLUTIONS; i++)
//     {   
//         printf("the %dth chromosome in solution buffer is: ", i+1);
//         for (int j = 0; j < NUM_GENES; j++)
//         {
//             printf("%d",solutionbuffer[i][j]);
//         }
//         printf("\n");
//     }


    //////---------------------------------------TEST CASE FOR BINARY SEARCH----------------------------------------------------------------
    // int pos = binarySearch(solutionpool_value,0,NUM_SOLUTIONS -1, 15);
    // printf("The first position is %d",pos);
    // printf("\n");

    // int pos1 = binarySearch(solutionpool_value,0,NUM_SOLUTIONS -1, 100);
    // printf("The second position is %d",pos1);
    // printf("\n");





// }