// genetic_algorithm.h

#ifndef GENETIC_ALGORITHM_H
#define GENETIC_ALGORITHM_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

// Constants
#define MAX_ITER 300
#define MUTATION_RATE 0.05
#define CROSS_PROB 0.75
#define THRESHOLD 0.75
#define MAX_OCCUR 272

// Struct definitions
typedef struct
{
    int *genes;
    double fitness;
} Chromosome;

typedef struct HashNode
{
    char *key;
    int count;
    struct HashNode *next;
} HashNode;

// Function prototypes
Chromosome *create_chromosome(int numbits);
Chromosome **initialize_population(int numbits, int solutionPoolSize);
// void initialize_population_data(Chromosome **population, int **tempPool, int numbits, int solutionPoolSize);
unsigned int hashFunction(char *str);
void insertOrUpdate(char *key, HashNode **hashTable);
void free_population(Chromosome **population, int solutionPoolSize);
void freeHashTable(HashNode **hashTable);
double evaluate(Chromosome *chrom, int numbits, float **W);                                   // Ensure W is defined correctly or passed appropriately
int *count(Chromosome **population, HashNode **hashTable, int numbits, int solutionPoolSize); //, int occur[]);
void selection(Chromosome **population, int numbits, int solutionPoolSize, float **W);        // Ensure W is defined or passed correctly
void crossover(Chromosome **population, int numbits, int solutionPoolSize);
void mutation(Chromosome **population, int numbits, int solutionPoolSize);
void runGeneticAlgor(Chromosome **population, int numbits, int solutionPoolSize, float **W);

#endif // GENETIC_ALGORITHM_H
