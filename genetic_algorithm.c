#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "genetic_algorithm.h"

#define MAX_ITER 300       // maximum generation (iteration)
#define MUTATION_RATE 0.05 // rate of mutation
#define CROSS_PROB 0.75    // probability of crossover
#define THRESHOLD 0.6      // stopping propoortion
#define MAX_OCCUR 272      // max key for each chromosome
#define TABLE_SIZE 1000    // size of hash table

// dynamically create a single chromsome.

Chromosome *create_chromosome(int numbits)
{
    Chromosome *chrom = (Chromosome *)malloc(sizeof(Chromosome));
    chrom->genes = (int *)malloc(numbits * sizeof(int));
    for (int j = 0; j < numbits; j++)
        chrom->genes[j] = rand() % 2; // 0 or 1
    chrom->fitness = 0.0;

    return chrom;
}

// create a population of chromesomes, return a population.

Chromosome **initialize_population(int numbits, int solutionPoolSize)
{
    // every population is a array of pointers, each pointers point to a struct, which is a single Chromsome.
    Chromosome **population = (Chromosome **)malloc(solutionPoolSize * sizeof(Chromosome *));
    for (int i = 0; i < solutionPoolSize; i++)
    {
        population[i] = create_chromosome(numbits);
    }
    return population;
}

unsigned int hashFunction(char *str)
{
    unsigned long hash = 5381;
    int c;
    while ((c = *str++))
        hash = ((hash << 5) + hash) + c;
    return hash % TABLE_SIZE;
}

void insertOrUpdate(char *key, HashNode **hashTable)
{
    unsigned int index = hashFunction(key);
    HashNode *node = hashTable[index];

    while (node != NULL)
    {
        if (strcmp(node->key, key) == 0)
        {
            node->count++;
            return;
        }
        node = node->next;
    }

    HashNode *newNode = (HashNode *)malloc(sizeof(HashNode));
    newNode->key = strdup(key);
    newNode->count = 1;
    newNode->next = hashTable[index];
    hashTable[index] = newNode;
}

// release the memroy that a population.

void free_population(Chromosome **population, int solutionPoolSize)
{
    for (int i = 0; i < solutionPoolSize; i++)
    {
        if (population[i] != NULL)
        {
            free(population[i]->genes); // free the genes array memory
            free(population[i]);        // free the struct memory
        }
    }
    free(population); // free the whole population memory
}

void freeHashTable(HashNode **hashTable)
{
    HashNode *current, *temp;
    for (int i = 0; i < TABLE_SIZE; ++i)
    {
        current = hashTable[i];
        while (current != NULL)
        {
            temp = current;
            current = current->next;
            free(temp->key);
            free(temp);
        }
    }
}

// x^TWx

float evaluate(Chromosome *chrom, int numbits, float **W) // modification, W need to be input or not?
{
    float value = 0.0;
    for (int i = 0; i < numbits; i++)
        for (int j = 0; j < numbits; j++)
            value += chrom->genes[i] * chrom->genes[j] * W[i][j];
    return value;
}

int *count(Chromosome **population, HashNode **hashTable, int numbits, int solutionPoolSize) //,int occur[])  //MAX_OCCUR need to be defined in advance
{
    char *geneKey = (char *)malloc((numbits + 1) * sizeof(char));
    int *flag = (int *)malloc(2 * sizeof(int));
    flag[0] = 0;

    for (int i = 0; i < solutionPoolSize; i++)
    {
        for (int j = 0; j < numbits; j++)
            geneKey[j] = '0' + population[i]->genes[j];

        geneKey[numbits] = '\0';
        insertOrUpdate(geneKey, hashTable);

        HashNode *node = hashTable[hashFunction(geneKey)];

        while (node != NULL && strcmp(node->key, geneKey) != 0)
            node = node->next;

        if (node != NULL && (float)node->count / solutionPoolSize >= THRESHOLD)
        {
            flag[0] = 1;
            flag[1] = i;
            free(geneKey);
            return flag;
        }
    }

    free(geneKey);
    return flag;
}

// Proportionate roulette wheel selection

void selection(Chromosome **population, int numbits, int solutionPoolSize, float **W) // need to add float **W as an imput, here for test conveneice only, we use the pre generated weight W.
{
    Chromosome **new_population = initialize_population(numbits, solutionPoolSize);
    float total_fitness = 0.0;
    float cum_prob[512];
    memset(cum_prob, 0, solutionPoolSize * sizeof(float));
    float max_fitness = 0.0;

    for (int i = 0; i < solutionPoolSize; i++)
    {
        population[i]->fitness = evaluate(population[i], numbits, W); // need W or not?

        // fine the maximum fitness
        if (population[i]->fitness > max_fitness)
            max_fitness = population[i]->fitness;
    }

    for (int i = 0; i < solutionPoolSize; i++)
    {

        // Since we want to minimize our energy, solution with low energy should have high fitness
        // We inverse the fitness by subtracting current fitness from the maximum one, which also ensures all fitness greater or equal to 0
        // +1 to ensure strctly greater than 0 (those with low fitness also have chance to be chosen)
        population[i]->fitness = max_fitness - population[i]->fitness + 1;
        total_fitness += population[i]->fitness;
    }

    for (int i = 0; i < solutionPoolSize; i++)
    {
        // cumulative probability
        cum_prob[i] = (i == 0) ? (population[i]->fitness / total_fitness) : (cum_prob[i - 1] + population[i]->fitness / total_fitness);
    }

    for (int i = 0; i < solutionPoolSize; i++)
    {
        float rand_num = (float)rand() / RAND_MAX;
        for (int j = 0; j < solutionPoolSize; j++)
        {
            // low fitness with low probability to be chosen
            if (rand_num <= cum_prob[j])
            {
                memcpy(new_population[i]->genes, population[j]->genes, numbits * sizeof(int));
                new_population[i]->fitness = population[j]->fitness;
                break;
            }
        }
    }

    // copy new population back to current population
    for (int i = 0; i < solutionPoolSize; i++)
    {
        memcpy(population[i]->genes, new_population[i]->genes, numbits * sizeof(int));
        population[i]->fitness = new_population[i]->fitness;
    }

    free_population(new_population, solutionPoolSize);
}

void crossover(Chromosome **population, int numbits, int solutionPoolSize)
{
    for (int i = 0; i < solutionPoolSize - 1; i += 2)
    {
        float rand_num = (float)rand() / RAND_MAX;

        if (rand_num < CROSS_PROB)
        {
            int crossover_point = rand() % numbits;

            // even chromosomes with odd ones to crossover
            for (int j = crossover_point; j < numbits; j++)
            {
                int temp = population[i]->genes[j];
                population[i]->genes[j] = population[i + 1]->genes[j];
                population[i + 1]->genes[j] = temp;
            }
        }
    }
}

void mutation(Chromosome **population, int numbits, int solutionPoolSize)
{
    for (int i = 0; i < solutionPoolSize; i++)
        for (int j = 0; j < numbits; j++)
        {
            float rand_num = (float)rand() / RAND_MAX;
            if (rand_num < MUTATION_RATE)
                population[i]->genes[j] = 1 - population[i]->genes[j]; // 1 to 0 or 0 to 1
        }
}

// Call to run the whole genetic algorithm, input a population, output the same length population with mutation.

void runGeneticAlgor(Chromosome **population, int numbits, int solutionPoolSize, float **W)
{
    HashNode *hashTable[TABLE_SIZE];
    int *flag;

    for (int i = 0; i < TABLE_SIZE; i++)
        hashTable[i] = NULL;

    flag = count(population, hashTable, numbits, solutionPoolSize);

    if (flag[0] == 1)
        return;

    selection(population, numbits, solutionPoolSize, W);
    crossover(population, numbits, solutionPoolSize);
    mutation(population, numbits, solutionPoolSize);

    free(flag);
    freeHashTable(hashTable);

    printf("\n");
}