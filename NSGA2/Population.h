#include "pch.h"
#include "Individual.h"

#pragma once

#ifndef POPULATION_H
#define POPULATION_H

using namespace std;

class Population
{
public:
	vector<Individual> individualSet;
	int populationSize;
	int geneSize;
	float crossover_prob;
	float mutation_prob;

	Population(int, int, float, float);
	void initialize();
	void clear();
	Population copy();
	Population copy_all();
	Population combination(Population);

	void evaluation();
};

#endif // !POPULATION_H