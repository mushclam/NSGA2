#pragma once
#include "pch.h"
#include "Population.h"
#include "Individual.h"

#ifndef GENETIC_ALGORITHM_H
#define GENETIC_ALGORITHM_H

typedef vector<Individual*> Front;

using namespace std;

class GeneticAlgorithm
{
public:
	void wheelSelection(Population&);
	void crossover(Population&);
	void mutation(Population&);
	void SBX(int, Population&);
	void PLM(int, float, float, Population&);
	vector<Front> fastNonDominatedSort(Population*);
	void crowdingDistanceAssignment(Front, float, float);

	void normalize();
	void associate();
	void niching();
};

#endif // !GENETIC_ALGORITHM_H
