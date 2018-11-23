#pragma once
#include "pch.h"

#ifndef INDIVIDUAL_H
#define INDIVIDUAL_H

using namespace std;

class Individual
{
public:
	using Comparator = bool(Individual::*)(Individual, Individual);

	vector<float> genes;
	vector<float> objectiveSet;
	vector<Individual*> dominatedSet;
	int dominatedCount;
	int rank;
	float distance;

	Individual();
	Individual(int);
	void evaluation();
	bool dominate(Individual);
};

#endif // !INDIVIDUAL_H
