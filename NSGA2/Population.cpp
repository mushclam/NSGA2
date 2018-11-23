#include "pch.h"
#include "Population.h"

Population::Population(int _popSize, int _geneSize, float probCrossover, float probMutation)
{
	populationSize = _popSize;
	geneSize = _geneSize;
	crossover_prob = probCrossover;
	mutation_prob = probMutation;
}

void Population::initialize()
{
	for (int i = 0; i < populationSize; i++)
	{
		individualSet.push_back(Individual(geneSize));
	}
}

void Population::clear()
{
	individualSet.clear();
}

Population Population::copy()
{
	return Population(populationSize, geneSize, crossover_prob, mutation_prob);
}

Population Population::copy_all()
{
	Population _tmpPop(populationSize, geneSize, crossover_prob, mutation_prob);
	_tmpPop.individualSet = individualSet;
	return _tmpPop;
}

Population Population::combination(Population q)
{
	Population _tmp(populationSize * 2, geneSize * 2, crossover_prob, 1 / (populationSize * 2));
	_tmp.individualSet = individualSet;

	for (auto ind : q.individualSet)
	{
		_tmp.individualSet.push_back(ind);
	}
	return _tmp;
}

void Population::evaluation()
{
	for (auto & ind : individualSet)
	{
		ind.evaluation();
	}
}