#include "pch.h"
#include "GeneticAlgorithm.h"

void GeneticAlgorithm::wheelSelection(Population &p)
{

}

void GeneticAlgorithm::crossover(Population &p)
{
	for (int i = 0; i < p.populationSize / 2; i++)
	{
		std::random_device rd;
		std::mt19937 rng(rd());
		std::uniform_real_distribution<> dist_real(0.0, 1.0);
		std::uniform_int_distribution<std::mt19937::result_type> dist_int(0, p.geneSize - 2);

		if (dist_real(rng) <= p.crossover_prob)
		{
			int pos[3];
			for (int j = 0; j < 3; j++)
			{
				pos[j] = dist_int(rng);
			}
			std::sort(pos, pos + 3);

			for (int j = pos[0] + 1; j < pos[1]; j++)
			{
				int aux = p.individualSet[i].genes[j];
				p.individualSet[i].genes[j] = p.individualSet[i + p.populationSize / 2].genes[j];
				p.individualSet[i + p.populationSize / 2].genes[j] = aux;
			}
			for (int j = pos[2] + 1; j < p.geneSize; j++)
			{
				int aux = p.individualSet[i].genes[j];
				p.individualSet[i].genes[j] = p.individualSet[i + p.populationSize / 2].genes[j];
				p.individualSet[i + p.populationSize / 2].genes[j] = aux;
			}
		}
	}
}

void GeneticAlgorithm::mutation(Population &p)
{
	for (int i = 0; i < p.populationSize; i++)
	{
		for (int j = 0; j < p.geneSize; j++)
		{
			std::random_device rd;
			std::mt19937 rng(rd());
			std::uniform_real_distribution<> dist(0, 1);

			if (dist(rng) < p.mutation_prob) {
				if (p.individualSet[i].genes[j] == 0) {
					p.individualSet[i].genes[j] = 1;
				}
				else {
					p.individualSet[i].genes[j] = 0;
				}
			}
		}
	}
}

void GeneticAlgorithm::SBX(int index, Population &p)
{
	vector<Individual> _tmpSet;
	std::random_device rd;
	std::mt19937 rng(rd());
	std::uniform_real_distribution<> dist_real(0.0, 1.0);

	for (int i = 0; i < p.populationSize / 2; i++)
	{
		Individual p1 = p.individualSet[i];
		Individual p2 = p.individualSet[i + (p.populationSize / 2)];

		Individual c1(p.geneSize);
		Individual c2(p.geneSize);

		for (int j = 0; j < p.geneSize; j++)
		{
			if (dist_real(rng) <= p.crossover_prob)
			{
				float u = dist_real(rng);
				float x_bar = (p1.genes[j] + p2.genes[j]) / 2;
				float b_bar;
				if (u <= 0.5)
				{
					b_bar = pow(2.0 * u, 1.0 / ((float)index + 1.0));
				}
				else
				{
					b_bar = pow(1.0 / 2.0 * (1.0 - u), 1.0 / ((float)index + 1.0));
				}
				c1.genes[j] = x_bar - b_bar * abs(p2.genes[j] - p1.genes[j]) / 2;
				c2.genes[j] = x_bar + b_bar * abs(p2.genes[j] - p1.genes[j]) / 2;
			}
		}
		_tmpSet.push_back(c1);
		_tmpSet.push_back(c2);
	}
	p.individualSet = _tmpSet;
}

void GeneticAlgorithm::PLM(int index, float lower, float upper, Population &p)
{
	std::random_device rd;
	std::mt19937 rng(rd());
	std::uniform_real_distribution<> dist_real(0, 1);

	for (int i = 0; i < p.individualSet.size(); i++)
	{
		Individual & ind = p.individualSet[i];
		for (int j = 0; j < p.geneSize; j++)
		{
			if (dist_real(rng) <= p.mutation_prob)
			{
				float u = dist_real(rng);
				if (u <= 0.5)
				{
					float delta_l = pow(2.0 * u, 1.0 / (1.0 + (float)index)) - 1.0;
					ind.genes[j] = ind.genes[j] + delta_l * (ind.genes[j] - lower);
				}
				else
				{
					float delta_r = 1.0 - pow(2.0 * (1.0 - u), 1.0 / (1.0 + (float)index));
					ind.genes[j] = ind.genes[j] + delta_r * (upper - ind.genes[j]);
				}
			}
		}
	}
}

vector<Front> GeneticAlgorithm::fastNonDominatedSort(Population * r)
{
	vector<Front> front;
	front.resize(1);

	for (int i = 0; i < r->populationSize; i++)
	{
		Individual * p = &r->individualSet[i];
		p->dominatedSet.clear();
		p->dominatedCount = 0;

		for (int j = 0; j < r->populationSize; j++)
		{
			if (i == j) { continue; }
			Individual * q = &r->individualSet[j];
			if (p->dominate(*q))
			{
				p->dominatedSet.push_back(q);
			}
			else if (q->dominate(*p))
			{
				p->dominatedCount++;
			}
		}

		if (p->dominatedCount == 0)
		{
			p->rank = 0;
			front[0].push_back(p);
		}
	}
	
	int i = 0;
	while (front[i].size() != 0)
	{
		vector<Individual*> _tmpSet;
		_tmpSet.resize(0);
		for (int j = 0; j < front[i].size(); j++)
		{
			Individual * p = front[i][j];
			for (int k = 0; k < p->dominatedSet.size(); k++)
			{
				Individual * q = p->dominatedSet[k];
				q->dominatedCount--;
				if (q->dominatedCount == 0)
				{
					q->rank = i + 1;
					_tmpSet.push_back(q);
				}
			}
		}
		i++;
		front.push_back(_tmpSet);
	}

	return front;
}

bool comparator_0(Individual * a, Individual * b)
{
	return a->objectiveSet[0] < b->objectiveSet[0];
}

bool comparator_1(Individual * a, Individual * b)
{
	return a->objectiveSet[1] < b->objectiveSet[1];
}

void GeneticAlgorithm::crowdingDistanceAssignment(Front f, float max, float min)
{
	int geneSize = 1;
	int objectiveSize = 2;
	Individual f_max(geneSize);
	Individual f_min(geneSize);
	f_max.genes[0] = max;
	f_min.genes[0] = min;
	f_max.evaluation();
	f_min.evaluation();

	int l = f.size();
	for (int i = 0; i < l; i++)
	{
		f[i]->distance = 0;
	}
	for (int m = 0; m < objectiveSize; m++)
	{
		if (m == 0)
		{
			sort(f.begin(), f.end(), comparator_0);
		}
		else if (m == 1)
		{
			sort(f.begin(), f.end(), comparator_1);
		}
		f[0]->distance = numeric_limits<float>::infinity();
		f[l-1]->distance = numeric_limits<float>::infinity();
		for (int j = 2; j < l - 1; j++)
		{
			f[j]->distance = f[j]->distance + (f[j + 1]->objectiveSet[m] - f[j - 1]->objectiveSet[m]) / (f_max.objectiveSet[m] - f_min.objectiveSet[m]);
		}
	}
}