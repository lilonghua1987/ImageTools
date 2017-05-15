#include "Swarm.h"
#include <algorithm>
#include <functional>

Swarm::Swarm()
{
	bestCost = std::numeric_limits<double>::max();

	leftP = nullptr;
	rightP = nullptr;
	topP = nullptr;
	bottomP = nullptr;
}


Swarm::~Swarm()
{
	for (auto particle:particles)
	{
		delete particle;
	}
}


int Swarm::SwarmSize()
{
	return particles.size();
}


double Swarm::CurrentBestCost()
{
	if(particles.empty()) return bestCost;
	double sumCost = 0.0;
	for(auto& particle:particles)
	{
		sumCost += particle->cost;
	}
	return sumCost/particles.size();
}


double Swarm::Iteration()
{
	for (auto& particle:particles)
	{
		particle->CalculateCost();
	}

	double currentCost = CurrentBestCost();

	return bestCost;
}


void Swarm::SortParticles()
{
	sort(particles.begin(), particles.end(), [](const Particle* p1, const Particle* p2) -> bool {return p1->cost < p2->cost; });
}
