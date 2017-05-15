#include "MyPSO.h"

MinSwarm::MinSwarm(int disp, int swarmID)
{
	this->maxPostion = disp;
	this->swarmID = swarmID;
}



void MinSwarm::AddParticle(CostFunction* function,const Point2D<nodeType>& point)
{
	particles.push_back(new MinParticle(function,this,point));

	if (leftP == nullptr) leftP = particles.at(particles.size() - 1);
	else
	{
		if (leftP->point.x > point.x) leftP = particles.at(particles.size() - 1);
	}

	if (rightP == nullptr) rightP = particles.at(particles.size() - 1);
	else
	{
		if (rightP->point.x < point.x) rightP = particles.at(particles.size() - 1);
	}

	if (topP == nullptr) topP = particles.at(particles.size() - 1);
	else
	{
		if (topP->point.y > point.y) topP = particles.at(particles.size() - 1);
	}

	if (bottomP == nullptr) bottomP = particles.at(particles.size() - 1);
	else
	{
		if (bottomP->point.y < point.y) bottomP = particles.at(particles.size() - 1);
	}
}