#include "MyPSO.h"
#include "CostFunction.h"

const Point2D<int> MinParticle::NEIGHBORS[] = {Point2D<int>(-1,0),Point2D<int>(1,0),Point2D<int>(0,-1),Point2D<int>(0,1)};/*{ Point2D<int>(1, 0), Point2D<int>(0, -1) };*/

MinParticle::MinParticle(CostFunction* function,Swarm* swarm, const Point2D<nodeType>& point)
{
	this->function = function;
	this->swarm = swarm;
	this->point = point;
	init();
	bestPosition = swarm->bestPosition;
}


MinParticle::~MinParticle()
{
}


void MinParticle::CalculateCost()
{
	//cost = dataCost() + smoothCost();
	cost = dataCost();
}


double MinParticle::dataCost()
{
	Point2D<nodeType> pd;
	double sum = 0.0;
	if (function->dispMap.at(point.y, point.x) == 0)
	{
		sum = function->truncate;
	}
	else
	{
		for (auto& p : NEIGHBORS)
		{
			pd = p + point;
			sum += function->dataCost(point, pd);
		}

		sum /= 4;
	}
	return sum;
}


double MinParticle::smoothCost()
{
	Point2D<nodeType> pd;
	double sum = 0.0;
	if (function->dispMap.at(point.y,point.x) == 0)
	{
		sum = function->truncate;
	} 
	else
	{
		for (auto& p:NEIGHBORS)
		{
			pd = p + point;
			sum += function->smoothCost(point,pd);
		}
	}
	return sum/4;
}


void MinParticle::init()
{
	CalculateCost();
}


