#pragma once

#include "Particle.h"

class Swarm
{
public:
	//SwarmID
	int swarmID;
	//------------------------------------------------------
	std::vector<Particle*> particles;

	//------------------------------------------------------
	Point2D<nodeType> bestPosition;

	//------------------------------------------------------
	double bestCost;

	int maxPostion;

	//
	Particle* leftP;
	Particle* rightP;
	Particle* topP;
	Particle* bottomP;

	//------------------------------------------------------
	/// 构造函数；
	Swarm();

	//------------------------------------------------------
	/// 析构函数；
	virtual ~Swarm();

	//------------------------------------------------------
	int SwarmSize();

	//------------------------------------------------------
	double CurrentBestCost();

	//------------------------------------------------------
	virtual double Iteration();

	void SortParticles();

};