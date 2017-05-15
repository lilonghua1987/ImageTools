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
	/// ���캯����
	Swarm();

	//------------------------------------------------------
	/// ����������
	virtual ~Swarm();

	//------------------------------------------------------
	int SwarmSize();

	//------------------------------------------------------
	double CurrentBestCost();

	//------------------------------------------------------
	virtual double Iteration();

	void SortParticles();

};