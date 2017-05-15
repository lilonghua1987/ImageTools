#pragma once
#ifndef MYPSO_H
#define MYPSO_H

#include "Swarm.h"

class CostFunction;

class MinParticle : public Particle
{
    private:
	    CostFunction* function;
		static const Point2D<int> NEIGHBORS[];
	public:
		//------------------------------------------------------
		/// 构造函数；
		MinParticle(CostFunction* function,Swarm* swarm, const Point2D<nodeType>& point);
		
		//------------------------------------------------------
		/// 析构函数；
		virtual ~MinParticle();

		//------------------------------------------------------
		void CalculateCost();

		double dataCost();
		double smoothCost();

protected:
	void init();

};



class MinSwarm : public Swarm
{
	public:
		//------------------------------------------------------
		MinSwarm(int disp, int swarmID);

		virtual ~MinSwarm(){};

		//------------------------------------------------------
		void AddParticle(CostFunction* function,const Point2D<nodeType>& point);
};

#endif