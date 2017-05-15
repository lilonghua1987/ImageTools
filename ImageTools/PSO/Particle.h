#pragma once

#include <omp.h>
#include <vector>
#include <time.h>
#include "../Image.h"
#include "../tools.h"

typedef int nodeType;

class Swarm;

/// 粒子类
class Particle
{
	public:
		Swarm* swarm;

		//Original Position
		Point2D<nodeType> point;

		double cost;

		Point2D<nodeType> bestPosition;

		//------------------------------------------------------
		/// 构造函数;
		Particle();

		Particle(Swarm* swarm, const Point2D<nodeType>& point, const Point2D<nodeType>& bestPosition, double cost = 0.0);

		//------------------------------------------------------
		/// 析构函数;
		virtual ~Particle();

		//------------------------------------------------------
		///	消费函数，计算粒子的消费值;
		virtual void CalculateCost() = 0;

		virtual double dataCost() = 0;
		virtual double smoothCost() = 0;
};