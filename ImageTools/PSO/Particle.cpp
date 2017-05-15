#include "Particle.h"
#include <cmath>
#include "Swarm.h"
#include <limits>

Particle::Particle():cost(0.0),swarm(nullptr){}

Particle::Particle(Swarm* swarm, const Point2D<nodeType>& point, const Point2D<nodeType>& bestPosition, double cost)
: swarm(swarm)
, point(point)
, bestPosition(bestPosition)
, cost(cost)
{

}


Particle::~Particle()
{
	//std::cout << "V.x = " << velocity.x << " V.y" << velocity.y << std::endl;
}