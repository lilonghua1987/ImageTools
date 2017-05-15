#pragma once

#include <vector>
#include "../Image.h"
#include "../Matrix.h"

class CostFunction
{
public:
	CostFunction(const Image& imL, const Image& imR, int disp = 50)
		:imL(imL)
		,imR(imR)
		,disp(disp)
		,truncate(40.0)
	{
	};
	virtual ~CostFunction(void){};

public:
	virtual double dataCost(const Point2D<nodeType>& pL, const Point2D<nodeType>& pRL) = 0;
	virtual double smoothCost(const Point2D<nodeType>& p, const Point2D<nodeType>& pd) = 0;

	virtual void RunMatch() = 0;

protected:
	Image imL,imR;
	const int disp;
public:
	const double truncate;
	Matrix<nodeType> dispMap;
};

