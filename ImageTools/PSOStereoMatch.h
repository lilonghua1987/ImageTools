#pragma once
#include "MyPSO.h"
#include "tools.h"
#include "../SLIC.h"

class PSOStereoMatch
{
public:
	PSOStereoMatch(const Image& imL, const Image& imR, int disp = 50);
	~PSOStereoMatch(void);


	double Cost(const Point2D<double>& pL,const Point2D<double>& pR);

	friend double CostFun(const Point2D<double>& pL, const vector<double>& pRL)
	{
		double cost = 100000.0;
		for (auto p:pRL)
		{
			if ((pL.x + p) < imL.width || (pL.x + p) >= 0)
			{
				double tc = Cost(pL,Point2D<double>(pL.x + p,pL.y));
				if(tc < cost) cost = tc;
			}
		}
		return cost;
	}

	void RunMatch();

private:
	Image imL,imR;
	int disp;
};

