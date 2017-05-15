#pragma once
#include "MyPSO.h"
#include "../SLIC.h"
#include "CostFunction.h"

class PSOStereoMatch : public CostFunction
{
public:
	PSOStereoMatch(const Image& imL, const Image& imR, int disp = 50);
	virtual ~PSOStereoMatch(void);

public:
	virtual double dataCost(const Point2D<nodeType>& pL, const Point2D<nodeType>& pR)
	{
		double d = 0.0,sum = 0.0;

		for (ulong i = 0; i < labL.channel; i++)
		{
			d = labL.at(pL.y,pL.x,i) - labR.at(pR.y,pR.x,i);
			sum += abs(d);
		}

		return tools::Min(truncate,sum);
	};

	virtual double smoothCost(const Point2D<nodeType>& p, const Point2D<nodeType>& pd)
	{
		if (pd < Point2D<nodeType>(imL.width,imL.height) && pd >= Point2D<nodeType>(0,0))
		{
			nodeType R = 1;

			nodeType d1 = dispMap.at(p.y,p.x);
			nodeType d2 = dispMap.at(pd.y,pd.x);
			unsigned long s1 = p.y * imL.width + p.x;
			unsigned long s2 = pd.y * imL.width + pd.x;

			if (d1 == d2) return 0;
			else
			{
				R = abs(d1 - d2);

				if(labelsL.at(s1) == labelsL.at(s2))
				{
					if(R > 1) R = 5;
				}else
				{
					/*R = 1;*/
					return 0;
				}
			}

			float d = 0.0f, sum = 0.0f;
			for (ulong i = 0; i < labL.channel; i++)
			{
				d = labL.at(p.y,p.x,i) - labL.at(pd.y,pd.x,i);
				sum += abs(d);
			}
			//return R*sum/labL.channel;

			if (sum < truncate/3)
				return R*10/labL.channel;
			else                          
				return R*20/labL.channel;
		}else
		{
			return 0;
		}
	};

	virtual void RunMatch();

private:
	Matrix<float> labL,labR;
	std::vector<int> labelsL,labelsR;
};

