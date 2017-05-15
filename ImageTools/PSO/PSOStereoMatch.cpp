#include "PSOStereoMatch.h"


PSOStereoMatch::PSOStereoMatch(const Image& imL, const Image& imR, int disp)
	:CostFunction(imL,imR,disp)
{
	dispMap = Matrix<nodeType>(imL.height,imL.width);
}


void PSOStereoMatch::RunMatch()
{
	int numlabelsL = 0, numlabelsR = 0;
	int m_spcount = 2200,m_compactness = 10;
	SLIC slicL;
	slicL.PerformSLICO_ForGivenK(imL, labelsL, numlabelsL, m_spcount, m_compactness);
	slicL.DrawContoursAroundSegments(imL, labelsL);
	/*slic.DrawContoursAroundSegmentsTwoColors(imL, labels);*/
	imL.save("temp/psoL.png");
	labL = slicL.lab;
	slicL.~SLIC();

	SLIC slicR;
	slicR.PerformSLICO_ForGivenK(imR, labelsR, numlabelsR, m_spcount, m_compactness);
	slicR.DrawContoursAroundSegments(imR, labelsR);
	labR = slicR.lab;
	imR.save("temp/psoR.png");
	slicR.~SLIC();

	std::vector<MinSwarm*> swarmsL(numlabelsL);

	for (size_t i = 0; i < labelsL.size(); i++)
	{
		if (swarmsL.at(labelsL.at(i)) == NULL) swarmsL.at(labelsL.at(i)) = new MinSwarm(disp, labelsL.at(i));
		swarmsL.at(labelsL.at(i))->AddParticle(this, Point2D<double>(i%imL.width, i / imL.width));
	}

	std::vector<MinSwarm*> swarmsR(numlabelsR);

	for (size_t i = 0; i < labelsR.size(); i++)
	{
		if (swarmsR.at(labelsR.at(i)) == NULL) swarmsR.at(labelsR.at(i)) = new MinSwarm(disp, labelsR.at(i));
		swarmsR.at(labelsR.at(i))->AddParticle(this, Point2D<double>(i % imL.width, i / imL.width));
	}

	std::vector<std::pair<MinSwarm*, MinSwarm*>> match;

	for (auto& swL:swarmsL)
	{
		double eL = 0.0;

		for (auto& p : swL->particles)
		{
			eL += p->cost;
		}

		eL /= swL->particles.size();

		eL += swL->leftP->cost;
		eL += swL->rightP->cost;
		eL += swL->topP->cost;
		eL += swL->bottomP->cost;

		double deltaE = std::numeric_limits<double>::max();

		MinSwarm* s = nullptr;

		for (auto& swR:swarmsR)
		{
			auto m = std::find_if(match.begin(), match.end(), [swR](const std::pair<MinSwarm*, MinSwarm*>& p){return (p.second == swR); });
			if (m != match.end()) continue;

			if ((swL->rightP->point.x - swR->rightP->point.x) <= abs(disp) && (swL->leftP->point.x >= swR->leftP->point.x) 
				&& (swL->topP->point.y < swR->bottomP->point.y) && (swL->bottomP->point.y > swR->topP->point.y))
			{
				double eR = 0.0;

				for (auto& p:swR->particles)
				{
					eR += p->cost;
				}

				eR /= swR->particles.size();

				eR += swR->leftP->cost;
				eR += swR->rightP->cost;
				eR += swR->topP->cost;
				eR += swR->bottomP->cost;

				double delta = eL - eR;

				if (deltaE > delta)
				{
					s = swR;
					deltaE = delta;
				}
			}
		}

		if (s) match.push_back(std::pair<MinSwarm*, MinSwarm*>(swL,s));
	}

	/*int count = 0,iter = 0,iterMax = 77;

	double sumEL = std::numeric_limits<double>::max();
	while (count < iterMax)
	{
		double sumE = 0.0;

        #pragma omp parallel for reduction(+:sumE)
		for (int i = 0; i < (int)swarmsL.size(); i++)
		{		
			sumE += swarmsL.at(i)->Iteration();
		}

		if(sumEL == sumE)
		{
			count++;
		}
		else
		{
			if(sumE < sumEL) sumEL = sumE;
			count = 0;
		}  
	}*/

	Image dispImg(dispMap.column, dispMap.row);

	for (auto& m:match)
	{
		int d = m.first->rightP->point.x - m.second->leftP->point.x;
		for (auto& p:m.first->particles)
		{
			dispImg.at<uchar>(p->point.y, p->point.x) = d * 4;
		}
	}

	

	/*for (ulong i = 0; i < dispMap.row; i++)
	{
		for (ulong j = 0; j < dispMap.column; j++)
		{
			dispImg.at<BYTE>(i, j) = (BYTE)dispMap.get(i,j);
		}
	}*/

	dispImg.save("temp/pso_venus.png");

	for (auto swarm : swarmsL)
	{
		delete swarm;
	}
}


PSOStereoMatch::~PSOStereoMatch(void)
{
}
