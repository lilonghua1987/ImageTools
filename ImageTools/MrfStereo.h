#pragma once

#include "mrf/mrf.h"
#include "mrf/ICM.h"
#include "mrf/GCoptimization.h"
#include "mrf/MaxProdBP.h"
#include "mrf/TRW-S.h"
#include "mrf/BP-S.h"
#include "Image.h"
#include "ImageTools.h"
#include "Matrix.h"
#include "tools.h"
#include "ColorConversion.h"

namespace mrf
{
	class Parameters
	{

	public:
		Parameters(int nD = 16, int birchfield = 1, int squaredDiffs = 1, 
			       int truncDiffs = 255, int MRFalg = 1, int smoothexp = 1, 
				   int smoothmax = 2, int lambda = 20, int gradThresh = -1,
				   int gradPenalty = 2)
				   :nD(nD)
				   ,birchfield(birchfield)
				   ,squaredDiffs(squaredDiffs)
				   ,truncDiffs(truncDiffs)
				   ,MRFalg(MRFalg)
				   ,smoothexp(smoothexp)
				   ,smoothmax(smoothmax)
				   ,lambda(lambda)
				   ,gradThresh(gradThresh)
				   ,gradPenalty(gradPenalty)
		{
		}

		~Parameters(){}

	public:
		// parameters controlled via command-line options:
		int nD;           // disparity levels (d = 0 .. nD-1)
		int birchfield;    // use Birchfield/Tomasi costs
		int squaredDiffs;  // use squared differences (absolute differences by default)
		int truncDiffs;  // truncated differences (before squaring), by default not
		int MRFalg;        // 0-ICM, 1-GC/expansion (default), 2-GC/swap, 3-TRWS, 4-BPS, 5-BPM, 9-all
		int smoothexp;     // exponent of smoothness term: 1 (default) or 2, i.e. L1 or L2 norm
		int smoothmax;     // maximum value of smoothness term (2 by default)
		int lambda;       // weight of smoothness term (20 by default)
		int gradThresh;   // intensity gradient cue threshold, by default none
		int gradPenalty;   // if grad < gradThresh, multiply smoothness cost by this
	};
}

enum algtypes {aICM, aExpansion, aSwap, aTRWS, aBPS, aBPM};

class MrfStereo
{
public:
	MrfStereo(const Image& imgL,const Image& imgR,const mrf::Parameters& paramet);
	~MrfStereo(void);

	void computeDSI(const Image& imL, const Image& imR, MRF::CostVal *&dsi, int nD, int birchfield, int squaredDiffs, int truncDiffs);
	void computeCues(const Image& img, MRF::CostVal *&hCue, MRF::CostVal *&vCue,int gradThresh, int gradPenalty);
	void WTA(MRF::CostVal *dsi, int nD, Image &disp);
	void getDisparities(MRF *mrf, Image &disp);
	void setDisparities(const Image& disp, MRF *mrf);

	void dataCost();

	friend void MrfStereoMatch(const Image& imL, const Image& imR);

private:
	Image imgL;   // reference frame
	Image imgR;  // matching image frame
	Matrix<float> labL,labR;
	Matrix<float> costRaw;

	mrf::Parameters paramet;
	MRF *mrf;
	DataCost *dcost;
	SmoothnessCost *scost;

	const int size;
	const double truncate;

private:
	double weight(const Matrix<float>& lab, const Point2D<long>& pL,const Point2D<long>& pR)
	{
		double rc = 5.0;
		double rp = (2 * size + 1)/2.0;

		double d = 0.0,dC = 0.0;
		for (ulong i = 0; i < lab.channel; i++)
		{
			d = lab.at(pL.y,pL.x,i) - lab.at(pR.y,pR.x,i);
			dC += d*d;
		}
		dC = sqrt(dC);
		double dP = tools::EuclideanDistance(pL,pR);
		dC/=rc;
		dP/=rp;
		return exp(-(dC+dP));
	};

	double rawCost(const Point2D<long>& pL,const Point2D<long>& pR)
	{
		double d = 0.0,sum = 0.0;

		for (ulong i = 0; i < labL.channel; i++)
		{
			d = labL.at(pL.y,pL.x,i) - labR.at(pR.y,pR.x,i);
			sum += abs(d);
		}

		return tools::Min(truncate,sum);
		//return sum;
	};
};
