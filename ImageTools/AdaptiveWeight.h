#pragma once

#include "ImageExpandTools.h"
#include "ColorConversion.h"
#include <vector>
#include <assert.h>

using namespace std;


class AdaptiveWeight
{
public:
	AdaptiveWeight(Image imgL,Image imgR);
	~AdaptiveWeight(void);

private:

	int T;
	int size;

	float rc;
	float rp;

	int dMin;
	int dMax;

public:
	Image imgL,imgR;
	/*vector<vector<Pixel<float> > > labImgL,labImgR;*/
	Matrix<float> greCost;

private:
	inline float getWeight(Matrix<Pixel<float> > labImg,int i,int j,int m,int n);
	inline float getTAD(const Image &imgL,const Image &imgR,int i,int j,int d);

public:
	void RawCost(const Image &imgL, const Image &imgR, Matrix<float> &rawCost, int range);

	void AdaptiveSupportWeightComputation(const Image &img,Matrix<float> &pSW,int mask);

	void SupportAggregation(Matrix<float> pSW_ref,Matrix<float> pSW_tar,const Matrix<float> pRawCost,int range,int mask);

	/*int computeByThread(int bX,int bY,int eX,int eY);*/

	void optimization();
};

