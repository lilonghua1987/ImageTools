#pragma once
#include "ImageProcess.h"

class sharpen
{
public:
	sharpen(void);
	~sharpen(void);

	double*  laplacianGrad(const Image &img);

	void printLaplacianGrad(double *grad,const Image &templateImg);

	Image  laplacian(const Image &img,int k);

	void CreatGauss(double sigma, double **pdKernel, int *pnWidowSize);

	Image GaussianSmooth(const Image &img,double sigma);

	void Grad(const Image &img,int *pGradX, int *pGradY, int *pMag);

	void NonmaxSuppress(int *pMag, int *pGradX, int *pGradY, Image &img);

	void EstimateThreshold(const Image &img,int *pMag,int *pThrHigh, int *pThrLow,double dRatHigh, double dRatLow);

	void Hysteresis(Image &img, int *pMag,double dRatLow, double dRatHigh);

	void TraceEdge(Image &img, int y, int x, int nThrLow,int *pMag);

	void Canny(const Image &img, double sigma, double dRatLow, double dRatHigh);
};

