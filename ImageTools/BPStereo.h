//===========================================================================
// For Example:
//BPStereo bp;
//	/*Image result = bp.bpMatch(ImageExpandTools::ImageToGrey(imgLeft.c_str()),ImageExpandTools::ImageToGrey(imgRight.c_str()));
//	ImageTools::SaveImage(result,"temp//bp.jpg");*/
//	/*Image result = bp.restore(ImageExpandTools::ImageToGrey(imgLeft.c_str()));
//	ImageTools::SaveImage(result,"temp//restore.jpg");*/
//===========================================================================

#pragma once

#include "ImageExpandTools.h"
#include "ImageProcess.h"

class BPStereo
{
public:
	BPStereo(int disp,int iter = 5,int levels = 5,float sigma = 0.7);
	~BPStereo(void);

	Image bpMatch(const Image& imL,const Image& imR);
	Image restore(const Image& img);
private:
	static const float discK; // truncation of discontinuity cost
	static const float dataK; // truncation of data cost
	static const float lambda; // weighting of data cost
	static const int infVal = 1000000; // large cost

	int disp; // number of possible disparities
	int iter; // number of BP iterations at each scale
	int levels; // number of scales
	float sigma; // amount to smooth the input images

private:
	// dt of 1d function
    void dt(float f[]);

	// compute message
    void msg(float s1[], float s2[], float s3[], float s4[], float dst[]);
	// computation of data costs
	Matrix<float> dataCosts(const Image& imL,const Image& imR);
	// computation of data costs
	Matrix<float> dataCosts(const Image& img);
	// belief propagation using checkerboard update scheme
	void bp(Matrix<float>& u, Matrix<float>& d,Matrix<float>& l, Matrix<float>& r,Matrix<float>& data,int iter);
	// generate output from current messages
	Image generate(Matrix<float>& u, Matrix<float>& d, Matrix<float>& l, Matrix<float>& r, Matrix<float>& data);
};

