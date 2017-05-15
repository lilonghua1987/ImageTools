#pragma once

#include "ImageExpandTools.h"
#include "ImageProcess.h"
#include <vector>

using namespace std;

struct GaussianPyramid
{
	vector<vector<Image*> > img;
	vector<double> factor;
};

struct DOGPyramid
{
	vector<vector<Image*> > img;
};


struct ExtremePointMatrix
{
	vector<vector<FeaturePoint>> layer;
};

struct ExtremeMatrix
{
	vector<ExtremePointMatrix> group;
};

class Sift
{
public:
	Sift(void);
	Sift(Image img);
	Sift(Image img, double sigma, int k, int octave, int sub_level, int S);
	~Sift(void);
	void init(double sigma, int k, int octave, int sub_level, int S);

	GaussianPyramid createGaussianPyramid(Image img);

	void printGaussianPyramid(GaussianPyramid pyramid);

	DOGPyramid createDOGPyramid(GaussianPyramid pyramid);

	void printDOGPyramid(DOGPyramid dog);

	ExtremeMatrix createExtremeMatrix(DOGPyramid dog);

	void printExtremeMatrix(ExtremeMatrix eMatrix);

private:
	Image img;
	double sigma;
	int k;
	int octave;
	int sub_level;
	int S;
};

