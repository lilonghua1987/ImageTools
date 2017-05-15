#pragma once

#include "ImageProcess.h"


struct AngularPoint
{
	int x;
	int y;
	double r;
	bool isAngular;
	bool isMax;
};

class Harris
{
private:
	int _IMG_HEIGHT_;
	int _IMG_WIDTH_;
public:
	Harris(void);
	Harris(Image img);
	~Harris(void);

	double**  difference(Image img);

	double* GausscianBlur2D(double *matrixData,int height,int width,int window,double sigma);

	double* GausscianSeparateBlur(double *matrixData,int height,int width,int window,double sigma);

	/*
	* k : �ǵ��������ֵ���ڴ�С
	* thresh : ��ֵ����
	*/
	AngularPoint** checkAngular(double**,int window,double thresh);

	void printAngular(AngularPoint** angular);

	vector<vector<AngularPoint>> getAngular(AngularPoint** angular);

	void printAngular(vector<vector<AngularPoint>>,Image img);

	ParallaxMatrix NCCByAngular(Image referenceImg,Image targetImg,int H,int W,int dispMax);
};

