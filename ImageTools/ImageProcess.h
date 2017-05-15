#pragma once

#include "ImageExpandTools.h"
#include "ColorConversion.h"
#include "tools.h"
#include <vector>
#include <algorithm>
#include<cmath>
#include<iostream>
#include <string>
#include<complex>

#ifndef PI
#define PI 3.14159265
#endif

#ifndef PATH
#define PATH "temp\\"
#endif

using namespace std;

struct GreyHistogram
{
	vector<int> colorValue;

	vector<int> colorTimes;

	vector<double> colorRate;

	long unsigned pixelNU;

};

struct MatrixMat
{
	vector<vector<double>> data;

	int width;

	int height;

	double operator *(MatrixMat matrix){
		double value=0;

		for(int h=0;h<height;h++){
			for(int w=0;w<width;w++){
				value+=data[h][w]*matrix.data[h][w];
			}
		}

		return value;
	}

};


struct FeaturePoint
{
	double pixel;
	int indexO;
	int indexS;
	int indexX;
	int indexY;
	int directX;
	int directY;
	bool min;

};


struct ParallaxPoint
{
	int parallax;

	double relation;
};


struct ParallaxMatrix
{
	vector<vector<ParallaxPoint>> value;

	int width;

	int height;
};


class ImageProcess
{

private:
	//int wp,hp;//用于快速傅立叶中的蝶形算法

public:
	ImageProcess(void);

	virtual ~ImageProcess(void);

	static Image imageReverse(const Image& img);

	static GreyHistogram getGreyHistogram(const Image& img);

	static void printHistogram(GreyHistogram histogram);

	static Image HistogramEqualization(const Image& img);

	static Image GausscianBlur(const Image& img,double sigma);

	static Image GausscianSeparateBlur(const Image& img,double sigma);

	static Matrix<float> GausscianSeparateBlur(const Matrix<float>& img,double sigma);

	static Image halfSize(const Image& img);

	static Image doubleSize(const Image& img);

	MatrixMat inverseMatrixThree(const MatrixMat& matrix);

	static FeaturePoint PixelInterpolation(FeaturePoint,Image img,vector<double>);

	static Image NCC(const Image& imgL,const Image& imgR,int W,int dispMax);

	static void printMapingPointMatrix(ParallaxMatrix,FREE_IMAGE_TYPE img_type,FREE_IMAGE_FORMAT fif);

	static void fft(const complex<double> *a, complex<double> *y, int power);

	static void ifft(const complex<double> y[], complex<double> a[], int power);

	static Image adjustImageSize(Image image,ulong &wp,ulong &hp);

	static complex<double> * fourier(Image image);

	static void inverseFourier(complex<double> *,int wp,int hp,FREE_IMAGE_TYPE img_type);

	static void printFourierData(complex<double> *,int wp,int hp);

	static void dft(const complex<double>*,complex<double>*,int);

	static void idft(const complex<double>*,complex<double>*,int);

	static complex<double> * dFourier(Image image);

	static void inverseDFourier(complex<double> *,ulong width,ulong height,FREE_IMAGE_TYPE img_type);

	static void printDFourierData(complex<double> *fd,ulong width,ulong height);

	static void parseDFourier(Image referenceImg,Image targetImg);

	static void parseFourier(Image referenceImg,Image targetImg);

	static Image BilinearInterpolation(const Image& sourceImg,ulong height,ulong width);

	static Image BilateralFilter(const Image& img,int w,double sigmaC,double sigmaS);
	static Matrix<float> BilateralFilterWithLab(const Image& img,int w,double sigmaC,double sigmaS);
	static void RecursiveBilateralFilter(const Image& img,Matrix<double>& imgFilter,double sigmaS,double sigmaR);

	static Image MeanShift(const Image& img, const float Spatial = 50.0f, const float Color = 50.0f, const unsigned int itre = 30);
};

