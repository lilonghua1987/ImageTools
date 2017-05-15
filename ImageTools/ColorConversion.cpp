#include "ColorConversion.h"


const double ColorConversion::epsilon = 0.008856;
const double ColorConversion::kappa = 903.3;

ColorConversion::ColorConversion(void)
{
}


ColorConversion::~ColorConversion(void)
{
}


// converting rgb values into CIELab values in a whole image
Matrix<Pixel<float> > ColorConversion::ColorSpaceConversionFromRGB2Lab(const Image &pRGB)
{

	assert(pRGB.channel == 3);

	Matrix<Pixel<float> > labImg(pRGB.height,pRGB.width);

	for(ulong i=0; i<pRGB.height; i++)
	{
		for(ulong j=0; j<pRGB.width; j++)
		{
			Pixel<int> point = pRGB.get<BYTE>(i,j);
			Pixel<double> lab;
			RGBtoLab(point.red, point.green, point.blue, lab.red, lab.green, lab.blue);		
			
			labImg.at(i,j) = lab;
		}
	}

	return labImg;
}



Matrix<float> ColorConversion::RGBToLab(const Image& img)
{

	assert(img.channel == 3);

	Matrix<float> labImg(img.height,img.width,img.channel);

	for(ulong i=0; i<img.height; i++)
	{
		for(ulong j=0; j<img.width; j++)
		{
			Pixel<double> point = img.get<double>(i,j);
			Pixel<double> lab;
			RGBtoLab(point.red, point.green, point.blue, lab.red, lab.green, lab.blue);	

			if(labImg.channel > 0) labImg.at(i,j,0) = (float)lab.red;
			if(labImg.channel > 1) labImg.at(i,j,1) = (float)lab.green;
			if(labImg.channel > 2) labImg.at(i,j,2) = (float)lab.blue;
			if(labImg.channel > 3) labImg.at(i,j,3) = (float)lab.alpha;
		}
	}

	return labImg;
}


void  ColorConversion::RGBToLab(const Image& img, Matrix<float>& dest)
{
	assert(img.data && dest.data);
	assert(img.channel == 3 && dest.channel == img.channel);
	assert(img.height > 0 && dest.row == img.height);
	assert(img.width > 0 && dest.column == img.width);

	for(ulong i=0; i<img.height; i++)
	{
		for(ulong j=0; j<img.width; j++)
		{
			Pixel<double> point = img.get<double>(i,j);
			Pixel<double> lab;
			RGBtoLab(point.red, point.green, point.blue, lab.red, lab.green, lab.blue);	

			if(dest.channel > 0) dest.at(i,j,0) = (float)lab.red;
			if(dest.channel > 1) dest.at(i,j,1) = (float)lab.green;
			if(dest.channel > 2) dest.at(i,j,2) = (float)lab.blue;
			if(dest.channel > 3) dest.at(i,j,3) = (float)lab.alpha;
		}
	}
}


void ColorConversion::LabtoRGB(const Matrix<float> &imgLab,Image &dest)
{

	for(ulong i=0; i<imgLab.row; i++)
	{
		for(ulong j=0; j<imgLab.column; j++)
		{
			Pixel<double> lab;

			if(imgLab.channel > 0) lab.red = imgLab.at(i,j,0);
			if(imgLab.channel > 1) lab.green = imgLab.at(i,j,1);
			if(imgLab.channel > 2) lab.blue = imgLab.at(i,j,2);

			Pixel<double> point;
			LabtoRGB(lab.red, lab.green, lab.blue,point.red,point.green,point.blue);
			dest.set(i,j,point);
		}
	}
}


Matrix<float> ColorConversion::RGBToLUV(const Image& img)
{
	assert(img.channel == 3);

	Matrix<float> lUVImg(img.height,img.width,img.channel);

	for(ulong i=0; i<img.height; i++)
	{
		for(ulong j=0; j<img.width; j++)
		{
			Pixel<double> point = img.get<double>(i,j);
			Pixel<double> luv;
			RGBtoLUV(point.red, point.green, point.blue, luv.red, luv.green, luv.blue);	

			if(lUVImg.channel > 0) lUVImg.at(i,j,0) = (float)luv.red;
			if(lUVImg.channel > 1) lUVImg.at(i,j,1) = (float)luv.green;
			if(lUVImg.channel > 2) lUVImg.at(i,j,2) = (float)luv.blue;
			if(lUVImg.channel > 3) lUVImg.at(i,j,3) = (float)luv.alpha;
		}
	}

	return lUVImg;
}



void  ColorConversion::RGBtoLUV(const Image& img, Matrix<float>& dest)
{
	assert(img.data && dest.data);
	assert(img.channel == 3 && dest.channel == img.channel);
	assert(img.height > 0 && dest.row == img.height);
	assert(img.width > 0 && dest.column == img.width);

	for(ulong i=0; i<img.height; i++)
	{
		for(ulong j=0; j<img.width; j++)
		{
			Pixel<double> point = img.get<double>(i,j);
			Pixel<double> luv;
			RGBtoLUV(point.red, point.green, point.blue, luv.red, luv.green, luv.blue);	

			if(dest.channel > 0) dest.at(i,j,0) = (float)luv.red;
			if(dest.channel > 1) dest.at(i,j,1) = (float)luv.green;
			if(dest.channel > 2) dest.at(i,j,2) = (float)luv.blue;
			if(dest.channel > 3) dest.at(i,j,3) = (float)luv.alpha;
		}
	}
}


void ColorConversion::LUVtoRGB(const Matrix<float> &luvMat,Image &dest)
{
	assert(dest.data && luvMat.data);
	assert(dest.channel == 3 && luvMat.channel == dest.channel);
	assert(dest.height > 0 && luvMat.row == dest.height);
	assert(dest.width > 0 && luvMat.column == dest.width);

	for(ulong i=0; i<luvMat.row; i++)
	{
		for(ulong j=0; j<luvMat.column; j++)
		{
			Pixel<double> luv;

			if(luvMat.channel > 0) luv.red = luvMat.at(i,j,0);
			if(luvMat.channel > 1) luv.green = luvMat.at(i,j,1);
			if(luvMat.channel > 2) luv.blue = luvMat.at(i,j,2);

			Pixel<double> point;
			LUVtoRGB(luv.red, luv.green, luv.blue,point.red,point.green,point.blue);
			dest.set(i,j,point);
		}
	}
}