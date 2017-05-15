#pragma once

#define _CRT_SECURE_NO_WARNINGS

#include "ImageTools.h"
#include "Log.h"
#include <math.h>
#include <assert.h>
#include "Matrix.h"
#include "tools.h"
#include <thread>

#ifndef INT_64
typedef long long INT_64;
#endif

class ImageExpandTools
{
public:
	ImageExpandTools(void);
	~ImageExpandTools(void);
	
	template <typename T1,typename T2>
	static double EuclideanDistance(const Pixel<T1> &p1, const Pixel<T2> &p2);
	template <typename T1,typename T2>
	static double EuclideanDistance(const Point2D<T1> &p1, const Point2D<T2> &p2);
	template <typename T1,typename T2>
	static double EuclideanDistance(const Point3D<T1> &p1, const Point3D<T2> &p2);

	static void DiffHorizontal(const Image& img, Matrix<double>& diff);
	static void DiffVertical(const Image& img, Matrix<double>& diff);

	static Matrix<Pixel<uchar> > Image2Mat(const Image &img);
	template<class T>
	static void Image2Mat(const Image &img,Matrix<T>& mat);
	static Matrix<float> Image2FloatMat(const Image &img);
	template<typename T>
	static Image Mat2Image(const Matrix<Pixel<T> > &mat);

	template<class T>
	static Image Mat2Image(const Matrix<T>& mat, uint scale = 0);

	static Image ImageToJPEG(const Image& src);

	template<typename T>
	static Matrix<T> IntegralImage(const Image& src);

	static Matrix<float> logImage(const Image& img);
	static Matrix<float> logImage(const Matrix<float>& mat);

	/*
	* 函数功能：将非标准图片转换为标准图片（主要是将高动态图转换为标准图片）
	* 参数： 
	*      src:输入图片
	*      EffectiveBits：每个通道的有效表示位数（EffectiveBits>0 && EffectiveBits<=16）
	* 返回值：返回变换后的图片
	*/
	static Image Image2Standard(const Image& src,int EffectiveBits);
	
	static Image InPaint(Image& src ,Image& mask ,uint r);

	static void reSizeImage(const string dir, const string fileExtension, const string saveDir, Point2D<int> top, Point2D<int> bottom);
	static void inPaintAllImages(const string dir, const string fileExtension, const string saveDir);

	static Matrix<float> loadFromFile(std::string fileName);

	static void imShow(const string& windowName, const Image& img);
};


template <typename T1,typename T2>
inline double ImageExpandTools::EuclideanDistance(const Pixel<T1> &p1, const Pixel<T2> &p2)
{
	double d1 = p1.red - p2.red;
	double d2 = p1.green - p2.green;
	double d3 = p1.blue - p2.blue;
	double d4 = p1.alpha - p2.alpha;
	return sqrt(d1*d1 + d2*d2 + d3*d3 + d4*d4);
}


template <typename T1,typename T2>
inline double ImageExpandTools::EuclideanDistance(const Point2D<T1> &p1, const Point2D<T2> &p2)
{
	double d1 = p1.x - p2.x;
	double d2 = p1.y - p2.y;
	return sqrt(d1*d1 + d2*d2);
}


template <typename T1,typename T2>
inline double ImageExpandTools::EuclideanDistance(const Point3D<T1> &p1, const Point3D<T2> &p2)
{
	double d1 = p1.x - p2.x;
	double d2 = p1.y - p2.y;
	double d3 = p1.z - p2.z;
	return sqrt(d1*d1 + d2*d2 + d3*d3);
}



template<class T>
void ImageExpandTools::Image2Mat(const Image &img,Matrix<T>& mat)
{
	assert((img.channel == 1) || (img.channel == 3) || (img.channel == 4));
	assert(mat.channel == img.channel);

	for (ulong i = 0; i < img.height; i++)
	{
		for (ulong j = 0; j < img.width; j++)
		{
			Pixel<T> p = img.getPixel<T>(i, j);
			if (img.channel > 0) mat.at(i, j, 0) = p.red;
			if (img.channel > 1) mat.at(i, j, 1) = p.green;
			if (img.channel > 2) mat.at(i, j, 2) = p.blue;
			if (img.channel > 3) mat.at(i, j, 3) = p.alpha;

			/*mat.at(i, j, c) = img.atVal<T>(i, j, c);*/
		}
	}

}



template<class T>
Image ImageExpandTools::Mat2Image(const Matrix<Pixel<T> > &mat)
{
	assert((mat.channel == 1) || (mat.channel == 3) || (mat.channel == 4));
	Image img(mat.column, mat.row, IMAGE_TYPE::FIT_BITMAP, mat.channel);

	for (ulong i = 0; i < mat.row; i++)
	{
		for (ulong j = 0; j < mat.column; j++)
		{
			ImageTools::setColorPixel(img,i, j,mat.get(i,j));
		}
	}

	ImageTools::NormalImage(img,0,255);

	return img;
}


template<class T>
Image ImageExpandTools::Mat2Image(const Matrix<T>& mat, uint scale)
{
	assert((mat.channel == 1) || (mat.channel == 3) || (mat.channel == 4));
	Image img(mat.column, mat.row, IMAGE_TYPE::FIT_BITMAP, (uchar)mat.channel);

	for (ulong i = 0; i < mat.row; i++)
	{
		for (ulong j = 0; j < mat.column; j++)
		{
			for (uchar c = 0; c < (uchar)mat.channel; c++)
			{
				if (scale) img.atVal<uchar>(i, j, c) = static_cast<uchar>(mat.at(i, j, c) * scale);
			}
		}
	}
	if (scale == 0) ImageTools::NormalImage(img, 0, 255);
	return img;
}


template<typename T>
Matrix<T> ImageExpandTools::IntegralImage(const Image& src)
{
	Matrix<T> result(src.height,src.width,src.channel);
	Pixel<BYTE> point;

	for (ulong r = 0; r < src.height; r++)
	{
		for (ulong c = 0; c < src.width; c++)
		{
			point = src.getPixel<BYTE>(r,c);

			if ((r == c) && (r == 0))
			{
				result.at(r,c,0) = (T)point.red;
				if(result.channel >1) result.at(r,c,1) = (T)point.green;
				if(result.channel >2) result.at(r,c,2) = (T)point.blue;
				if(result.channel >3) result.at(r,c,3) = (T)point.alpha;
			} 
			else
			{
				if (c == 0)
				{
					T t1 = (T)result.at(r-1,c,0);
					result.at(r,c,0) = t1 + (T)point.red;

					if (result.channel >1)
					{
						T t2 = (T)result.at(r-1,c,1);
						result.at(r,c,1) = t2 + (T)point.green;
					}

					if (result.channel >2)
					{
						T t3 = (T)result.at(r-1,c,2);
						result.at(r,c,2) = t3 + (T)point.blue;
					}

					if (result.channel >3)
					{
						T t4 = (T)result.at(r-1,c,3);
						result.at(r,c,3) = t4 + (T)point.alpha;
					}
				} 
				else
				{
					if (r == 0)
					{
						T t1 = (T)result.at(r,c-1,0);
						result.at(r,c,0) = t1 + (T)point.red;

						if (result.channel >1)
						{
							T t2 = (T)result.at(r,c-1,1);
							result.at(r,c,1) = t2 + (T)point.green;
						}

						if (result.channel >2)
						{
							T t3 = (T)result.at(r,c-1,2);
							result.at(r,c,2) = t3 + (T)point.blue;
						}

						if (result.channel >3)
						{
							T t4 = (T)result.at(r,c-1,3);
							result.at(r,c,3) = t4 + (T)point.alpha;
						}
					} 
					else
					{
						T t1 = (T)result.at(r,c-1,0);
						T t2 = (T)result.at(r-1,c,0);
						T t3 = (T)result.at(r-1,c-1,0);
						result.at(r,c,0) = t1 + t2 - t3 + (T)point.red;

						if (result.channel >1)
						{
							t1 = (T)result.at(r,c-1,1);
							t2 = (T)result.at(r-1,c,1);
							t3 = (T)result.at(r-1,c-1,1);
							result.at(r,c,1) = t1 + t2 - t3 + (T)point.green;
						}

						if (result.channel >2)
						{
							t1 = (T)result.at(r,c-1,2);
							t2 = (T)result.at(r-1,c,2);
							t3 = (T)result.at(r-1,c-1,2);
							result.at(r,c,2) = t1 + t2 - t3 + (T)point.blue;
						}

						if (result.channel >3)
						{
							t1 = (T)result.at(r,c-1,3);
							t2 = (T)result.at(r-1,c,3);
							t3 = (T)result.at(r-1,c-1,3);
							result.at(r,c,3) = t1 + t2 - t3 + (T)point.alpha;
						}
					}
				}
			}
		}
	}

	return result;
}