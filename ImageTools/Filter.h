#pragma once

#include "Image.h"
#include "Matrix.h"

namespace Filter
{
	template<class T>
	void sobel3x3(const Matrix<T>& src, Matrix<T>& outH, Matrix<T>& outV)
	{
		for (unsigned long i = 1; i < src.row - 1; i++)
		{
			for (unsigned long j = 1; j < src.column - 1; j++)
			{
				for (unsigned long c = 0; c < src.channel; c++)
				{
					outH.at(i,j,c) = src.at(i-1,j-1,c) - src.at(i+1,j-1,c) + src.at(i-1,j+1,c) - src.at(i+1,j+1,c) 
						             + 2 * (src.at(i,j,c) - src.at(i,j,c));

					outV.at(i,j,c) = src.at(i-1,j-1,c) - src.at(i-1,j+1,c) + src.at(i+1,j-1,c) - src.at(i+1,j+1,c) 
						+ 2 * (src.at(i,j-1,c) - src.at(i,j+1,c));
				}
			}
		}
	}

	template<class T>
	void sobel3x3(const Image& src, Matrix<T>& outH, Matrix<T>& outV)
	{
		for (unsigned long i = 1; i < src.height - 1; i++)
		{
			for (unsigned long j = 1; j < src.width - 1; j++)
			{
				for (unsigned char c = 0; c < src.channel; c++)
				{
					outV.at(i,j,c) = src.atVal<BYTE>(i+1,j-1,c) - src.atVal<BYTE>(i-1,j-1,c) + src.atVal<BYTE>(i+1,j+1,c) - src.atVal<BYTE>(i-1,j+1,c) 
						           + 2 * (src.atVal<BYTE>(i+1,j,c) - src.atVal<BYTE>(i-1,j,c));

					outH.at(i,j,c) = src.atVal<BYTE>(i-1,j-1,c) - src.atVal<BYTE>(i-1,j+1,c) + src.atVal<BYTE>(i+1,j-1,c) - src.atVal<BYTE>(i+1,j+1,c) 
						           + 2 * (src.atVal<BYTE>(i,j-1,c) - src.atVal<BYTE>(i,j+1,c));
				}
			}
		}
	}

	Matrix<double> boxFilter(const Matrix<double>& src, uint r);

	class GuidedFilter
	{
	public:
		GuidedFilter();
		GuidedFilter(const Matrix<double>& I, uint r, double eps);

		~GuidedFilter();

		void guidedFilterPre(const Matrix<double>& I, uint r, double eps);
		Matrix<double> guidedFilterRun(const Matrix<double>& p);
	public:
		Matrix<double> I;
		uint r;
		double eps;

	private:
		Matrix<double> N;
		Matrix<double> partR, partG, partB;
		Matrix<double> meanR, meanG, meanB;
		Matrix<Matrix<double>> invSigma;
	};

	Matrix<double> weightedMedianFilter(const Matrix<double>& disp, const Matrix<double>& gImg, uint r, uint disNu, double eps = 0.01);

	Matrix<double> medianFilter(const Matrix<double>& src, uint r = 3);
}