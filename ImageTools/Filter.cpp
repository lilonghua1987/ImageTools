#include "Filter.h"
#include <fstream>


namespace Filter
{
	Matrix<double> boxFilter(const Matrix<double>& src, uint r)
	{
		if (!src.data) throw(new std::exception("boxFilter: the boxFilter data mustn't be empty !"));
		if (r < 0) throw(new std::exception("boxFilter: the box radius must be greater than zeros !"));

		Matrix<double> mDst(src.row,src.column,src.channel);

		Matrix<double> cum(src.cumulativeRow());

		for (ulong v = 0; v < src.row; v++)
		{
			for (ulong u = 0; u < src.column; u++)
			{
				for (ulong c = 0; c < src.channel; c++)
				{
					if (v <= r)
					{
						mDst.at(v, u, c) = cum.at(v + r, u, c);
					}
					else if (v < src.row - r)
					{
						mDst.at(v, u, c) = cum.at(v + r, u, c) - cum.at(v - r - 1, u, c);
					}
					else
					{
						mDst.at(v, u, c) = cum.at(cum.row - 1, u, c) - cum.at(v - r - 1, u, c);
					}
				}
			}
		}

		cum = mDst.cumulativeCol();

		for (ulong u = 0; u < cum.column; u++)
		{
			for (ulong v = 0; v < cum.row; v++)
			{
				for (ulong c = 0; c < cum.channel; c++)
				{
					if (u <= r)
					{
						mDst.at(v, u, c) = cum.at(v, u + r, c);
					}
					else if (u < cum.column - r)
					{
						mDst.at(v, u, c) = cum.at(v, u + r, c) - cum.at(v, u - r - 1, c);
					}
					else
					{
						mDst.at(v, u, c) = cum.at(v, cum.column - 1, c) - cum.at(v, u - r - 1, c);
					}
				}
			}
		}

		return mDst;
	}

	GuidedFilter::GuidedFilter()
		: I(Matrix<double>())
		, r(0)
		, eps(0)
		, N(Matrix<double>())
		, meanR(Matrix<double>())
		, meanG(Matrix<double>())
		, meanB(Matrix<double>())
	{}

	GuidedFilter::GuidedFilter(const Matrix<double>& I, uint r, double eps)
		: I(I)
		, r(r)
		, eps(eps)
		, N(Matrix<double>())
		, meanR(Matrix<double>())
		, meanG(Matrix<double>())
		, meanB(Matrix<double>())
	{
		guidedFilterPre(I,r,eps);
	}

	GuidedFilter::~GuidedFilter(){}

	void GuidedFilter::guidedFilterPre(const Matrix<double>& I, uint r, double eps)
	{
		assert(I.channel == 3);

		N = boxFilter(Matrix<double>::ones(I.row,I.column),r);

		if (!partR.data && I.channel > 0) partR = Matrix<double>::getChannel(I, 0);
		if (!partG.data && I.channel > 1) partG = Matrix<double>::getChannel(I, 1);
		if (!partB.data && I.channel > 2) partB = Matrix<double>::getChannel(I, 2);

		if (!meanR.data && I.channel > 0) meanR = boxFilter(partR, r).dotDivsion(N);
		if (!meanG.data && I.channel > 1) meanG = boxFilter(partG, r).dotDivsion(N);
		if (!meanB.data && I.channel > 2) meanB = boxFilter(partB, r).dotDivsion(N);

		Matrix<double> meanRR = boxFilter(partR.dotMultip(partR), r).dotDivsion(N) - meanR.dotMultip(meanR);
		Matrix<double> meanRG = boxFilter(partR.dotMultip(partG), r).dotDivsion(N) - meanR.dotMultip(meanG);
		Matrix<double> meanRB = boxFilter(partR.dotMultip(partB), r).dotDivsion(N) - meanR.dotMultip(meanB);
		Matrix<double> meanGG = boxFilter(partG.dotMultip(partG), r).dotDivsion(N) - meanG.dotMultip(meanG);
		Matrix<double> meanGB = boxFilter(partG.dotMultip(partB), r).dotDivsion(N) - meanG.dotMultip(meanB);
		Matrix<double> meanBB = boxFilter(partB.dotMultip(partB), r).dotDivsion(N) - meanB.dotMultip(meanB);

		Matrix<double> sigma(3,3);
		invSigma = Matrix<Matrix<double>>(I.row,I.column);

		Matrix<double> e = Matrix<double>::eye(3) * eps;

		for (ulong v = 0; v < I.row; v++)
		{
			for (ulong u = 0; u < I.column; u++)
			{
				sigma.at(0, 0) = meanRR.at(v, u);
				sigma.at(0, 1) = meanRG.at(v, u);
				sigma.at(0, 2) = meanRB.at(v, u);

				sigma.at(1, 0) = meanRG.at(v, u);
				sigma.at(1, 1) = meanGG.at(v, u);
				sigma.at(1, 2) = meanGB.at(v, u);

				sigma.at(2, 0) = meanRB.at(v, u);
				sigma.at(2, 1) = meanGB.at(v, u);
				sigma.at(2, 2) = meanBB.at(v, u);

				sigma = sigma + e;

				invSigma.at(v, u) = Matrix<double>::inverse(sigma);
				/*if (sigma.inv()) invSigma.at(v, u) = sigma;
				else throw(new std::exception("guidedFilterPre: the sigma inv is error !"));*/
			}
		}

		/*std::ofstream f("temp/invSigma.txt");
		f << invSigma;*/
	}

	Matrix<double> GuidedFilter::guidedFilterRun(const Matrix<double>& p)
	{
		assert(p.channel == 1 && p.row == I.row && p.column == I.column);

		Matrix<double> meanP = boxFilter(p, r).dotDivsion(N);

		Matrix<double> covPR = boxFilter(partR.dotMultip(p), r).dotDivsion(N) - meanR.dotMultip(meanP);
		Matrix<double> covPG = boxFilter(partG.dotMultip(p), r).dotDivsion(N) - meanG.dotMultip(meanP);
		Matrix<double> covPB = boxFilter(partB.dotMultip(p), r).dotDivsion(N) - meanB.dotMultip(meanP);

		Matrix<double> cov(1, 3);
		Matrix<double> temp(p.row,p.column,3);

		for (ulong v = 0; v < p.row; v++)
		{
			for (ulong u = 0; u < p.column; u++)
			{
				cov.at(0, 0) = covPR.at(v, u);
				cov.at(0, 1) = covPG.at(v, u);
				cov.at(0, 2) = covPB.at(v, u);

				cov = cov * invSigma.at(v, u);

				for (uint i = 0; i < 3; i++)
				{
					temp.at(v, u, i) = cov.at(0,i);
				}
			}
		}

		/*std::ofstream f("temp/temp.txt");
		f << temp;*/

		Matrix<double> pR = Matrix<double>::getChannel(temp, 0);
		Matrix<double> pG = Matrix<double>::getChannel(temp, 1);
		Matrix<double> pB = Matrix<double>::getChannel(temp, 2);

		Matrix<double> b = meanP - pR.dotMultip(meanR) - pG.dotMultip(meanG) - pB.dotMultip(meanB);
		return (boxFilter(pR, r).dotMultip(partR) + boxFilter(pG, r).dotMultip(partG) + boxFilter(pB, r).dotMultip(partB) + boxFilter(b, r)).dotDivsion(N);
	}

	Matrix<double> weightedMedianFilter(const Matrix<double>& disp, const Matrix<double>& gImg, uint r, uint disNu, double eps)
	{
		Matrix<double> dispOut(disp.row,disp.column);
		Matrix<double> imgAccum(disp.row, disp.column);

		GuidedFilter filter(gImg,r,eps);

		for (size_t d = 1; d <= disNu; d++)
		{
			imgAccum = imgAccum + filter.guidedFilterRun(disp.equal(d));

			for (ulong v = 0; v < disp.row; v++)
			{
				for (ulong u = 0; u < disp.column; u++)
				{
					if (imgAccum.at(v, u) > 0.5 && dispOut.at(v, u) == 0) dispOut.at(v, u) = d;
				}
			}
		}

		return dispOut;
	}

	Matrix<double> medianFilter(const Matrix<double>& src, uint r)
	{
		// temporary memory
		Matrix<double> D_temp(src);
		Matrix<double> out(src);

		Matrix<double> vals(1, r * 2 + 1);

		// first step: horizontal median filter
		for (ulong u = r; u < src.column - r; u++)
		{
			for (ulong v = r; v < src.row - r; v++)
			{
				if (src.at(v, u) >= 0)
				{
					int j = 0;
					for (ulong u2 = u - r; u2 <= u + r; u2++)
					{
						double temp = src.at(v, u2);
						int i = j - 1;
						while (i >= 0 && vals.at(0, i) > temp)
						{
							vals.at(0, i + 1) = vals.at(0, i);
							i--;
						}
						vals.at(0, i + 1) = temp;
						j++;
					}
					D_temp.at(v, u) = vals.at(0, r);
				}
				else
				{
					D_temp.at(v, u) = src.at(v, u);
				}

			}
		}

		// second step: vertical median filter
		for (ulong u = r; u < src.column - r; u++)
		{
			for (ulong v = r; v < src.row - r; v++)
			{
				if (src.at(v, u) >= 0)
				{
					int j = 0;
					for (ulong v2 = v - r; v2 <= v + r; v2++)
					{
						double temp = D_temp.at(v2, u);
						int i = j - 1;
						while (i >= 0 && vals.at(0, i) > temp)
						{
							vals.at(0, i + 1) = vals.at(0, i);
							i--;
						}
						vals.at(0, i + 1) = temp;
						j++;
					}
					out.at(v, u) = vals.at(0, r);
				}
				else {
					out.at(v, u) = out.at(v, u);
				}
			}
		}

		return out;
	}
}
