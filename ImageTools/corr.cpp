/* corr.cpp */
/* Vladimir Kolmogorov (vnk@cs.cornell.edu), 2001-2003. */

#include "StereoMatch.h"

/************************************************************/
/************************************************************/
/************************************************************/

inline double StereoMatch::CORR_dataPenalty(const Point2D<long>& l, const Point2D<long>& r)
{
	if (params.sub_pixel)
		return params.denominator*dataPenaltySubpixel(l, r);
	else
		return params.denominator*dataPenalty(l, r);
}

/************************************************************/
/************************************************************/
/************************************************************/

#define DEFAULT_V 100

/*
  compute correlation values for disparity d
  using 1D correlation window (along horizontal lines)
*/
void StereoMatch::CORR_hor(const Point2D<long>& d, Matrix<long>& v_hor)
{
	Point2D<long> p, pd;
	int w_size, w2_size;
	int x_start, x_finish;
	double v, v_sum;

	w_size = params.corr_size;
	w2_size = 2*w_size + 1;

	for (p.y=0; p.y<imgSize.y; p.y++)
	{
		pd.y = p.y + d.y;
		if (pd.y<0 || pd.y>=imgSize.y)
		{
			for (p.x=0; p.x<imgSize.x; p.x++) v_hor.at(p.y,p.x) = w2_size*DEFAULT_V;
			continue;
		}

		if (d.x > 0)
		{
			x_start = 0;
			x_finish = imgSize.x - 1 - d.x;
			if (x_start > x_finish) x_finish = x_start - 1;
		}
		else
		{
			x_start = 0 - d.x;
			x_finish = imgSize.x - 1;
			if (x_start > x_finish) x_start = x_finish + 1;
		}
		/* now we can access L(p) and R(p+d) if x_start<=p.x<=x_finish */

		v_sum = w2_size*DEFAULT_V;
		for (p.x=x_start; p.x<x_start+w_size; p.x++)
		{
			pd.x = p.x + d.x;
			v = CORR_dataPenalty(p, pd);
			v_sum += v;
			v_sum -= DEFAULT_V;
		}
		for (; p.x<x_start+w2_size; p.x++)
		{
			pd.x = p.x + d.x;
			v = CORR_dataPenalty(p, pd);
			v_sum += v;
			v_sum -= DEFAULT_V;
			v_hor.at(p.y,p.x-w_size) = (long)v_sum;
		}
		for (; p.x<=x_finish; p.x++)
		{
			pd.x = p.x + d.x;
			v = CORR_dataPenalty(p, pd);
			v_sum += v;
			v = CORR_dataPenalty(Point2D<long>(p.x-w2_size, p.y), Point2D<long>(pd.x-w2_size, pd.y));
			v_sum -= v;
			v_hor.at(p.y,p.x-w_size) = (long)v_sum;
		}
		for (; p.x<=x_finish+w_size; p.x++)
		{
			pd.x = p.x + d.x;
			v_sum += DEFAULT_V;
			v = CORR_dataPenalty(Point2D<long>(p.x-w2_size, p.y), Point2D<long>(pd.x-w2_size, pd.y));
			v_sum -= v;
			v_hor.at(p.y, p.x - w_size) = (long)v_sum;
		}

		if (d.x > 0)
		{
			v_sum = v_hor.at(p.y,x_finish);
			for (p.x=x_finish+1; p.x<imgSize.x; p.x++)
			{
				if (p.x-w_size-1>=0 && p.x-w_size-1+d.x<imgSize.x)
				{
					v = CORR_dataPenalty(Point2D<long>(p.x-w_size-1, p.y), Point2D<long>(p.x-w_size-1+d.x, pd.y));
					v_sum -= v;
				}
				else v_sum -= DEFAULT_V;
				v_sum += DEFAULT_V;
				v_hor.at(p.y, p.x) = (long)v_sum;
			}
		}
		else
		{
			v_sum = v_hor.at(p.y,x_start);
			for (p.x=x_start-1; p.x>=0; p.x--)
			{
				if (p.x+w_size+1<imgSize.x && p.x+w_size+1+d.x>=0)
				{
					v = CORR_dataPenalty(Point2D<long>(p.x+w_size+1, p.y), Point2D<long>(p.x+w_size+1+d.x, pd.y));
					v_sum -= v;
				}
				else v_sum -= DEFAULT_V;
				v_sum += DEFAULT_V;
				v_hor.at(p.y, p.x) = (long)v_sum;
			}
		}
	}
}

/************************************************************/
/************************************************************/
/************************************************************/

/*
  compute correlation values
  using full 2D correlation window
*/
void StereoMatch::CORR_full(Matrix<long>& v_hor, Matrix<long>& v_full)
{
	Point2D<long> p;
	int w_size, w2_size;
	int v, v_sum;

	w_size = params.corr_size;
	w2_size = 2*w_size + 1;

	for (p.x=0; p.x<imgSize.x; p.x++)
	{
		v_sum = w2_size*w2_size*DEFAULT_V;
		for (p.y=0; p.y<w_size; p.y++)
		{
			v = v_hor.at(p.y,p.x);
			v_sum += v;
			v_sum -= w2_size*DEFAULT_V;
		}
		for (; p.y<w2_size; p.y++)
		{
			v = v_hor.at(p.y,p.x);
			v_sum += v;
			v_sum -= w2_size*DEFAULT_V;
			v_full.at(p.y-w_size,p.x) = v_sum;
		}
		for (; p.y<imgSize.y; p.y++)
		{
			v = v_hor.at(p.y,p.x);
			v_sum += v;
			v = v_hor.at( p.y-w2_size,p.x);
			v_sum -= v;
			v_full.at(p.y-w_size,p.x) = v_sum;
		}
		for (; p.y<imgSize.y+w_size; p.y++)
		{
			v_sum += w2_size*DEFAULT_V;
			v = v_hor.at( p.y-w2_size,p.x);
			v_sum -= v;
			v_full.at(p.y-w_size,p.x) = v_sum;
		}
	}
}

/************************************************************/
/************************************************************/
/************************************************************/

void StereoMatch::CORR()
{
	Point2D<long> p, d, pd;
	int w_size, w2_size;
	Matrix<Matrix<long>> V_FULL(1,disp_size);
	Matrix<long> v_full, v_hor;
	int v;

	w_size = params.corr_size;
	w2_size = 2*w_size + 1;

	if (w_size < 0)
	{
		fprintf(stderr, "Error in CORR: wrong parameter!\n");
		exit(1);
	}

	printf("CORR:  corr_size = %d, data_cost = L%d\n",
		w_size, params.data_cost == Parameters::L1 ? 1 : 2);

		for (d.x=disp_base; d.x<=disp_max; d.x++)
		{
			V_FULL.at(0,d.x-disp_base) = Matrix<long>(imgSize.y,imgSize.x);
		}
		v_hor  = Matrix<long>(imgSize.y,imgSize.x);

		/* computing correlation values */
		for (d.x = disp_base; d.x <= disp_max; d.x++)
		{
			/* computing v_hor */
			CORR_hor(d, v_hor);

			/* computing v_full */
			v_full = V_FULL.at(0, d.x - disp_base);
			CORR_full(v_hor, v_full);

#ifdef NOT_DEFINED
			/* sanity check */
			{
				for (p.x = im_base.x; p.x <= im_max.x; p.x++)
				{
					Point2D<long> q, qd;
					int v_sum = 0;
					for (q.x = p.x-w_size; q.x <= p.x+w_size; q.x++)
					{
						qd = q + d;
						if (q.x >= im_base.x && q <= im_max.x && qd.x >= im_base.x && qd.x <= im_max.x)
						{
							v = CORR_value(q, qd);
							v_sum += v;
						}
						else v_sum += DEFAULT_V;
					}
						
				}
					
			}
#endif
		}

		/* finding disparity with minimum correlation value */
		for (p.x = 0; p.x<imgSize.x; p.x++)
		{
			Point2D<long> d_min;
			int v_min = -1;

			for (d.x = disp_base; d.x <= disp_max; d.x++)
			{
				v_full = V_FULL.at(0, d.x - disp_base);
				v = v_full.at(p.y, p.x);
				if (v_min < 0 || v < v_min)
				{
					v_min = v;
					d_min = d;
				}
			}
			leftX.at(p.y, p.x) = d_min.x;

		}
		std::cout << "y = " << p.y << std::endl;

		uniqueFlag = false;
}
