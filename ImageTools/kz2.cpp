/* kz2.cpp */
/* Vladimir Kolmogorov (vnk@cs.cornell.edu), 2001-2003. */

#include "StereoMatch.h"

/************************************************************/
/************************************************************/
/************************************************************/

inline double StereoMatch::KZ2_dataPenalty(const Point2D<long>& l, const Point2D<long>& r)
{
	return params.denominator * ((params.sub_pixel ? dataPenaltySubpixel(l, r) : dataPenalty(l, r)) - params.K);
}

inline double StereoMatch::KZ2_smoothPenalty2(const Point2D<long>& p, const Point2D<long>& np, long disp)
{
	double d, d_max = 0;
	/*int r = inSegment(labelsL, p, np, imgSize.x) + inSegment(labelsR, Point2D<long>(p.x + disp, p.y), Point2D<long>(np.x + disp, np.y), imgSize.x);
	if (r == 0) r = 1;*/
	for (auto i = 0; i< imgL.channel; i++)
	{
		//d = abs(labL.at(p.y,p.x,i) - labL.at(np.y,np.x,i));
		d = abs(imgL.atVal<uchar>(p.y, p.x, i) - imgL.atVal<uchar>(np.y, np.x, i));
		if (d_max<d) d_max = d;
		//d = abs(labR.at(p.y, p.x + disp, i) - labR.at(np.y, (np.x + disp), i));
		d = abs(imgR.atVal<uchar>(p.y, p.x + disp, i) - imgR.atVal<uchar>(np.y, (np.x + disp), i));
		if (d_max<d) d_max = d;
	}

	//if (r > 5) return params.lambda1;

	if (d_max<params.I_threshold2) return params.lambda1 / 1;
	else                           return params.lambda2 / 1;
}

/************************************************************/
/************************************************************/
/************************************************************/

/* computes current energy */
double StereoMatch::KZ2_ComputeEnergy()
{
	int k;
	Point2D<long> p, np;

	E = 0;

	for (p.y=0; p.y<imgSize.y; p.y++)
	for (p.x=0; p.x<imgSize.x; p.x++)
	{
		long d = leftX.at(p.y,p.x);

		if (d != OCCLUDED) E += KZ2_dataPenalty(p, Point2D<long>(p.x + d,p.y));

		for (k=0; k<NEIGHBOR_NUM; k++)
		{
			np = p + NEIGHBORS[k];

			if (np>=Point2D<long>(0,0) && np<imgSize)
			{
				long nd = leftX.at(np.y,np.x);
				if (d == nd) continue;
				if (d != OCCLUDED && Point2D<long>(np.x + d, np.y) >= Point2D<long>(0, 0) && Point2D<long>(np.x + d, np.y) < imgSize)
					E += KZ2_smoothPenalty2(p, np, d);
				if (nd != OCCLUDED && Point2D<long>(p.x + nd, p.y) >= Point2D<long>(0, 0) && Point2D<long>(p.x + nd, p.y) < imgSize)
					E += KZ2_smoothPenalty2(p, np, nd);
			}
		}
	}

	return E;
}

/************************************************************/
/************************************************************/
/************************************************************/
#define VAR_ACTIVE     ((Energy::Var)0)
#define VAR_NONPRESENT ((Energy::Var)1)
#define IS_VAR(var) ((unsigned)var>1)

#define KZ2_ALPHA_SINK
/*
	if KZ2_ALPHA_SINK is defined then interpretation of a cut is as in the paper:
	for assignments in A^0:
		SOURCE means 1
		SINK   means 0
	for assigments in A^{\alpha}:
		SOURCE means 0
		SINK   means 1

	if KZ2_ALPHA_SINK is not defined then SOURCE and SINK are swapped
*/
#ifdef KZ2_ALPHA_SINK
	#define ADD_TERM1(var, E0, E1) add_term1(var, E0, E1)
	#define ADD_TERM2(var1, var2, E00, E01, E10, E11) add_term2(var1, var2, E00, E01, E10, E11)
	#define VALUE0 0
	#define VALUE1 1
#else
	#define ADD_TERM1(var, E0, E1) add_term1(var, E1, E0)
	#define ADD_TERM2(var1, var2, E00, E01, E10, E11) add_term2(var1, var2, E11, E10, E01, E00)
	#define VALUE0 1
	#define VALUE1 0
#endif

/* computes the minimum a-expansion configuration */
void StereoMatch::KZ2_Expand(long a)
{
	Point2D<long> p, pd, pa, np;
	long d = 0, nd = 0;
	double E_old, delta;
	Energy::Var var_0, var_a, nvar_0, nvar_a;
	int k;
	
	Energy *e = new Energy;

	/* initializing */
	for (p.y=0; p.y<imgSize.y; p.y++)
	for (p.x=0; p.x<imgSize.x; p.x++)
	{
		d = leftX.at(p.y,p.x);
		pd.x = p.x + d; pd.y = p.y;
		if (a == d)
		{
			ptrImL.at(p.y,p.x) = VAR_ACTIVE;
			ptrImR.at(p.y,p.x) = VAR_ACTIVE;

			e -> add_constant(KZ2_dataPenalty(p, pd));
			continue;
		}

		if (d != OCCLUDED)
		{
			ptrImL.at(p.y,p.x) = var_0 = e -> add_variable();
			e -> ADD_TERM1(var_0, KZ2_dataPenalty(p, pd), 0);
		}
		else ptrImL.at(p.y,p.x) = VAR_NONPRESENT;

		pa.x = p.x + a; pa.y = p.y;
		if (pa>=Point2D<long>(0,0) && pa<imgSize)
		{
			ptrImR.at(p.y,p.x) = var_a = e -> add_variable();
			e -> ADD_TERM1(var_a, 0, KZ2_dataPenalty(p, pa));
		}
		else ptrImR.at(p.y,p.x) = VAR_NONPRESENT;
	}

	for (p.y=0; p.y<imgSize.y; p.y++)
	for (p.x=0; p.x<imgSize.x; p.x++)
	{
		d = leftX.at(p.y, p.x);
		var_0 = ptrImL.at(p.y,p.x);
		var_a = ptrImR.at(p.y,p.x);

		/* adding smoothness */
		for (k=0; k<NEIGHBOR_NUM; k++)
		{
			np = p + NEIGHBORS[k];
			if ( ! ( np>=Point2D<long>(0,0) && np<imgSize ) ) continue;
			nd = leftX.at(np.y, np.x);
			nvar_0 = ptrImL.at(np.y,np.x);
			nvar_a = ptrImR.at(np.y,np.x);

			/* disparity a */
			if (var_a!=VAR_NONPRESENT && nvar_a!=VAR_NONPRESENT)
			/* p+a and np+a are inside the right image */
			{
				delta = KZ2_smoothPenalty2(p, np, a);

				if (var_a != VAR_ACTIVE)
				{
					if (nvar_a != VAR_ACTIVE) e -> ADD_TERM2(var_a, nvar_a, 0, delta, delta, 0);
					else                      e -> ADD_TERM1(var_a, delta, 0);
				}
				else
				{
					if (nvar_a != VAR_ACTIVE) e -> ADD_TERM1(nvar_a, delta, 0);
					else                      {}
				}
			}

			/* disparity d (unless it was checked before) */
			if (IS_VAR(var_0) && Point2D<long>(np.x + d, np.y) >= Point2D<long>(0, 0) && Point2D<long>(np.x + d, np.y) < imgSize)
			{
				delta = KZ2_smoothPenalty2(p, np, d);
				if (d == nd) e -> ADD_TERM2(var_0, nvar_0, 0, delta, delta, 0);
				else e -> ADD_TERM1(var_0, delta, 0);
			}

			/* direction nd (unless it was checked before) */
			if (IS_VAR(nvar_0) && d != nd && Point2D<long>(p.x + nd, np.y) >= Point2D<long>(0, 0) && Point2D<long>(p.x + nd, np.y) < imgSize)
			{
				delta = KZ2_smoothPenalty2(p, np, nd);
				e -> ADD_TERM1(nvar_0, delta, 0);
			}
		}

		/* adding hard constraints in the left image */
		if (IS_VAR(var_0) && var_a!=VAR_NONPRESENT)
			e -> ADD_TERM2(var_0, var_a, 0, infinityP, 0, 0);

		/* adding hard constraints in the right image */
		d = rightX.at(p.y, p.x);
		if (d != OCCLUDED)
		{
			var_0 = ptrImL.at(p.y, (p.x + d));
			if (var_0 != VAR_ACTIVE)
			{
				pa.x = p.x - a; pa.y = p.y;
				if (pa>=Point2D<long>(0,0) && pa<imgSize)
				{
					var_a = ptrImR.at(pa.y,pa.x);
					e -> ADD_TERM2(var_0, var_a, 0, infinityP, 0, 0);
				}
			}
		}
	}

	E_old = E;
	E = e -> minimize();

	if (E < E_old)
	{
		for (p.y=0; p.y<imgSize.y; p.y++)
		for (p.x=0; p.x<imgSize.x; p.x++)
		{
			rightX.at(p.y,p.x) = OCCLUDED;
		}

		for (p.y=0; p.y<imgSize.y; p.y++)
		for (p.x=0; p.x<imgSize.x; p.x++)
		{
			var_0 = ptrImL.at(p.y,p.x);
			var_a = ptrImR.at(p.y,p.x);
			if ( (IS_VAR(var_0) && e->get_var(var_0)==VALUE0) ||
			     (var_0==VAR_ACTIVE) )
			{
				d = leftX.at(p.y, p.x);
				pd.x = p.x + d; pd.y = p.y;
				rightX.at(pd.y,pd.x) = -d;
			}
			else if (IS_VAR(var_a) && e->get_var(var_a)==VALUE1)
			{
				pa.x = p.x + a; pa.y = p.y;
				leftX.at(p.y,p.x) = a;
				rightX.at(pa.y,pa.x) = -a;
			}
			else leftX.at(p.y,p.x) = OCCLUDED;
		}
	}

	delete e;
}

