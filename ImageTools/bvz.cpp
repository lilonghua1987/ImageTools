/* bvz.cpp */
/* Vladimir Kolmogorov (vnk@cs.cornell.edu), 2001-2003. */

#include "StereoMatch.h"

/************************************************************/
/************************************************************/
/************************************************************/

inline double StereoMatch::BVZ_dataPenalty(const Point2D<long>& p, long d)
{
	Point2D<long> pd(p);

	if (d == OCCLUDED) return params.denominator*params.occlusion_penalty;
	pd.x += d;
	if (!(pd>=Point2D<long>(0,0) && pd<imgSize)) return 1000;

	register double v = 0;
	if (params.sub_pixel)
	{
		v = params.denominator*dataPenaltySubpixel(p, pd);
	} 
	else
	{
		v = params.denominator*dataPenalty(p, pd);
	}
	return v;
}

/************************************************************/
/************************************************************/
/************************************************************/

/* computes current energy */
double StereoMatch::BVZ_ComputeEnergy()
{
	int k;
	Point2D<long> p, q;

	E = 0;

	for (p.y=0; p.y<imgSize.y; p.y++)
	for (p.x=0; p.x<imgSize.x; p.x++)
	{
		long d = leftX.at(p.y, p.x);

		E += BVZ_dataPenalty(p, d);

		for (k=0; k<NEIGHBOR_NUM; k++)
		{
			q = p + NEIGHBORS[k];

			if (q >= Point2D<long>(0,0) && q < imgSize)
			{
				long dq = leftX.at(q.y, q.x);
				E += smoothPenalty(labL,labelsL,p, q, d, dq);
			}
		}
	}

	return E;
}

/************************************************************/
/************************************************************/
/************************************************************/

#define VAR_ACTIVE ((Energy::Var)0)

#define BVZ_ALPHA_SINK
/*
	if BVZ_ALPHA_SINK is defined then interpretation of a cut is as follows:
		SOURCE means initial label
		SINK   means new label \alpha

	if BVZ_ALPHA_SINK is not defined then SOURCE and SINK are swapped
*/
#ifdef BVZ_ALPHA_SINK
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

void StereoMatch::BVZ_Expand(long a)
{
	Point2D<long> p, q;
	long d = 0, dq = 0;
	Energy::Var var, qvar;
	double E_old, E00, E0a, Ea0;
	int k;

	/* node_vars stores variables corresponding to nodes */

	Energy *e = new Energy;

	/* initializing */
	for (p.y=0; p.y<imgSize.y; p.y++)
	for (p.x=0; p.x<imgSize.x; p.x++)
	{
		d = leftX.at(p.y, p.x);
		if (a == d)
		{
			ptrImL.at(p.y,p.x) = VAR_ACTIVE;
			e -> add_constant(BVZ_dataPenalty(p, d));
		}
		else
		{
			ptrImL.at(p.y,p.x) = var = e -> add_variable();
			e -> ADD_TERM1(var, BVZ_dataPenalty(p, d), BVZ_dataPenalty(p, a));
		}
	}

	for (p.y=0; p.y<imgSize.y; p.y++)
	for (p.x=0; p.x<imgSize.x; p.x++)
	{
		d = leftX.at(p.y, p.x);
		var = ptrImL.at(p.y,p.x);

		/* smoothness term */
		for (k=0; k<NEIGHBOR_NUM; k++)
		{
			q = p + NEIGHBORS[k];
			if ( ! ( q>=Point2D<long>(0,0) && q<imgSize ) ) continue;
			qvar = ptrImL.at(q.y,q.x);
			dq = leftX.at(q.y, q.x);

			if (var != VAR_ACTIVE && qvar != VAR_ACTIVE)
				E00 = smoothPenalty(labL, labelsL,p, q, d, dq);
			if (var != VAR_ACTIVE)
				E0a = smoothPenalty(labL, labelsL,p, q, d, a);
			if (qvar != VAR_ACTIVE)
				Ea0 = smoothPenalty(labL, labelsL,p, q, a, dq);

			if (var != VAR_ACTIVE)
			{
				if (qvar != VAR_ACTIVE) e -> ADD_TERM2(var, qvar, E00, E0a, Ea0, 0);
				else                    e -> ADD_TERM1(var, E0a, 0);
			}
			else
			{
				if (qvar != VAR_ACTIVE) e -> ADD_TERM1(qvar, Ea0, 0);
				else                    {}
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
			var = ptrImL.at(p.y,p.x);
			if (var!=VAR_ACTIVE && e->get_var(var)==VALUE1)
			{
				leftX.at(p.y,p.x) = a;
			}
		}
	}

	delete e;
}
