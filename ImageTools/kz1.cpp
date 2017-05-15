#include "StereoMatch.h"

/*
if d1 and d2 are two disparities for a pixel p in the left image,
then is_blocked(d1, d2) is true if (p,d2) is closer to the camera
than (p,d1)
*/
#define is_blocked(d1, d2) (KZ1_visibility && (d1 > d2)) /* true only for stereo */

#define VAR_ACTIVE ((Energy::Var)0)

#define KZ1_ALPHA_SINK
/*
if KZ1_ALPHA_SINK is defined then interpretation of a cut is as follows:
SOURCE means initial label
SINK   means new label \alpha

if KZ1_ALPHA_SINK is not defined then SOURCE and SINK are swapped
*/
#ifdef KZ1_ALPHA_SINK
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

inline double StereoMatch::KZ1_dataPenalty(const Point2D<long>& l, const Point2D<long>& r)
{
	register double v = 0;
	if (params.sub_pixel)
	{
		v = params.denominator*dataPenaltySubpixel(l, r);
	} 
	else
	{
		v = params.denominator*dataPenalty(l, r);
		//v = params.denominator*dataCost(l, r);
	}
	v -= params.K;
	if (v>0) v = 0;
	return v;
}


double StereoMatch::KZ1_ComputeEnergy()
{
	Point2D<long> p, q;
	long d = 0, dq = 0;

	E = 0;

	for (p.y=0; p.y<imgSize.y; p.y++)
		for (p.x=0; p.x<imgSize.x; p.x++)
		{
			/* left image and data penalty */
			d = leftX.at(p.y,p.x);

			if (d == OCCLUDED) E += KZ1_OCCLUSION_PENALTY;
			else
			{
				q.x = p.x + d; q.y = p.y;
				if (q >= Point2D<long>(0,0) && q < imgSize)
				{
					dq = rightX.at(q.y,q.x);
					if (d == -dq)
					{
						E += KZ1_dataPenalty(p, q);
					}

					if(!RUNTIMECHECK) //release do that
					{
						/* check visibility constaint */
						if (dq != OCCLUDED && is_blocked(-dq, d))
						{
							fprintf(stderr, "KZ1: Visibility constraint is violated!\n");
							exit(1);
						}
					}
				}
			}

			for (auto k=0; k<NEIGHBOR_NUM; k++)
			{
				q = p + NEIGHBORS[k];

				if (q >= Point2D<long>(0,0) && q < imgSize)
				{
					dq = leftX.at(q.y,q.x);
					E += smoothPenalty(labL,labelsL,p, q, d, dq);
				}
			}


			/* right image */
			d = rightX.at(p.y,p.x);

			if (d == OCCLUDED)
			{
				E += KZ1_OCCLUSION_PENALTY;
			}
			else if(!RUNTIMECHECK) //release do that
			{
				/* check visibility constaint */
				q.x = p.x + d; q.y = p.y;
				if (q >= Point2D<long>(0,0) && q < imgSize)
				{
					dq = leftX.at(q.y,q.x);
					if (dq != OCCLUDED && is_blocked(dq, -d))
					{
						fprintf(stderr, "KZ1: Visibility constraint is violated!\n");
						exit(1);
					}
				}
			}

			for (auto k=0; k<NEIGHBOR_NUM; k++)
			{
				q = p + NEIGHBORS[k];

				if (q >= Point2D<long>(0,0) && q < imgSize)
				{
					dq = rightX.at(q.y,q.x);
					E += smoothPenalty(labR,labelsR,p, q, d, dq);
				}
			}
		}

		return E;
}


/* computes the minimum a-expansion configuration */
void StereoMatch::KZ1_Expand(long a)
{
	Point2D<long> p;
	Energy::Var var, qvar;
	double E_old, delta, E00, E0a, Ea0, Eaa;

	Energy *e = new Energy;

	/* initializing */
	for (p.y = 0; p.y < imgSize.y; p.y++)
	{
		for (p.x=0; p.x < imgSize.x; p.x++)
		{
			long d = leftX.at(p.y,p.x);
			if (a == d) ptrImL.at(p.y,p.x) = VAR_ACTIVE;
			else
			{
				ptrImL.at(p.y,p.x) = var = e -> add_variable();
				if (d == OCCLUDED) e -> ADD_TERM1(var, KZ1_OCCLUSION_PENALTY, 0);
			}

			d = rightX.at(p.y,p.x);
			if (a == -d) ptrImR.at(p.y,p.x) = VAR_ACTIVE;
			else
			{
				ptrImR.at(p.y,p.x) = var = e -> add_variable();
				if (d == OCCLUDED) e -> ADD_TERM1(var, KZ1_OCCLUSION_PENALTY, 0);
			}
		}
	}

	for (p.y=0; p.y<imgSize.y; p.y++)
	{
		for (p.x=0; p.x<imgSize.x; p.x++)
		{
			/* data and visibility terms */
			long d = leftX.at(p.y,p.x);
			long dq = 0;
			var = ptrImL.at(p.y,p.x);
			Point2D<long> q;

			if (d != a && d != OCCLUDED)
			{
				q.x = p.x + d; q.y = p.y;
				if (q >= Point2D<long>(0,0) && q < imgSize)
				{
					qvar = ptrImR.at(q.y,q.x);
					dq = rightX.at(q.y,q.x);
					if (d == -dq)
					{
						delta = (is_blocked(a, d)) ? infinityP : 0;
						e -> ADD_TERM2(var, qvar, KZ1_dataPenalty(p, q), delta, delta, 0);
					}
					else if (is_blocked(a, d))
					{
						e -> ADD_TERM2(var, qvar, 0, infinityP, 0, 0);
					}
				}
			}

			q.x = p.x + a; q.y = p.y;
			if (q >= Point2D<long>(0,0) && q < imgSize)
			{
				qvar = ptrImR.at(q.y,q.x);
				dq = rightX.at(q.y,q.x);

				E0a = (is_blocked(d, a)) ? infinityP : 0;
				Ea0 = (is_blocked(-dq, a)) ? infinityP : 0;
				Eaa = KZ1_dataPenalty(p, q);

				if (var != VAR_ACTIVE)
				{
					if (qvar != VAR_ACTIVE) e -> ADD_TERM2(var, qvar, 0, E0a, Ea0, Eaa);
					else                    e -> ADD_TERM1(var, E0a, Eaa);
				}
				else
				{
					if (qvar != VAR_ACTIVE) e -> ADD_TERM1(qvar, Ea0, Eaa);
					else                    e -> add_constant(Eaa);
				}
			}

			/* left smoothness term */
			for (auto k=0; k<NEIGHBOR_NUM; k++)
			{
				q = p + NEIGHBORS[k];
				if ( ! ( q >= Point2D<long>(0,0) && q < imgSize ) ) continue;
				qvar = ptrImL.at(q.y,q.x);
				dq = leftX.at(q.y,q.x);

				if (var != VAR_ACTIVE && qvar != VAR_ACTIVE)
					E00 = smoothPenalty(labL,labelsL,p, q, d, dq);
				if (var != VAR_ACTIVE)
					E0a = smoothPenalty(labL,labelsL,p, q, d, a);
				if (qvar != VAR_ACTIVE)
					Ea0 = smoothPenalty(labL,labelsL,p, q, a, dq);

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

			/* right smoothness term */
			d = rightX.at(p.y,p.x);
			var = ptrImR.at(p.y,p.x);
			for (auto k=0; k<NEIGHBOR_NUM; k++)
			{
				q = p + NEIGHBORS[k];
				if ( ! ( q >= Point2D<long>(0,0) && q < imgSize ) ) continue;
				qvar = ptrImR.at(q.y,q.x);
				dq = rightX.at(q.y,q.x);

				if (var != VAR_ACTIVE && qvar != VAR_ACTIVE)
					E00 = smoothPenalty(labR,labelsR,p, q, d, dq);
				if (var != VAR_ACTIVE)
					E0a = smoothPenalty(labR,labelsR,p, q, d, -a);
				if (qvar != VAR_ACTIVE)
					Ea0 = smoothPenalty(labR,labelsR,p, q, -a, dq);

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

			/* visibility term */
			if (d != OCCLUDED && is_blocked(a, -d))
			{
				q.x = p.x + d; q.y = p.y;
				if (q >= Point2D<long>(0,0) && q < imgSize)
				{
					if (d != -leftX.at(q.y,q.x))
						e -> ADD_TERM2(var, ptrImL.at(q.y,q.x),0, infinityP, 0, 0);
				}
			}
		}
	}

	E_old = E;
	E = e -> minimize();

	if (E < E_old)
	{
		for (p.y=0; p.y<imgSize.y; p.y++)
		{
			for (p.x=0; p.x<imgSize.x; p.x++)
			{
				var = ptrImL.at(p.y,p.x);

				if (var != VAR_ACTIVE && e->get_var(var)==VALUE1)
				{
					leftX.at(p.y,p.x) = a;
				}

				var = ptrImR.at(p.y,p.x);

				if (var != VAR_ACTIVE && e->get_var(var)==VALUE1)
				{
					rightX.at(p.y,p.x) = -a;
				}
			}
		}
	}

	delete e;
}