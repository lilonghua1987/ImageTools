#pragma once

#include "ImageExpandTools.h"
#include "ColorConversion.h"
#include "tools.h"
#include "energy.h"
#include <vector>
#include <omp.h>
#include "SLIC.h"

extern int gcd(int a, int b);

struct Parameters
{
	/********** data term for CORR, KZ1, KZ2, BVZ **********/
		/*
			if sub_pixel is true then the data term is computed as described in

			Stan Birchfield and Carlo Tomasi
			"A pixel dissimilarity measure that is insensitive to image sampling"
			PAMI 20(4):401-406, April 98

			with one distinction: intensity intervals for a pixels
			are computed from 4 neighbors rather than 2.
		*/
		bool			sub_pixel;
		enum { L1, L2, BT } data_cost;
		int				denominator; /* data term is multiplied by denominator.  */
									 /* Equivalent to using lambda1/denominator, */
									 /* lambda2/denominator, K/denominator       */

		/********** smoothness term for KZ1, KZ2, BVZ **********/
		int				I_threshold;  /* intensity threshold for KZ1 and BVZ */
		int				I_threshold2; /* intensity threshold for KZ2 */
		int				interaction_radius; /* 1 for Potts, >1 for truncated linear */
		int				lambda1, lambda2;

		/********** penalty for an assignment being inactive for KZ1, KZ2 **********/
		int				K;

		/********** occlusion penalty for BVZ (usually INFINITY) **********/
		int				occlusion_penalty;

		/********** iteration parameters for KZ1, KZ2, BVZ **********/
		int				iter_max;
		bool			randomize_every_iteration;

		/********** correlation window for CORR **********/
		int				corr_size;

		/***/
		int             scale;
		/**/
		std::string     saveDir;
};


typedef enum
{
	METHOD_KZ1,
	METHOD_KZ2,
	METHOD_BVZ,
} Method;


/* (half of) the neighborhood system
   the full neighborhood system is edges in NEIGHBORS
   plus reversed edges in NEIGHBORS */
const struct Point2D<long> NEIGHBORS[] = {Point2D<long>(-1,0),Point2D<long>(1,0),Point2D<long>(0,-1),Point2D<long>(0,1)};/*{ Point2D<long>(1, 0), Point2D<long>(0, -1) };*/
const static int NEIGHBOR_NUM = (sizeof(NEIGHBORS) / sizeof(Point2D<long>));


class StereoMatch
{
public:
	StereoMatch(const Image& imgL,const Image& imgR);
	~StereoMatch(void);
	/*friend class Parameters;*/

public:
	void SetDispRange(int _disp_base, int _disp_max);
	double GetK();
	void setParameters(const Parameters& params);
	void SaveXLeft(const std::string& fileName, bool flage = false); /* if flag is TRUE then larger */
	void SaveYLeft(const std::string& fileName, bool flage = false); /* disparities are brighter    */
	void SaveScaledXLeft(const std::string& fileName, bool flag = false); /* if flag is TRUE then larger */
	void SaveScaledYLeft(const std::string& fileName, bool flag = false); /* disparities are brighter    */
	void SaveXRight(const std::string& fileName, bool flage = false); /* if flag is TRUE then larger */
	void KZ1();
	void Run_KZ_BVZ(Method method);

	void CORR();
	void KZ2();
	void BVZ();

	/* data penalty functions for CORR, KZ1, KZ2, BVZ */
	double dataPenalty(const Point2D<long>& l, const Point2D<long>& r);
	double dataBTCost(const Point3D<long>& l, const Point3D<long>& r);
	double dataPenaltySubpixel(const Point2D<long>& l, const Point2D<long>& r);

	/* smoothness penalty functions for KZ1, BVZ */
	double  smoothPenalty(const Matrix<float>& img,const std::vector<int>& labels,const Point2D<long>& p, const Point2D<long>& np, long disp, long ndisp);

	inline int inSegment(const std::vector<int>& labels, const Point2D<long>& p, const Point2D<long>& np, long width)
	{
		return abs(labels.at(p.y*width + p.x) - labels.at(np.y*width + np.x));
	};

	/*************** CORR ALGORITHM *************/
	double		CORR_dataPenalty(const Point2D<long>& l, const Point2D<long>& r);
	void		CORR_hor(const Point2D<long>& d, Matrix<long>& v_hor);
	void		CORR_full(Matrix<long>& v_hor, Matrix<long>& v_full);

	/**************** KZ1 ALGORITHM *************/
	double		KZ1_dataPenalty(const Point2D<long>& l, const Point2D<long>& r);
	double		KZ1_ComputeEnergy();			/* computes current energy */
	void		KZ1_Expand(long a);			/* computes the minimum a-expansion configuration */
	bool		KZ1_visibility;		/* defined only for stereo - then visibility constraint is enforced */

	/**************** KZ2 ALGORITHM *************/
	double		KZ2_dataPenalty(const Point2D<long>& l, const Point2D<long>& r);
	double		KZ2_smoothPenalty2(const Point2D<long>& p, const Point2D<long>& np, long disp);
	double		KZ2_ComputeEnergy();			/* computes current energy */
	void		KZ2_Expand(long a);			/* computes the minimum a-expansion configuration * /

	/**************** BVZ ALGORITHM *************/
	double		BVZ_dataPenalty(const Point2D<long>& p, long d);
	double		BVZ_ComputeEnergy();			/* computes current energy */
	void		BVZ_Expand(long a);			/* computes the minimum a-expansion configuration */

	void swapImage();
	void crossCheck();
	void fillOcclusion();
	void makeUnique();

	void ANCC(float theta, int w);

private:
	void initParams();
	void InitSubPixel();
	void setSubPixel(const Image& Im, Image& ImMin, Image& ImMax);
	void RawCosts();
	double dataCost(const Point2D<long>& pL, const Point2D<long>& pR);

public:
	static const int OCCLUDED = 255;
	static const int infinityP = 10000	;	/* infinite capacity */
	static const int CUTOFF = 1000;
	static const int KZ1_OCCLUSION_PENALTY = 1000;

private:
	RunTimer timer;
	int disp_base, disp_max, disp_size;	/* range of disparities */
	Parameters params;
	bool			uniqueFlag;	/* true if current configuration is unique */
	/* (each pixel corresponds to at most one pixel in the other image */

	/********* INTERNAL VARIABLES **************/
	double				E;					/* current energy */
	Matrix<Energy::Var> ptrImL, ptrImR;	/* used for storing variables corresponding to nodes */

	/************************************************************************/
	/*       segmentation parameters                                                               */
	/************************************************************************/
	int numLabels;
	std::vector<int> labelsL,labelsR;
	int mSpcount,mCompactness;
	Matrix<float> labL,labR;

public:
	Image imgL;   // reference frame
	Image imgR;  // matching image frame
	/*Image segmL;
	Image segmR;*/  /* segmentation images */
	Image imgLMin,imgLMax,imgRMin,imgRMax;
	Matrix<long> leftX;
	Matrix<long> rightX;
	Matrix<float> costRaw;
	Matrix<float> smoothCost;
	Matrix<int> dispMat;
	Point2D<long> imgSize;

	const int size;
	const double truncate;

	private:
		double weight(const Matrix<float>& lab, const Point2D<long>& pL,const Point2D<long>& pR)
		{
			double rc = 5.0;
			double rp = (2 * size + 1)/2.0;

			double d = 0.0,dC = 0.0;
			for (ulong i = 0; i < lab.channel; i++)
			{
				d = lab.at(pL.y,pL.x,i) - lab.at(pR.y,pR.x,i);
				dC += d*d;
			}
			dC = sqrt(dC);
			double dP = tools::EuclideanDistance(pL,pR);
			dC/=rc;
			dP/=rp;
			return exp(-(dC+dP));
		};

		double rawCost(const Point2D<long>& pL,const Point2D<long>& pR)
		{
			double d = 0.0,sum = 0.0;

			for (ulong i = 0; i < labL.channel; i++)
			{
				d = labL.at(pL.y,pL.x,i) - labR.at(pR.y,pR.x,i);
				sum += abs(d);
			}

			return tools::Min(truncate,sum);
		};
};


extern void runStereoMatch(const Image& imgL,const Image& imgR, const std::string& dir);
