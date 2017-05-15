#pragma once
#include <vector>
#include "ImageTools.h"
#include "Filter.h"
//#include "Delaunay/Delaunay3D.h"
#include "tetgen.h"

class Elas
{
public:
	enum setting {ROBOTICS,MIDDLEBURY};

	// parameter settings
	struct parameters {
		int width;
		int height;
		int disp_min;               // min disparity
		int disp_max;               // max disparity
		float   support_threshold;      // max. uniqueness ratio (best vs. second best support match)
		int support_texture;        // min texture for support points
		int candidate_stepsize;     // step size of regular grid on which support points are matched
		int incon_window_size;      // window size of inconsistent support point check
		int incon_threshold;        // disparity similarity threshold for support point to be considered consistent
		int incon_min_support;      // minimum number of consistent support points
		bool    add_corners;            // add support points at image corners with nearest neighbor disparities
		int grid_size;              // size of neighborhood for additional support point extrapolation
		float   beta;                   // image likelihood parameter
		float   gamma;                  // prior constant
		float   sigma;                  // prior sigma
		float   sradius;                // prior sigma radius
		int match_texture;          // min texture for dense matching
		int lr_threshold;           // disparity threshold for left/right consistency check
		float   speckle_sim_threshold;  // similarity threshold for speckle segmentation
		int speckle_size;           // maximal size of a speckle (small speckles get removed)
		int ipol_gap_width;         // interpolate small gaps (left<->right, top<->bottom)
		bool    filter_median;          // optional median filter (approximated)
		bool    filter_adaptive_mean;   // optional adaptive mean filter (approximated)
		bool    postprocess_only_left;  // saves time by not postprocessing the right image
		bool    subsampling;            // saves time by only computing disparities for each 2nd pixel
		// note: for this option D1 and D2 must be passed with size
		//       width/2 x height/2 (rounded towards zero)

		// constructor
		parameters (setting s=ROBOTICS) {

			// default settings in a robotics environment
			// (do not produce results in half-occluded areas
			//  and are a bit more robust towards lighting etc.)
			if (s==ROBOTICS) {
				disp_min              = 0;
				disp_max              = 255;
				support_threshold     = 0.85f;
				support_texture       = 10;
				candidate_stepsize    = 5;
				incon_window_size     = 5;
				incon_threshold       = 5;
				incon_min_support     = 5;
				add_corners           = 0;
				grid_size             = 20;
				beta                  = 0.02f;
				gamma                 = 3;
				sigma                 = 1;
				sradius               = 2;
				match_texture         = 1;
				lr_threshold          = 2;
				speckle_sim_threshold = 1;
				speckle_size          = 200;
				ipol_gap_width        = 3;
				filter_median         = 0;
				filter_adaptive_mean  = 1;
				postprocess_only_left = 1;
				subsampling           = 0;

				// default settings for middlebury benchmark
				// (interpolate all missing disparities)
			} else {
				disp_min              = 0;
				disp_max              = 255;
				support_threshold     = 0.95f;
				support_texture       = 10;
				candidate_stepsize    = 5;
				incon_window_size     = 5;
				incon_threshold       = 5;
				incon_min_support     = 5;
				add_corners           = 1;
				grid_size             = 20;
				beta                  = 0.02f;
				gamma                 = 5;
				sigma                 = 1;
				sradius               = 3;
				match_texture         = 0;
				lr_threshold          = 2;
				speckle_sim_threshold = 1;
				speckle_size          = 200;
				ipol_gap_width        = 5000;
				filter_median         = 1;
				filter_adaptive_mean  = 0;
				postprocess_only_left = 0;
				subsampling           = 0;
			}
		}
	};

	struct Triangle 
	{
		int c1,c2,c3;
		float   t1a,t1b,t1c,t1d;
		float   t2a,t2b,t2c,t2d;
		Triangle(int c1,int c2,int c3):c1(c1),c2(c2),c3(c3){}
	};

public:
	Elas(const parameters& param);
	~Elas(void);

	// matching function
	// inputs: pointers to left (I1) and right (I2) intensity image (uint8, input)
	//         pointers to left (D1) and right (D2) disparity image (float, output)
	//         note: D1 and D2 must be allocated before (bytes per line = width)
	//               if subsampling is not active their size is width x height,
	//               otherwise width/2 x height/2 (rounded towards zero)
	void process(const Image& imgL,const Image& imgR,Matrix<float> &D1,Matrix<float> &D2);

private:
	parameters param;

	void createDescriptor(const Matrix<short>& inH, const Matrix<short>& inV,Matrix<short>& dest,const Image& src);

	inline unsigned getAddressOffsetGrid(const int& x,const int& y,const int& d,const int& width,const int& disp_num) {
		return (y*width+x)*disp_num+d;
	}

	// support point functions
	void removeInconsistentSupportPoints(Matrix<short>& D_can);
	void removeRedundantSupportPoints(Matrix<short>& D_can,int redun_max_dist, int redun_threshold, bool vertical);

	void addCornerSupportPoints (std::vector<Point3D<int>> &p_support);
	inline short computeMatchingDisparity (int u,int v,const Matrix<short>& imLS, const Matrix<short>& imRS, bool right_image);
	std::vector<Point3D<int>> computeSupportMatches(const Matrix<short>& imLS, const Matrix<short>& imRS);

	// triangulation & grid
	std::vector<Triangle> computeDelaunayTriangulation (const std::vector<Point3D<int>>& p_support,bool right_image);

	void computeDelaunayTera(const std::vector<Point3D<int>>& p_support,bool right_image);

	void computeDisparityPlanes (const std::vector<Point3D<int>>& p_support,std::vector<Triangle> &tri,bool right_image);
	void planeFitting(const std::vector<Point3D<int>>& p_support,std::vector<Triangle> &tri);
	void createGrid (const std::vector<Point3D<int>>& p_support,Matrix<int>& dispGrid,bool right_image);

	// matching
	inline void updatePosteriorMinimum (const Matrix<short>& imLS, const Matrix<short>& imRS, int u, int v, int ud,int d, float w,float &min_val,int &min_d)
	{
		float val = w;
		for (unsigned long c = 0; c < imLS.channel; c++)
		{
			val += abs(imLS.at(v,u,c) - imRS.at(v,ud,c));
		}
		if (val<min_val) 
		{
			min_val = val;
			min_d   = d;
		}
	};

	inline void updatePosteriorMinimum (const Matrix<short>& imLS, const Matrix<short>& imRS, int u, int v, int ud,int d,float &min_val,int &min_d)
	{
		float val = 0;
		for (unsigned long c = 0; c < imLS.channel; c++)
		{
			val += abs(imLS.at(v,u,c) - imRS.at(v,ud,c));
		}

		if (val < min_val)
		{
			min_val = val;
			min_d   = d;
		}
	};

	inline void findMatch (int u,int v,float plane_a,float plane_b,float plane_c,
		const Matrix<int>& dispGrid,const Matrix<short>& imLS, const Matrix<short>& imRS,
		const Matrix<float>& P,int plane_radius,bool valid,bool right_image,Matrix<float> &D);
	void computeDisparity (const std::vector<Point3D<int>>& p_support,const std::vector<Triangle>& tri,const Matrix<int>& dispGrid,
		const Matrix<short>& imLS, const Matrix<short>& imRS,bool right_image,Matrix<float> &D);
	void computePlaneDisparity (const std::vector<Point3D<int>>& p_support,const std::vector<Triangle>& tri,const Matrix<int>& dispGrid,
		const Matrix<short>& imLS, const Matrix<short>& imRS,bool right_image,Matrix<float> &D);

	// L/R consistency check
	void leftRightConsistencyCheck (Matrix<float> &D1,Matrix<float> &D2);

	// postprocessing
	void removeSmallSegments (Matrix<float> &D);
	void gapInterpolation (Matrix<float> &D);

	// optional postprocessing
	void adaptiveMean (Matrix<float> &D);
	void median (Matrix<float> &D);
	void invalidD(Matrix<float> &D);
	void computeDelaunayTera(Matrix<float> &D);
	void drawDelaunay(const std::vector<Point3D<int>>& p_support, const std::vector<Triangle>& tri, const Image& img, bool rImg);
	inline void findMatch (int u,int v,const Matrix<int>& dispGrid,const Matrix<short>& imLS, const Matrix<short>& imRS,
		int midP, int plane_radius,bool valid,bool right_image,Matrix<float> &D);
};


extern int ElasStereoMatch(const std::string& imL, const std::string& imR);