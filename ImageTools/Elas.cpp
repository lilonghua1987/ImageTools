#include "Elas.h"
#include "triangle.h"
#include "tools.h"
#include "math/plane.h"
#include "meshIO.h"
#include "ImageExpandTools.h"


Elas::Elas(const parameters& param)
	:param(param)
{
}


Elas::~Elas(void)
{
}


void Elas::createDescriptor(const Matrix<short>& inH, const Matrix<short>& inV,Matrix<short>& dest,const Image& src)
{
	for (unsigned long v = 2; v < inH.row - 2; v++)
	{
		for (unsigned long u = 2; u < inH.column - 2; u++)
		{
			dest.at(v,u,0) = inH.at((v - 2),u);
			dest.at(v,u,1) = inH.at((v - 1),u - 2);
			dest.at(v,u,2) = inH.at((v - 1),u);
			dest.at(v,u,3) = inH.at((v - 1),u + 2);
			dest.at(v,u,4) = inH.at(v,u - 1);
			dest.at(v,u,5) = inH.at(v,u);
			dest.at(v,u,6) = inH.at(v,u);
			dest.at(v,u,7) = inH.at(v,u + 1);
			dest.at(v,u,8) = inH.at((v + 1),u - 2);
			dest.at(v,u,9) = inH.at((v + 1),u);
			dest.at(v,u,10) = inH.at((v + 1),u + 2);
			dest.at(v,u,11) = inH.at((v + 2),u);
			dest.at(v,u,12) = inV.at((v - 1),u);
			dest.at(v,u,13) = inV.at(v ,u - 1);
			dest.at(v,u,14) = inV.at(v ,u + 1);
			dest.at(v,u,15) = inV.at((v + 1),u);
			/*dest.at(v,u,16) = src.atVal<uchar>(v,u,0);
			dest.at(v,u,17) = src.atVal<uchar>(v,u,1);
			dest.at(v,u,18) = src.atVal<uchar>(v,u,2);*/
		}
	}
}


std::vector<Point3D<int>> Elas::computeSupportMatches(const Matrix<short>& imLS, const Matrix<short>& imRS)
{
	int D_candidate_stepsize = param.candidate_stepsize;
	if (param.subsampling)
		D_candidate_stepsize += D_candidate_stepsize%2;

	// create matrix for saving disparity candidates
	int D_can_width  = (int)ceil(param.width/(D_candidate_stepsize*1.0f));
	int D_can_height = (int)ceil(param.height/(D_candidate_stepsize*1.0f));
	Matrix<short> dScan(D_can_height,D_can_width);

	// for all point candidates in image 1 do
	for (int u_can = 0; u_can < D_can_width; u_can++) 
	{
		//dScan.at(0,u_can) = -1;
		int u = u_can*D_candidate_stepsize;
		for (int v_can = 0; v_can < D_can_height; v_can++) 
		{
			int v = v_can*D_candidate_stepsize;

			// initialize disparity candidate to invalid
			//dScan.at(v_can,0) = -1;
			dScan.at(v_can,u_can) = -1;

			//For robustness we impose consistency
			// find forwards
			short d = computeMatchingDisparity(u,v,imLS,imRS,false);
			if (d >= 0) 
			{
				// find backwards
				short d2 = computeMatchingDisparity(u-d,v,imLS,imRS,true);
				if (d2 >= 0 && abs(d-d2) <= param.lr_threshold)
					dScan.at(v_can,u_can) = d;
			}
		}
	}

	std::ofstream dScanFile ("temp/dScanMatrix.txt");
	if (dScanFile.is_open()) 
	{
		dScanFile << dScan;
		dScanFile.close();
	}

	// remove inconsistent support points
	removeInconsistentSupportPoints(dScan);

	// remove support points on straight lines, since they are redundant
	// this reduces the number of triangles a little bit and hence speeds up
	// the triangulation process
	removeRedundantSupportPoints(dScan,5,1,true);
	removeRedundantSupportPoints(dScan,5,1,false);

	// move support points from image representation into a vector representation
	std::vector<Point3D<int>> p_support;
	for (int u_can = 1; u_can < D_can_width; u_can++)
		for (int v_can = 1; v_can < D_can_height; v_can++)
			if (dScan.at(v_can,u_can) >= 0)
				p_support.push_back(Point3D<int>(u_can*D_candidate_stepsize,v_can*D_candidate_stepsize,dScan.at(v_can,u_can)));

	// if flag is set, add support points in image corners
	// with the same disparity as the nearest neighbor support point
	/*if (param.add_corners)
		addCornerSupportPoints(p_support);*/

	// return support point vector
	return p_support; 
}


std::vector<Elas::Triangle> Elas::computeDelaunayTriangulation(const std::vector<Point3D<int>>& p_support,bool right_image)
{
	// input/output structure for triangulation
	struct triangulateio in, out;

	// inputs
	in.numberofpoints = p_support.size();
	in.pointlist = (float*)malloc(in.numberofpoints*2*sizeof(float));

	int k = 0;
	for (std::size_t i = 0; i < p_support.size(); i++)
	{
		in.pointlist[k++] = (float)(right_image ? (p_support[i].x-p_support[i].z) : p_support[i].x);
		in.pointlist[k++] = (float)p_support[i].y;
	}

	in.numberofpointattributes = 0;
	in.pointattributelist      = NULL;
	in.pointmarkerlist         = NULL;
	in.numberofsegments        = 0;
	in.numberofholes           = 0;
	in.numberofregions         = 0;
	in.regionlist              = NULL;

	// outputs
	out.pointlist              = NULL;
	out.pointattributelist     = NULL;
	out.pointmarkerlist        = NULL;
	out.trianglelist           = NULL;
	out.triangleattributelist  = NULL;
	out.neighborlist           = NULL;
	out.segmentlist            = NULL;
	out.segmentmarkerlist      = NULL;
	out.edgelist               = NULL;
	out.edgemarkerlist         = NULL;

	// do triangulation (z=zero-based, n=neighbors, Q=quiet, B=no boundary markers)
	char parameters[] = "zQB";
	triangulate(parameters, &in, &out, NULL);

	// put resulting triangles into vector tri
	std::vector<Triangle> tri;
	k=0;
	for (int i=0; i<out.numberoftriangles; i++) 
	{
		tri.push_back(Triangle(out.trianglelist[k],out.trianglelist[k+1],out.trianglelist[k+2]));
		k+=3;
	}

	// free memory used for triangulation
	free(in.pointlist);
	free(out.pointlist);
	free(out.trianglelist);

	// return triangles
	return tri;
}


void Elas::computeDelaunayTera(const std::vector<Point3D<int>>& p_support,bool right_image)
{
	tetgenio in, out;
	// inputs
	in.numberofpoints = p_support.size();
	in.pointlist = new REAL[in.numberofpoints * 3];

	int k = 0;
	for (std::size_t i = 0; i < p_support.size(); i++)
	{
		in.pointlist[k++] = right_image ? (p_support[i].x-p_support[i].z) : p_support[i].x;
		in.pointlist[k++] = p_support[i].y;
		in.pointlist[k++] = p_support[i].z;
	}

	// Output the PLC to files 'barin.node' and 'barin.poly'.
	in.save_nodes("temp/barin");
	in.save_poly("temp/barin");

	// Tetrahedralize the PLC. Switches are chosen to read a PLC (p),
	// do quality mesh generation (q) with a specified quality bound
	// (1.414), and apply a maximum volume constraint (a0.1).
	tetrahedralize(new tetgenbehavior, &in, &out);
	// Output mesh to files 'barout.node', 'barout.ele' and 'barout.face'.
	out.save_nodes("temp/barout");
	out.save_elements("temp/barout");
	out.save_faces("temp/barout");
	out.save_faces2smesh("temp/barout");
}


void Elas::computeDisparityPlanes(const std::vector<Point3D<int>>& p_support,std::vector<Triangle> &tri,bool right_image)
{
	// init matrices
	Matrix<float> A(3,3);
	Matrix<float> b(3,1);

	// for all triangles do
	for (std::size_t i = 0; i < tri.size(); i++) 
	{
		// get triangle corner indices
		int c1 = tri[i].c1;
		int c2 = tri[i].c2;
		int c3 = tri[i].c3;

		// compute matrix A for linear system of left triangle
		A.at(0, 0) = (float)p_support[c1].x;
		A.at(1, 0) = (float)p_support[c2].x;
		A.at(2, 0) = (float)p_support[c3].x;
		A.at(0, 1) = (float)p_support[c1].y; A.at(0, 2) = 1;
		A.at(1, 1) = (float)p_support[c2].y; A.at(1, 2) = 1;
		A.at(2, 1) = (float)p_support[c3].y; A.at(2, 2) = 1;

		// compute vector b for linear system (containing the disparities)
		b.at(0, 0) = (float)p_support[c1].z;
		b.at(1, 0) = (float)p_support[c2].z;
		b.at(2, 0) = (float)p_support[c3].z;

		// on success of gauss jordan elimination
		if (b.solve(A))
		{
			// grab results from b
			tri[i].t1a = b.at(0,0);
			tri[i].t1b = b.at(1,0);
			tri[i].t1c = b.at(2,0);
		} 
		else 
		{
			// otherwise: invalid
			tri[i].t1a = 0;
			tri[i].t1b = 0;
			tri[i].t1c = 0;
		}

		// compute matrix A for linear system of right triangle
		A.at(0, 0) = (float)(p_support[c1].x - p_support[c1].z);
		A.at(1, 0) = (float)(p_support[c2].x - p_support[c2].z);
		A.at(2, 0) = (float)(p_support[c3].x - p_support[c3].z);
		A.at(0, 1) = (float)p_support[c1].y; A.at(0, 2) = 1;
		A.at(1, 1) = (float)p_support[c2].y; A.at(1, 2) = 1;
		A.at(2, 1) = (float)p_support[c3].y; A.at(2, 2) = 1;

		// on success of gauss jordan elimination
		if (b.solve(A))
		{
			// grab results from b
			tri[i].t2a = b.at(0,0);
			tri[i].t2b = b.at(1,0);
			tri[i].t2c = b.at(2,0);
		} 
		else 
		{
			// otherwise: invalid
			tri[i].t2a = 0;
			tri[i].t2b = 0;
			tri[i].t2c = 0;
		}
		
	} 
}


void Elas::planeFitting(const std::vector<Point3D<int>>& p_support,std::vector<Triangle> &tri)
{
	// init matrices
	math::Plane3f planeL(math::Vec3f(0,0,0),0),planeR(math::Vec3f(0,0,0),0);
	math::Vec3f p1(0,0,0),p2(0,0,0),p3(0,0,0);

	// for all triangles do
	for (std::size_t i = 0; i < tri.size(); i++) 
	{
		// get triangle corner indices
		int c1 = tri[i].c1;
		int c2 = tri[i].c2;
		int c3 = tri[i].c3;

		// compute matrix A for linear system of left triangle
		p1 = math::Vector<float, 3>((float)p_support[c1].x, (float)p_support[c1].y, (float)p_support[c1].z);
		p2 = math::Vector<float, 3>((float)p_support[c2].x, (float)p_support[c2].y, (float)p_support[c2].z);
		p3 = math::Vector<float, 3>((float)p_support[c3].x, (float)p_support[c3].y, (float)p_support[c3].z);

		planeL = math::Plane3f(p1,p2,p3);

		tri[i].t1a = planeL.n(0);
		tri[i].t1b = planeL.n(1);
		tri[i].t1c = planeL.n(2);
		tri[i].t1d = planeL.d;

		/*std::cout << p1(0) * planeL.n(0) + p1(1) * planeL.n(1) + p1(2) * planeL.n(2) << std::endl;
		std::cout << p2(0) * planeL.n(0) + p2(1) * planeL.n(1) + p2(2) * planeL.n(2) << std::endl;
		std::cout << p3(0) * planeL.n(0) + p3(1) * planeL.n(1) + p3(2) * planeL.n(2) << std::endl;*/

		// compute matrix A for linear system of right triangle
		p1 = math::Vector<float, 3>((float)(p_support[c1].x - p_support[c1].z), (float)p_support[c1].y, (float)p_support[c1].z);
		p2 = math::Vector<float, 3>((float)(p_support[c2].x - p_support[c2].z), (float)p_support[c2].y, (float)p_support[c2].z);
		p3 = math::Vector<float, 3>((float)(p_support[c3].x - p_support[c3].z), (float)p_support[c3].y, (float)p_support[c3].z);

		planeR = math::Plane3f(p1,p2,p3);

		tri[i].t2a = planeR.n(0);
		tri[i].t2b = planeR.n(1);
		tri[i].t2c = planeR.n(2);
		tri[i].t2d = planeR.d;

		/*std::cout << p1(0) * planeR.n(0) + p1(1) * planeR.n(1) + p1(2) * planeR.n(2) << std::endl;
		std::cout << p2(0) * planeR.n(0) + p2(1) * planeR.n(1) + p2(2) * planeR.n(2) << std::endl;
		std::cout << p3(0) * planeR.n(0) + p3(1) * planeR.n(1) + p3(2) * planeR.n(2) << std::endl;*/
	} 
}


void Elas::createGrid(const std::vector<Point3D<int>>& p_support,Matrix<int>& dispGrid,bool right_image)
{
	Matrix<int> temp1(dispGrid.row,dispGrid.column,param.disp_max+1);
	Matrix<int> temp2(dispGrid.row,dispGrid.column,param.disp_max+1);

	// for all support points do
	for(std::size_t i = 0; i < p_support.size(); i++)
	{
		// compute disparity range to fill for this support point
		int x_curr = p_support[i].x;
		int y_curr = p_support[i].y;
		int d_curr = p_support[i].z;
		int d_min  = max(d_curr-1,0);
		int d_max  = min(d_curr+1,param.disp_max);

		int x = (int)floor((float)(x_curr - (right_image ? d_curr : 0)) / (float)param.grid_size);
		int y = (int)floor((float)y_curr / (float)param.grid_size);

		// fill disparity grid helper
		for (int d = d_min; d <= d_max; d++) 
		{
			// point may potentially lay outside (corner points)
			if (x >= 0 && x < (int)dispGrid.column && y >= 0 && y < (int)dispGrid.row)
			{
				temp1.at(y,x,(ulong)d) = d;
			}
		}
	}
	
	for (int v = 0; v < (int)temp1.row; v++)
	{
		for (int u = 0; u < (int)temp1.column; u++)
		{
			for (int c = 0; c < (int)temp1.channel; c++)
			{
				for (int i = -1; i <= 1; i++)
				{
					for (int j = -1; j <= 1; j++)
					{
						if ((v + i) >= 0 && (v + i) < (int)temp1.row && (u + j) >= 0 && (u + j) < (int)temp1.column)
						{
							int d = temp1.at((v + i), (u + j), c);
							if (d > 0 && std::find(temp2.ptr(v, u), temp2.ptr(v, u, c), d) >= 0)
							{
								temp2.at(v, u, c) = d;
								break;
							}
						}
					}
				}
				/*temp2.at(v-1,u-1,c) = temp1.at(v-2,u-2,c) | temp1.at(v-2,u-1,c) | temp1.at(v-2,u,c)
					                  | temp1.at(v-1,u-2,c) | temp1.at(v-1,u-1,c) | temp1.at(v-1,u,c)
									  | temp1.at(v,u-2,c) | temp1.at(v,u-1,c) | temp1.at(v,u,c);*/
			}
		}
	}
	

	// for all grid positions create disparity grid
	for (unsigned long x = 0; x < dispGrid.column; x++)
	{
		for (unsigned long y = 0; y < dispGrid.row; y++)
		{
			// start with second value (first is reserved for count)
			int curr_ind = 1;

			// for all disparities do
			//for (int d = 0; d <= param.disp_max; d++)
			//{
			//	// if yes => add this disparity to current cell
			//	if (temp2.at(y,x,d) > 0)
			//	{
			//		dispGrid.at(y,x,curr_ind) = d;
			//		curr_ind++;
			//	}
			//}

			for (int c = 0; c < (int)temp2.channel; c++)
			{
				int d = temp2.at(y, x, c);
				if (d > 0)
				{
					dispGrid.at(y, x, curr_ind) = d;
					curr_ind++;
				}
			}

			// finally set number of indices
			dispGrid.at(y,x,0) = curr_ind-1;
		}
	}
}


void Elas::computeDisparity(const std::vector<Point3D<int>>& p_support,const std::vector<Triangle>& tri,const Matrix<int>& dispGrid, const Matrix<short>& imLS, const Matrix<short>& imRS,bool right_image,Matrix<float> &D)
{
	// number of disparities
	const int disp_num  = dispGrid.channel-1;

	// init disparity image to -10
	for (unsigned long r = 0; r < D.row; r++)
	{
		for (unsigned long c = 0; c < D.column; c++)
		{
			for (unsigned long d = 0; d < D.channel; d++)
			{
				D.at(r,c,d) = -10;
			}
		}
	}

	// pre-compute prior 
	float two_sigma_squared = 2*param.sigma*param.sigma;
	Matrix<float> P(disp_num,1);
	for (int delta_d = 0; delta_d < disp_num; delta_d++)
		P.at(delta_d,0) = (-log(param.gamma+exp(-delta_d*delta_d/two_sigma_squared))+log(param.gamma))/param.beta;

	int plane_radius = (int)max((float)ceil(param.sigma*param.sradius),2.0f);

	// for all triangles do
	for (uint i=0; i<tri.size(); i++)
	{
		// get plane parameters
		float plane_a = right_image ? tri[i].t2a : tri[i].t1a;
		float plane_b = right_image ? tri[i].t2b : tri[i].t1b;
		float plane_c = right_image ? tri[i].t2c : tri[i].t1c;
		float plane_d = right_image ? tri[i].t1a : tri[i].t2a;

		// triangle corners
		int c1 = tri[i].c1;
		int c2 = tri[i].c2;
		int c3 = tri[i].c3;

		// sort triangle corners wrt. u (ascending)    
		std::vector<Point2D<float>> pTri(3);
		pTri.at(0).y = (float)p_support[c1].y;
		pTri.at(1).y = (float)p_support[c2].y;
		pTri.at(2).y = (float)p_support[c3].y;

		pTri.at(0).x = (float)(p_support[c1].x - (right_image ? p_support[c1].z : 0));
		pTri.at(1).x = (float)(p_support[c2].x - (right_image ? p_support[c2].z : 0));
		pTri.at(2).x = (float)(p_support[c3].x - (right_image ? p_support[c3].z : 0));

		sort(pTri.begin(),pTri.end(),[](const Point2D<float>& p1,const Point2D<float>& p2){return p1.x < p2.x;});

		// rename corners
		Point2D<float> A = pTri.at(0);
		Point2D<float> B = pTri.at(1);
		Point2D<float> C = pTri.at(2);

		// compute straight lines connecting triangle corners
		float AB_a = 0; float AC_a = 0; float BC_a = 0;
		if ((int)(A.x)!=(int)(B.x)) AB_a = (A.y-B.y)/(A.x-B.x);
		if ((int)(A.x)!=(int)(C.x)) AC_a = (A.y-C.y)/(A.x-C.x);
		if ((int)(B.x)!=(int)(C.x)) BC_a = (B.y-C.y)/(B.x-C.x);
		float AB_b = A.y-AB_a*A.x;
		float AC_b = A.y-AC_a*A.x;
		float BC_b = B.y-BC_a*B.x;

		// a plane is only valid if itself and its projection
		// into the other image is not too much slanted
		bool valid = (abs(plane_a) < 0.7) && (abs(plane_d) < 0.7);

		// first part (triangle corner A->B)
		if ((int)(A.x) != (int)(B.x))
		{
			for (int u = max((int)A.x,0); u < min((int)B.x,param.width); u++)
			{
				if (!param.subsampling || u%2==0)
				{
					int v_1 = (int)(AC_a*(float)u+AC_b);
					int v_2 = (int)(AB_a*(float)u+AB_b);
					for (int v = min(v_1,v_2); v < max(v_1,v_2); v++)
					{
						if (!param.subsampling || v%2==0)
						{
							findMatch(u,v,plane_a,plane_b,plane_c,dispGrid,imLS,imRS,P,plane_radius,valid,right_image,D);
						}
					}
				}
			}
		}

		// second part (triangle corner B->C)
		if ((int)(B.x)!=(int)(C.x))
		{
			for (int u = max((int)B.x,0); u < min((int)C.x,param.width); u++)
			{
				if (!param.subsampling || u%2==0)
				{
					int v_1 = (int)(AC_a*(float)u+AC_b);
					int v_2 = (int)(BC_a*(float)u+BC_b);
					for (int v = min(v_1,v_2); v < max(v_1,v_2); v++)
					{
						if (!param.subsampling || v%2==0)
						{
							findMatch(u,v,plane_a,plane_b,plane_c,dispGrid,imLS,imRS,P,plane_radius,valid,right_image,D);
						}
					}
				}
			}
		}

	}
}


void Elas::computePlaneDisparity(const std::vector<Point3D<int>>& p_support,const std::vector<Triangle>& tri,const Matrix<int>& dispGrid, const Matrix<short>& imLS, const Matrix<short>& imRS,bool right_image,Matrix<float> &D)
{
	// number of disparities
	const int disp_num  = dispGrid.channel-1;

	// init disparity image to -10
	for (unsigned long r = 0; r < D.row; r++)
	{
		for (unsigned long c = 0; c < D.column; c++)
		{
			for (unsigned long d = 0; d < D.channel; d++)
			{
				D.at(r,c,d) = -10;
			}
		}
	}

	// pre-compute prior 
	/*float two_sigma_squared = 2*param.sigma*param.sigma;
	Matrix<float> P(disp_num,1);
	for (int delta_d = 0; delta_d < disp_num; delta_d++)
		P.at(delta_d,0) = (-log(param.gamma+exp(-delta_d*delta_d/two_sigma_squared))+log(param.gamma))/param.beta;*/

	int plane_radius = (int)max((float)ceil(param.sigma*param.sradius),2.0f);

	// for all triangles do
	for (uint i=0; i<tri.size(); i++)
	{
		// get plane parameters
		float plane_a = right_image ? tri[i].t2a : tri[i].t1a;
		float plane_b = right_image ? tri[i].t2b : tri[i].t1b;
		float plane_c = right_image ? tri[i].t2c : tri[i].t1c;
		float plane_d = right_image ? tri[i].t2d : tri[i].t1d;

		math::Plane3f plane(math::Vec3f(plane_a,plane_b,plane_c),plane_d);

		// triangle corners
		int c1 = tri[i].c1;
		int c2 = tri[i].c2;
		int c3 = tri[i].c3;

		// sort triangle corners wrt. u (ascending)    
		std::vector<Point2D<float>> pTri(3);
		pTri.at(0).y = (float)p_support[c1].y;
		pTri.at(1).y = (float)p_support[c2].y;
		pTri.at(2).y = (float)p_support[c3].y;

		pTri.at(0).x = (float)(p_support[c1].x - (right_image ? p_support[c1].z : 0));
		pTri.at(1).x = (float)(p_support[c2].x - (right_image ? p_support[c2].z : 0));
		pTri.at(2).x = (float)(p_support[c3].x - (right_image ? p_support[c3].z : 0));

		sort(pTri.begin(),pTri.end(),[](const Point2D<float>& p1,const Point2D<float>& p2){return p1.x < p2.x;});

		Point2D<float> uR(pTri.at(0).x,pTri.at(2).x);

		sort(pTri.begin(),pTri.end(),[](const Point2D<float>& p1,const Point2D<float>& p2){return p1.y < p2.y;});

		Point2D<float> vR(pTri.at(0).y,pTri.at(2).y);

		math::Vec3f A((float)p_support[c1].x, (float)p_support[c1].y, (float)p_support[c1].z);
		math::Vec3f B((float)p_support[c2].x, (float)p_support[c2].y, (float)p_support[c2].z);
		math::Vec3f C((float)p_support[c3].x, (float)p_support[c3].y, (float)p_support[c3].z);

		for (int u = max((int)uR.x,0); u <= min((int)uR.y,param.width-1); u++)
		{
			for (int v = max((int)vR.x,0); v <= min((int)vR.y,param.height-1); v++)
			{			
				if (plane_c == 0)
				{
					//D.at(v,u) = p_support[c1].z;
					findMatch(u,v,dispGrid,imLS,imRS,p_support[c1].z,plane_radius,true,right_image,D);
				} 
				else
				{
					math::Vec3f pD((float)u, (float)v, (plane_d - plane_a * u - plane_b * v) / plane_c);
					if (plane.PointinTriangle(A,B,C,pD) || plane.point_dist(pD) < 0.001f)
					{
						//D.at(v,u) = P(2);
						findMatch(u,v,dispGrid,imLS,imRS,(int)pD(2),plane_radius,true,right_image,D);
					}
					else
					{
						math::Vec3f qD = plane.projection(pD);

						if (plane.PointinTriangle(A,B,C,qD) || plane.point_dist(qD) < 0.0001f)
						{
							//D.at(v,u) = pD(2);
							findMatch(u,v,dispGrid,imLS,imRS,(int)pD(2),plane_radius,true,right_image,D);
						}
						else std::cout << plane.point_dist(qD) << std::endl;
					}

					//D.at(v,u) = (plane_d - plane_a * u - plane_b * v) / plane_c;
				}
			}
		}

	}
}


inline void Elas::findMatch(int u,int v,float plane_a,float plane_b,float plane_c, const Matrix<int>& dispGrid,const Matrix<short>& imLS,const Matrix<short>& imRS,const Matrix<float>& P,int plane_radius,bool valid,bool right_image,Matrix<float> &D)
{
	// get image width and height
	const int disp_num    = dispGrid.channel - 1;
	// descriptor window_size
	int window_size = 2;

	// check if u is ok
	if (u < window_size || u >= param.width-window_size)
		return;

	int vb = max(min(v,param.height-3),2);

	// does this patch have enough texture?
	int sum = 0;
	for (unsigned long i = 0; i < imRS.channel; i++)
	{
		sum += abs(right_image ? imRS.at(vb,u,i):imLS.at(vb,u,i));
	}
	sum = abs(sum - (int)imRS.channel * 128 * 4);
	if (sum < param.match_texture)
		return;

	// compute disparity, min disparity and max disparity of plane prior
	int d_plane     = (int)ceil(plane_a*(float)u+plane_b*(float)v+plane_c);
	int d_plane_min = max((int)ceil(d_plane-plane_radius),0);
	int d_plane_max = min((int)ceil(d_plane+plane_radius),disp_num-1);

	// get grid pointer
	int  grid_x    = (int)floor((float)u/(float)param.grid_size);
	int  grid_y    = (int)floor((float)v/(float)param.grid_size); 
	int  num_grid  = dispGrid.at(grid_y,grid_x);

	float min_val = std::numeric_limits<float>::max();
	int min_d   = -1;

	//此处处理视差不在拟合得出的视差范围内
	for (int i = 1; i <= num_grid; i++) 
	{
		int d_curr = dispGrid.at(grid_y,grid_x,i);
		if (d_curr < d_plane_min || d_curr > d_plane_max)
		{
			int u_warp = right_image ? (u+d_curr) : (u-d_curr);
			if (u_warp < window_size || u_warp >= param.width-window_size)
				continue;
			if(right_image)
				updatePosteriorMinimum(imRS,imLS,u,vb,u_warp,d_curr,min_val,min_d);
			else
				updatePosteriorMinimum(imLS,imRS,u,vb,u_warp,d_curr,min_val,min_d);
		}
	}

	for (int d_curr = d_plane_min; d_curr <= d_plane_max; d_curr++)
	{
		int u_warp = right_image ? (u+d_curr) : (u-d_curr);
		if (u_warp < window_size || u_warp >= param.width-window_size)
			continue;
		if(right_image)
			updatePosteriorMinimum(imRS,imLS,u,vb,u_warp,d_curr,valid?P.at(abs(d_curr-d_plane),0):0,min_val,min_d);
		else
			updatePosteriorMinimum(imLS,imRS,u,vb,u_warp,d_curr,valid?P.at(abs(d_curr-d_plane),0):0,min_val,min_d);
	}

	// set disparity value
	if (min_d >= 0) D.at(v, u) = (float)min_d; // MAP value (min neg-Log probability)
	else          D.at(v,u) = -1;    // invalid disparity
}


inline void Elas::findMatch(int u,int v, const Matrix<int>& dispGrid,const Matrix<short>& imLS,const Matrix<short>& imRS,int midP, int plane_radius,bool valid,bool right_image,Matrix<float> &D)
{
	// get image width and height
	const int disp_num    = dispGrid.channel - 1;
	// descriptor window_size
	int window_size = 2;

	// check if u is ok
	if (u < window_size || u >= param.width-window_size)
		return;

	int vb = max(min(v,param.height-3),2);

	// does this patch have enough texture?
	int sum = 0;
	for (unsigned long i = 0; i < imRS.channel; i++)
	{
		sum += abs(right_image ? imRS.at(vb,u,i):imLS.at(vb,u,i));
	}
	sum = abs(sum - (int)imRS.channel * 128 * 4);
	if (sum < param.match_texture)
		return;

	// compute disparity, min disparity and max disparity of plane prior
	int d_plane_min = max((int)ceil(midP-plane_radius),0);
	int d_plane_max = min((int)ceil(midP+plane_radius),disp_num-1);

	// get grid pointer
	int  grid_x    = (int)floor((float)u/(float)param.grid_size);
	int  grid_y    = (int)floor((float)v/(float)param.grid_size); 
	int  num_grid  = dispGrid.at(grid_y,grid_x);

	float min_val = std::numeric_limits<float>::max();
	int min_d   = -1;

	//此处处理视差不在拟合得出的视差范围内
	int dMean = 0;
	for (int i = 1; i <= num_grid; i++) 
	{
		int d_curr = dispGrid.at(grid_y,grid_x,i);
		if (d_curr < d_plane_min || d_curr > d_plane_max)
		{
			int u_warp = right_image ? (u+d_curr) : (u-d_curr);
			if (u_warp < window_size || u_warp >= param.width-window_size)
				continue;

			if(right_image)
				updatePosteriorMinimum(imRS,imLS,u,vb,u_warp,d_curr,min_val,min_d);
			else
				updatePosteriorMinimum(imLS,imRS,u,vb,u_warp,d_curr,min_val,min_d);
		}
		dMean += d_curr;
	}

	float dMR = (-log(param.gamma + exp(-dMean*dMean / (2 * param.sigma*param.sigma))) + log(param.gamma)) / param.beta;

	for (int d_curr = d_plane_min; d_curr <= d_plane_max; d_curr++)
	{
		int u_warp = right_image ? (u+d_curr) : (u-d_curr);
		if (u_warp < window_size || u_warp >= param.width-window_size)
			continue;
		if(right_image)
			updatePosteriorMinimum(imRS, imLS, u, vb, u_warp, d_curr, dMR, min_val, min_d);
		else
			updatePosteriorMinimum(imLS, imRS, u, vb, u_warp, d_curr, dMR, min_val, min_d);
	}

	// set disparity value
	if (min_d >= 0) D.at(v, u) = (float)min_d; // MAP value (min neg-Log probability)
	else          D.at(v,u) = -1;    // invalid disparity
}



void Elas::leftRightConsistencyCheck(Matrix<float> &D1,Matrix<float> &D2)
{
	// make a copy of both images
	Matrix<float> D1_copy(D1),D2_copy(D2);
	// for all image points do
	for (unsigned long u = 0; u < D1.column; u++)
	{
		for (unsigned long v = 0; v < D1.row; v++) 
		{
			// compute disparity value
			float d1 = D1_copy.at(v,u);
			float d2 = D2_copy.at(v,u);
			float u_warp_1 = 0,u_warp_2 = 0;

			if (param.subsampling) {
				u_warp_1 = (float)u-d1/2;
				u_warp_2 = (float)u+d2/2;
			} else {
				u_warp_1 = (float)u-d1;
				u_warp_2 = (float)u+d2;
			}


			// check if left disparity is valid
			if (d1 >= 0 && u_warp_1 >= 0 && u_warp_1 < (float)D1.column) 
			{       
				// if check failed
				if (fabs(D2_copy.at(v,(ulong)u_warp_1)-d1) > param.lr_threshold)
					D1.at(v,u) = -10; // set invalid
			}
			else
			{
				D1.at(v,u) = -10; // set invalid
			}

			// check if right disparity is valid
			if (d2 >= 0 && u_warp_2 >= 0 && u_warp_2 < (float)D1.column) 
			{       
				// if check failed
				if (fabs(D1_copy.at(v,(unsigned)u_warp_2)-d2) > param.lr_threshold)
					D2.at(v,u) = -10; // set invalid			
			} 
			else
			{
				D2.at(v,u) = -10; // set invalid
			}
		}
	}
}


void Elas::removeSmallSegments(Matrix<float> &D)
{
	int D_speckle_size = param.subsampling ? (int)(sqrt((float)param.speckle_size) * 2) : param.speckle_size;

	// allocate memory on heap for dynamic programming arrays
	Matrix<int> D_done(D.row,D.column),seg_list_u(D.row,D.column),seg_list_v(D.row,D.column);
	int u_neighbor[4];
	int v_neighbor[4];

	// for all pixels do
	for (unsigned long u = 0; u < D.column; u++) 
	{
		for (unsigned long v = 0; v < D.row; v++) 
		{
			// if this pixel has not already been processed
			if (D_done.at(v,u) == 0)
			{
				// init segment list (add first element
				// and set it to be the next element to check)
				seg_list_u.at(0,0) = u;
				seg_list_v.at(0,0) = v;
				int seg_list_count  = 1;
				int seg_list_curr   = 0;

				// add neighboring segments as long as there
				// are none-processed pixels in the seg_list;
				// none-processed means: seg_list_curr<seg_list_count
				while (seg_list_curr < seg_list_count) 
				{
					// get current position from seg_list
					int u_seg_curr = seg_list_u.at(seg_list_curr/seg_list_u.column,seg_list_curr%seg_list_u.column);
					int v_seg_curr = seg_list_v.at(seg_list_curr/seg_list_v.column,seg_list_curr%seg_list_v.column);

					// fill list with neighbor positions
					u_neighbor[0] = u_seg_curr-1; v_neighbor[0] = v_seg_curr;
					u_neighbor[1] = u_seg_curr+1; v_neighbor[1] = v_seg_curr;
					u_neighbor[2] = u_seg_curr;   v_neighbor[2] = v_seg_curr-1;
					u_neighbor[3] = u_seg_curr;   v_neighbor[3] = v_seg_curr+1;

					// for all neighbors do
					for (int i = 0; i < 4; i++)
					{
						// check if neighbor is inside image
						if (u_neighbor[i] >= 0 && v_neighbor[i] >= 0 && u_neighbor[i] < (int)D.column && v_neighbor[i] < (int)D.row) 
						{
							// check if neighbor has not been added yet and if it is valid
							if (D_done.at(v_neighbor[i],u_neighbor[i]) == 0 && D.at(v_neighbor[i],u_neighbor[i]) >= 0) 
							{
								// is the neighbor similar to the current pixel
								// (=belonging to the current segment)
								if (fabs(D.at(v_seg_curr,u_seg_curr) - D.at(v_neighbor[i],u_neighbor[i])) <= param.speckle_sim_threshold) 
								{

									// add neighbor coordinates to segment list
									seg_list_u.at(seg_list_count/seg_list_u.column,seg_list_count%seg_list_u.column) = u_neighbor[i];
									seg_list_v.at(seg_list_count/seg_list_v.column,seg_list_count%seg_list_v.column) = v_neighbor[i];
									seg_list_count++;            

									// set neighbor pixel in I_done to "done"
									// (otherwise a pixel may be added 2 times to the list, as
									//  neighbor of one pixel and as neighbor of another pixel)
									D_done.at(v_neighbor[i],u_neighbor[i]) = 1;
								}
							}

						} 
					}

					// set current pixel in seg_list to "done"
					seg_list_curr++;

					// set current pixel in I_done to "done"
					D_done.at(v_seg_curr,u_seg_curr) = 1;

				} // end: while (seg_list_curr<seg_list_count)

				// if segment NOT large enough => invalidate pixels
				if (seg_list_count < D_speckle_size) 
				{

					// for all pixels in current segment invalidate pixels
					for (int i=0; i<seg_list_count; i++) 
					{
						D.at(seg_list_v.at(i/seg_list_v.column,i%seg_list_v.column),seg_list_u.at(i/seg_list_u.column,i%seg_list_u.column)) = -10;
					}
				}
			} 
		}
	}
}


void Elas::gapInterpolation(Matrix<float> &D)
{
	int D_ipol_gap_width = param.subsampling ? (param.ipol_gap_width/2+1) : param.ipol_gap_width;

	// discontinuity threshold
	float discon_threshold = 3.0;

	// declare loop variables
	int count,v_first,v_last,u_first,u_last;

	// 1. Row-wise:
	// for each row do
	for (unsigned long v = 0; v < D.row; v++) 
	{
		// init counter
		count = 0;

		// for each element of the row do
		for (unsigned long u = 0; u < D.column; u++) 
		{
			// if disparity valid
			if (D.at(v,u) >= 0) 
			{
				// check if speckle is small enough
				if (count >= 1 && count <= D_ipol_gap_width) 
				{
					// first and last value for interpolation
					u_first = u-count;
					u_last  = u-1;

					// if value in range
					if (u_first > 0 && u_last < (int)D.column-1) 
					{
						// compute mean disparity
						float d1 = D.at(v,u_first-1);
						float d2 = D.at(v,u_last+1);
						float d_ipol = min(d1,d2);
						if (abs(d1-d2) < discon_threshold) d_ipol = (d1+d2)/2;

						// set all values to d_ipol
						for (int u_curr = u_first; u_curr <= u_last; u_curr++)
							D.at(v,u_curr) = d_ipol;
					}
				}
				// reset counter
				count = 0;

				// otherwise increment counter
			} 
			else
			{
				count++;
			}
		}

		// if full size disp map requested
		if (param.add_corners) 
		{
			// extrapolate to the left
			for (ulong u = 0; u < D.column; u++) 
			{
				// if disparity valid
				if (D.at(v,u) >= 0) 
				{
					for (ulong u2 = (ulong)max((int)u - D_ipol_gap_width, 0); u2 < u; u2++)
						 D.at(v,u2) = D.at(v,u);
					break;
				}
			}

			// extrapolate to the right
			for (int u = (int)D.column-1; u >= 0; u--) 
			{
				// if disparity valid
				if (D.at(v,u) >= 0) 
				{
					for (int u2 = u; u2 <= min(u+D_ipol_gap_width,(int)D.column-1); u2++)
						D.at(v,u2) = D.at(v,u);
					break;
				}
			}
		}
	}

	// 2. Column-wise:
	// for each column do
	for (unsigned long u = 0; u < D.column; u++)
	{
		// init counter
		count = 0;

		// for each element of the column do
		for (unsigned long v = 0; v < D.row; v++) 
		{
			// if disparity valid
			if (D.at(v,u) >= 0) 
			{
				// check if gap is small enough
				if (count >= 1 && count <= D_ipol_gap_width) 
				{
					// first and last value for interpolation
					v_first = v-count;
					v_last  = v-1;

					// if value in range
					if (v_first > 0 && v_last < (int)D.row-1) 
					{
						// compute mean disparity
						float d1 = D.at(v_first-1,u);
						float d2 = D.at(v_last+1,u);
						float d_ipol = min(d1,d2);
						if (fabs(d1-d2) < discon_threshold) d_ipol = (d1+d2)/2;                           
						// set all values to d_ipol
						for (int v_curr = v_first; v_curr <= v_last; v_curr++)
							D.at(v_curr,u) = d_ipol;
					}
				}

				// reset counter
				count = 0;

				// otherwise increment counter
			} 
			else 
			{
				count++;
			}
		}

		// if full size disp map requested
		if (param.add_corners) 
		{
			// extrapolate to the top
			for (ulong v = 0; v < D.row; v++) 
			{
				// if disparity valid
				if (D.at(v,u) >= 0) 
				{
					for (ulong v2 = (ulong)max((int)v - D_ipol_gap_width, 0); v2 < v; v2++)
						D.at(v2,u) = D.at(v,u);
					break;
				}
			}

			// extrapolate to the bottom
			for (int v = (int)D.row-1; v >= 0; v--) 
			{
				// if disparity valid
				if (D.at(v,u) >= 0) 
				{
					for (int v2 = v; v2 <= min(v+D_ipol_gap_width,(int)D.row-1); v2++)
						D.at(v2,u) = D.at(v,u);
					break;
				}
			}
		}
	}
}


void Elas::adaptiveMean(Matrix<float> &D)
{
	// allocate temporary memory
	Matrix<float> D_copy(D),D_tmp(D.row,D.column);

	// zero input disparity maps to -10 (this makes the bilateral
	// weights of all valid disparities to 0 in this region)
	for (unsigned long v = 0; v < D.row; v++)
	{
		for (unsigned long u = 0; u < D.column; u++)
		{
			for (unsigned long c = 0; c < D.channel; c++)
			{
				if (D.at(v,u,c) < 0)
				{
					D_copy.at(v,u,c) = -10;
					D_tmp.at(v,u,c) = -10;
				}
			}
		}
	}

	__m128 xconst0 = _mm_set1_ps(0);
	__m128 xconst4 = _mm_set1_ps(4);
	__m128 xval,xweight1,xweight2,xfactor1,xfactor2;

	float *val     = (float *)_mm_malloc(8*sizeof(float),16);
	float *weight  = (float*)_mm_malloc(4*sizeof(float),16);
	float *factor  = (float*)_mm_malloc(4*sizeof(float),16);

	// set absolute mask
	__m128 xabsmask = _mm_set1_ps(0x7FFFFFFF);

	// when doing subsampling: 4 pixel bilateral filter width
	if (param.subsampling)
	{
		// horizontal filter
		for (unsigned long v = 3; v < D.row-3; v++) 
		{
			// init
			for (int u = 0; u < 3; u++)
				val[u] = D_copy.at(v,u);

			// loop
			for (unsigned long u = 3; u < D.column; u++)
			{
				// set
				float val_curr = D_copy.at(v,u-1);
				val[u%4] = D_copy.at(v,u);

				xval     = _mm_load_ps(val);      
				xweight1 = _mm_sub_ps(xval,_mm_set1_ps(val_curr));
				xweight1 = _mm_and_ps(xweight1,xabsmask);
				xweight1 = _mm_sub_ps(xconst4,xweight1);
				xweight1 = _mm_max_ps(xconst0,xweight1);
				xfactor1 = _mm_mul_ps(xval,xweight1);

				_mm_store_ps(weight,xweight1);
				_mm_store_ps(factor,xfactor1);

				float weight_sum = weight[0]+weight[1]+weight[2]+weight[3];
				float factor_sum = factor[0]+factor[1]+factor[2]+factor[3];

				if (weight_sum>0) 
				{
					float d = factor_sum/weight_sum;
					if (d>=0) D_tmp.at(v,u-1) = d;
				}
			}
		}

		// vertical filter
		for (unsigned long u = 3; u < D.column-3; u++) 
		{
			// init
			for (int v=0; v<3; v++)
				val[v] = D_tmp.at(v,u);

			// loop
			for (unsigned long v = 3; v < D.row; v++)
			{
				// set
				float val_curr = D_tmp.at(v-1,u);
				val[v%4] = D_tmp.at(v,u);

				xval     = _mm_load_ps(val);      
				xweight1 = _mm_sub_ps(xval,_mm_set1_ps(val_curr));
				xweight1 = _mm_and_ps(xweight1,xabsmask);
				xweight1 = _mm_sub_ps(xconst4,xweight1);
				xweight1 = _mm_max_ps(xconst0,xweight1);
				xfactor1 = _mm_mul_ps(xval,xweight1);

				_mm_store_ps(weight,xweight1);
				_mm_store_ps(factor,xfactor1);

				float weight_sum = weight[0]+weight[1]+weight[2]+weight[3];
				float factor_sum = factor[0]+factor[1]+factor[2]+factor[3];

				if (weight_sum>0)
				{
					float d = factor_sum/weight_sum;
					if (d>=0) D.at(v-1,u) = d;
				}
			}
		}

		// full resolution: 8 pixel bilateral filter width
	} 
	else 
	{
		// horizontal filter
		for (unsigned long v = 3; v < D.row; v++) 
		{
			// init
			for (int u=0; u<7; u++)
				val[u] = D_copy.at(v,u);

			// loop
			for (unsigned long u = 7; u < D.column-3; u++) 
			{
				// set
				float val_curr = D_copy.at(v,u-3);
				val[u%8] = D_copy.at(v,u);

				xval     = _mm_load_ps(val);      
				xweight1 = _mm_sub_ps(xval,_mm_set1_ps(val_curr));
				xweight1 = _mm_and_ps(xweight1,xabsmask);
				xweight1 = _mm_sub_ps(xconst4,xweight1);
				xweight1 = _mm_max_ps(xconst0,xweight1);
				xfactor1 = _mm_mul_ps(xval,xweight1);

				xval     = _mm_load_ps(val+4);      
				xweight2 = _mm_sub_ps(xval,_mm_set1_ps(val_curr));
				xweight2 = _mm_and_ps(xweight2,xabsmask);
				xweight2 = _mm_sub_ps(xconst4,xweight2);
				xweight2 = _mm_max_ps(xconst0,xweight2);
				xfactor2 = _mm_mul_ps(xval,xweight2);

				xweight1 = _mm_add_ps(xweight1,xweight2);
				xfactor1 = _mm_add_ps(xfactor1,xfactor2);

				_mm_store_ps(weight,xweight1);
				_mm_store_ps(factor,xfactor1);

				float weight_sum = weight[0]+weight[1]+weight[2]+weight[3];
				float factor_sum = factor[0]+factor[1]+factor[2]+factor[3];

				if (weight_sum>0) 
				{
					float d = factor_sum/weight_sum;
					if (d>=0) D_tmp.at(v,u-3) = d;
				}
			}
		}

		// vertical filter
		for (unsigned long u = 3; u < D.column-3; u++)
		{
			// init
			for (int v=0; v<7; v++)
				val[v] = D_tmp.at(v,u);

			// loop
			for (unsigned long v = 7; v < D.row; v++)
			{
				// set
				float val_curr = D_tmp.at(v-3,u);
				val[v%8] = D_tmp.at(v,u);

				xval     = _mm_load_ps(val);      
				xweight1 = _mm_sub_ps(xval,_mm_set1_ps(val_curr));
				xweight1 = _mm_and_ps(xweight1,xabsmask);
				xweight1 = _mm_sub_ps(xconst4,xweight1);
				xweight1 = _mm_max_ps(xconst0,xweight1);
				xfactor1 = _mm_mul_ps(xval,xweight1);

				xval     = _mm_load_ps(val+4);      
				xweight2 = _mm_sub_ps(xval,_mm_set1_ps(val_curr));
				xweight2 = _mm_and_ps(xweight2,xabsmask);
				xweight2 = _mm_sub_ps(xconst4,xweight2);
				xweight2 = _mm_max_ps(xconst0,xweight2);
				xfactor2 = _mm_mul_ps(xval,xweight2);

				xweight1 = _mm_add_ps(xweight1,xweight2);
				xfactor1 = _mm_add_ps(xfactor1,xfactor2);

				_mm_store_ps(weight,xweight1);
				_mm_store_ps(factor,xfactor1);

				float weight_sum = weight[0]+weight[1]+weight[2]+weight[3];
				float factor_sum = factor[0]+factor[1]+factor[2]+factor[3];

				if (weight_sum>0) 
				{
					float d = factor_sum/weight_sum;
					if (d>=0) D.at(v-3,u) = d;
				}
			}
		}
	}

	// free memory
	_mm_free(val);
	_mm_free(weight);
	_mm_free(factor);
}


void Elas::median(Matrix<float> &D)
{
	// temporary memory
	Matrix<float> D_temp(D);

	int window_size = 3;

	Matrix<float> vals(1,window_size*2+1);

	// first step: horizontal median filter
	for (ulong u = window_size; u < D.column-window_size; u++) 
	{
		for (ulong v = window_size; v < D.row-window_size; v++) 
		{
			if (D.at(v,u) >= 0) 
			{    
				int j = 0;
				for (ulong u2 = u-window_size; u2 <= u+window_size; u2++) 
				{
					float temp = D.at(v,u2);
					int i = j-1;
					while (i >= 0 && vals.at(0,i) > temp) 
					{
						vals.at(0,i+1) = vals.at(0,i);
						i--;
					}
					vals.at(0,i+1) = temp;
					j++;
				}
				D_temp.at(v,u) = vals.at(0,window_size);
			} 
			else 
			{
				D_temp.at(v,u) = D.at(v,u);
			}

		}
	}

	// second step: vertical median filter
	for (ulong u = window_size; u < D.column-window_size; u++) 
	{
		for (ulong v = window_size; v < D.row-window_size; v++) 
		{
			if (D.at(v,u) >= 0) 
			{
				int j = 0;
				for (ulong v2 = v-window_size; v2 <= v+window_size; v2++) 
				{
					float temp = D_temp.at(v2,u);
					int i = j-1;
					while (i >= 0 && vals.at(0,i) > temp) 
					{
						vals.at(0,i+1) = vals.at(0,i);
						i--;
					}
					vals.at(0,i+1) = temp;
					j++;
				}
				D.at(v,u) = vals.at(0,window_size);
			} else {
				D.at(v,u) = D.at(v,u);
			}
		}
	}
}


void Elas::invalidD(Matrix<float> &D)
{
	for (unsigned long u = 0; u < D.column; u++) 
	{
		for (unsigned long v = 0; v < D.row; v++) 
		{
			for (unsigned long c = 0; c < D.channel; c++)
			{
				if (D.at(v,u,c) < 0) D.at(v,u,c) = 0;
			}
		}
	}
}


void Elas::computeDelaunayTera(Matrix<float> &D)
{
	tetgenio in, out;
	// inputs
	in.numberofpoints = D.column * D.row;
	in.pointlist = new REAL[in.numberofpoints * 3];

	for (unsigned long u = 0; u < D.column; u++) 
	{
		for (unsigned long v = 0; v < D.row; v++) 
		{
			for (unsigned long c = 0; c < D.channel; c++)
			{
				in.pointlist[v * D.column + u] = u;
				in.pointlist[v * D.column + u + 1] = v;
				in.pointlist[v * D.column + u + 2] = D.at(v,u,c);
			}
		}
	}

	// Output the PLC to files 'barin.node' and 'barin.poly'.
	in.save_nodes("temp/dispin");
	in.save_poly("temp/dispin");

	// Tetrahedralize the PLC. Switches are chosen to read a PLC (p),
	// do quality mesh generation (q) with a specified quality bound
	// (1.414), and apply a maximum volume constraint (a0.1).
	tetrahedralize(new tetgenbehavior, &in, &out);
	// Output mesh to files 'barout.node', 'barout.ele' and 'barout.face'.
	out.save_nodes("temp/dispout");
	out.save_elements("temp/dispout");
	out.save_faces("temp/dispout");
	out.save_faces2smesh("temp/dispout");
}


void Elas::removeInconsistentSupportPoints(Matrix<short>& D_can)
{
	// for all valid support points do
	for (int u_can = 0; u_can < (int)D_can.column; u_can++) 
	{
		for (int v_can = 0; v_can < (int)D_can.row; v_can++) 
		{
			short d_can = D_can.at(v_can,u_can);
			if (d_can >= 0) 
			{
				// compute number of other points supporting the current point
				int support = 0;
				for (int u_can_2 = u_can-param.incon_window_size; u_can_2 <= u_can+param.incon_window_size; u_can_2++) 
				{
					for (int v_can_2 = v_can-param.incon_window_size; v_can_2 <= v_can+param.incon_window_size; v_can_2++) 
					{
						if (u_can_2 >= 0 && v_can_2 >= 0 && u_can_2 < (int)D_can.column && v_can_2 < (int)D_can.row) 
						{
							short d_can_2 = D_can.at(v_can_2,u_can_2);
							if (d_can_2 >= 0 && abs(d_can-d_can_2) <= param.incon_threshold)
								support++;
						}
					}
				}

				// invalidate support point if number of supporting points is too low
				if (support < param.incon_min_support)
					D_can.at(v_can,u_can) = -1;
			}
		}
	}
}



void Elas::removeRedundantSupportPoints(Matrix<short>& D_can,int redun_max_dist, int redun_threshold, bool vertical)
{
	// parameters
	int redun_dir_u[2] = {0,0};
	int redun_dir_v[2] = {0,0};
	if (vertical) 
	{
		redun_dir_v[0] = -1;
		redun_dir_v[1] = +1;
	}
	else 
	{
		redun_dir_u[0] = -1;
		redun_dir_u[1] = +1;
	}

	// for all valid support points do
	for (int u_can = 0; u_can < (int)D_can.column; u_can++)
	{
		for (int v_can = 0; v_can < (int)D_can.row; v_can++) 
		{
			short d_can = D_can.at(v_can,u_can);
			if (d_can >= 0) 
			{
				// check all directions for redundancy
				bool redundant = true;
				for (int i=0; i<2; i++)
				{
					// search for support
					int u_can_2 = u_can;
					int v_can_2 = v_can;
					bool support = false;
					for (int j = 0; j < redun_max_dist; j++)
					{
						u_can_2 += redun_dir_u[i];
						v_can_2 += redun_dir_v[i];
						if (u_can_2 < 0 || v_can_2 < 0 || u_can_2 >= (int)D_can.column || v_can_2 >= (int)D_can.row)
							break;
						short d_can_2 = D_can.at(v_can_2,u_can_2);
						if (d_can_2 >= 0 && abs(d_can-d_can_2) <= redun_threshold)
						{
							support = true;
							break;
						}
					}

					// if we have no support => point is not redundant
					if (!support)
					{
						redundant = false;
						break;
					}
				}

				// invalidate support point if it is redundant
				if (redundant)
					D_can.at(v_can,u_can) = -1;
			}
		}
	}
}



void Elas::addCornerSupportPoints(std::vector<Point3D<int>> &p_support)
{
	// list of border points
	std::vector<Point3D<int>> p_border;
	p_border.push_back(Point3D<int>(0,0,0));
	p_border.push_back(Point3D<int>(0,param.height-1,0));
	p_border.push_back(Point3D<int>(param.width-1,0,0));
	p_border.push_back(Point3D<int>(param.width-1,param.height-1,0));

	// find closest d
	for (std::size_t i=0; i<p_border.size(); i++) 
	{
		//int best_dist = std::numeric_limits<int>::max();
		int best_dist = 0x7fffffff;
		for (std::size_t j=0; j<p_support.size(); j++) 
		{
			int du = p_border[i].x-p_support[j].x;
			int dv = p_border[i].y-p_support[j].y;
			int curr_dist = du*du+dv*dv;
			if (curr_dist < best_dist)
			{
				best_dist = curr_dist;
				p_border[i].z = p_support[j].z;
			}
		}
	}

	// for right image
	p_border.push_back(Point3D<int>(p_border[2].x+p_border[2].z,p_border[2].y,p_border[2].z));
	p_border.push_back(Point3D<int>(p_border[3].x+p_border[3].z,p_border[3].y,p_border[3].z));

	// add border points to support points
	for (std::size_t i=0; i<p_border.size(); i++)
		p_support.push_back(p_border[i]);
}


inline short Elas::computeMatchingDisparity(int u, int v,const Matrix<short>& imLS, const Matrix<short>& imRS, bool right_image)
{
	const int u_step      = 2;
	const int v_step      = 2;
	const int window_size = 2;

	// check if we are inside the image region
	if (u >= window_size && u <= ((int)imLS.column-window_size-1) && v >= window_size && v <= ((int)imLS.row-window_size-1))
	{
		// we require at least some texture
		double sum = 0;
		for (unsigned long i = 0; i < imRS.channel; i++)
		{
			if (right_image) sum += abs(imRS.at(v,u,i));
			else sum += abs(imLS.at(v,u,i));
		}
		sum = abs(sum - (int)imRS.channel * 128 * 4);
		if (sum < param.support_texture) return -1;

		// best match
		//int min_1_E = std::numeric_limits<int>::max();
		double min_1_E = std::numeric_limits<int>::max();
		short min_1_d = -1;
		//int min_2_E = std::numeric_limits<int>::max();
		double min_2_E = std::numeric_limits<int>::max();
		short min_2_d = -1;

		// get valid disparity range
		int disp_min_valid = max(param.disp_min,0);
		int disp_max_valid = param.disp_max;
		if (!right_image) disp_max_valid = min(param.disp_max,u-window_size);
		else              disp_max_valid = min(param.disp_max,(int)imLS.column-u-window_size);

		// assume, that we can compute at least 10 disparities for this pixel
		if (disp_max_valid - disp_min_valid < 10)
			return -1;

		// for all disparities do
		for (short d = disp_min_valid; d <= disp_max_valid; d++)
		{
			//int s1 = 1,s2 = 1;
			sum = 0;
			// compute match energy at this disparity
			for (unsigned long c = 0; c < imLS.channel; c++)
			{
				/*for (int h = -window_size; h <= window_size; h++)
				{
					for (int w = -window_size; w <= window_size; w++)
					{
						if (right_image)
						    sum += abs(imRS.at(v + h,u + w,c) - imLS.at(v + h,u + w + d,c));
						else
							sum += abs(imLS.at(v + h,u + w,c) - imRS.at(v + h,u + w - d,c));
					}
				}*/
				if (right_image)
				{
					sum += abs(imRS.at(v - v_step,u - u_step,c) - imLS.at(v - v_step,u - u_step + d,c));
					sum += abs(imRS.at(v - v_step,u + u_step,c) - imLS.at(v - v_step,u + u_step + d,c));
					sum += abs(imRS.at(v + v_step,u - u_step,c) - imLS.at(v + v_step,u - u_step + d,c));
					sum += abs(imRS.at(v + v_step,u + u_step,c) - imLS.at(v + v_step,u + u_step + d,c));
					sum += abs(imRS.at(v,u,c) - imLS.at(v,u + d,c));

					/*s1 += abs(imRS.at(v,u,c));
					s2 += abs(imLS.at(v,u + d,c));*/
				} 
				else
				{
					sum += abs(imLS.at(v - v_step,u - u_step,c) - imRS.at(v - v_step,u - u_step - d,c));
					sum += abs(imLS.at(v - v_step,u + u_step,c) - imRS.at(v - v_step,u + u_step - d,c));
					sum += abs(imLS.at(v + v_step,u - u_step,c) - imRS.at(v + v_step,u - u_step - d,c));
					sum += abs(imLS.at(v + v_step,u + u_step,c) - imRS.at(v + v_step,u + u_step - d,c));
					sum += abs(imLS.at(v,u,c) - imRS.at(v,u - d,c));

					/*s1 += abs(imLS.at(v,u,c));
					s2 += abs(imRS.at(v,u - d,c));*/
				}
			}

			/*sum /= (s1 + s2);*/

			// best + second best match
			if (sum < min_1_E) 
			{
				min_1_E = sum;
				min_1_d = d;
			} 
			else if (sum < min_2_E) 
			{
				min_2_E = sum;
				min_2_d = d;
			}
		}

		// check if best and second best match are available and if matching ratio is sufficient
		if (min_1_d >= 0 && min_2_d >= 0 && min_1_E < param.support_threshold * min_2_E)
			return min_1_d;
		else
			return -1;

	} else
		return -1;
}


void Elas::drawDelaunay(const std::vector<Point3D<int>>& p_support, const std::vector<Triangle>& tri, const Image& img, bool rImg)
{
	Image imgDraw(img);
	Pixel<uchar> color;
	Point2D<int> bPoint,ePoint;

	for (std::size_t i = 0; i < tri.size(); i++) 
	{
		// get triangle corner indices
		int c1 = tri[i].c1;
		int c2 = tri[i].c2;
		int c3 = tri[i].c3;

		color.randColor();

		bPoint.x = rImg ? (p_support[c1].x - p_support[c1].z) : p_support[c1].x;
		bPoint.y = p_support[c1].y;
		ePoint.x = rImg ? (p_support[c2].x - p_support[c2].z) : p_support[c2].x;
		ePoint.y = p_support[c2].y;
		imgDraw.drawLine(bPoint,ePoint,color);

		bPoint.x = rImg ? (p_support[c2].x - p_support[c2].z) : p_support[c2].x;
		bPoint.y = p_support[c2].y;
		ePoint.x = rImg ? (p_support[c3].x - p_support[c3].z) : p_support[c3].x;
		ePoint.y = p_support[c3].y;
		imgDraw.drawLine(bPoint,ePoint,color);

		bPoint.x = rImg ? (p_support[c3].x - p_support[c3].z) : p_support[c3].x;
		bPoint.y = p_support[c3].y;
		ePoint.x = rImg ? (p_support[c1].x - p_support[c1].z) : p_support[c1].x;
		ePoint.y = p_support[c1].y;
		imgDraw.drawLine(bPoint,ePoint,color);

	} 

	if (rImg)
	   imgDraw.save("temp/imgRT.png");
	else
	   imgDraw.save("temp/imgLT.png");
}


void Elas::process(const Image& imgL,const Image& imgR,Matrix<float> &D1,Matrix<float> &D2)
{
	Matrix<short> vSL(imgL.height,imgL.width,imgL.channel);
	Matrix<short> uSL(imgL.height,imgL.width,imgL.channel);
	Filter::sobel3x3(ImageTools::ImageToGrey(imgL),uSL,vSL);
	Matrix<short> imLS(imgL.height,imgL.width,16);
	createDescriptor(uSL,vSL,imLS,imgL);
	//imLS.abs();

	Matrix<short> vSR(imgR.height,imgR.width,imgR.channel);
	Matrix<short> uSR(imgR.height,imgR.width,imgR.channel);
	Filter::sobel3x3(ImageTools::ImageToGrey(imgR),uSR,vSR);
	Matrix<short> imRS(imgR.height,imgR.width,16);
	createDescriptor(uSR,vSR,imRS,imgR);
	//imRS.abs();

	std::vector<Point3D<int>> p_support = computeSupportMatches(imLS,imRS);

	// if not enough support points for triangulation
	if (p_support.size() < 3) 
	{
		cout << "ERROR: Need at least 3 support points!" << endl;
		return;
	}

	Image matchImg(imgL.width*2,imgL.height,imgL.imgType,3);

	for (unsigned long v = 0; v < matchImg.height; v++)
	{
		for (unsigned long u = 0; u < matchImg.width; u++)
		{
			Pixel<uchar> p;

			if (u >= imgL.width)
				p = imgR.getPixel<uchar>(v,u%imgL.width);
			else
				p = imgL.getPixel<uchar>(v,u%imgL.width);

			matchImg.setPixel(v,u,p);
		}
	}

	for (std::size_t s = 0; s < p_support.size(); s++)
	{
		if(p_support.at(s).z >0 && p_support.at(s).y < (int)(imgL.height - 1) && p_support.at(s).x < (int)(imgL.width - 1) && p_support.at(s).y > 1 && p_support.at(s).x > 1)
		{
			for (int i = -1; i <= 1; i++)
			{
				for (int j = -1; j <= 1; j++)
				{
					matchImg.setPixel(p_support.at(s).y + i,p_support.at(s).x + j,Pixel<uchar>(255,0,0));
					matchImg.setPixel(p_support.at(s).y + i,p_support.at(s).x + imgL.width - p_support.at(s).z + j,Pixel<uchar>(255,0,0));
				}
			}
		}
	}

	matchImg.save("temp/matchImg.png");

	//将所有support points 作为顶点构成德劳内三角
	std::vector<Triangle> tri_1 = computeDelaunayTriangulation(p_support,false);
	std::vector<Triangle> tri_2 = computeDelaunayTriangulation(p_support,true);

	drawDelaunay(p_support,tri_1,imgL,false);
	drawDelaunay(p_support,tri_2,imgR,true);


	//将德劳内三角的三个顶点和顶点的视差值构成三维点
	//三个顶点确定的平面
	/*computeDisparityPlanes(p_support,tri_1,false);
	computeDisparityPlanes(p_support,tri_2,true);*/
	planeFitting(p_support,tri_1);
	planeFitting(p_support,tri_2);

	// allocate memory for disparity grid
	int gridW   = (int)ceil((float)imgL.width/(float)param.grid_size);
	int gridH = (int)ceil((float)imgL.height/(float)param.grid_size);
	Matrix<int> dispGrid1(gridH,gridW,param.disp_max+2);
	Matrix<int> dispGrid2(gridH,gridW,param.disp_max+2);

	createGrid(p_support,dispGrid1,false);
	createGrid(p_support,dispGrid2,true);

	/*std::ofstream dispGrid1File ("temp/dispGrid1Matrix.txt");
	if (dispGrid1File.is_open()) 
	{
		dispGrid1File << dispGrid1;
		dispGrid1File.close();
	}

	std::ofstream dispGrid2File ("temp/dispGrid2Matrix.txt");
	if (dispGrid2File.is_open()) 
	{
		dispGrid2File << dispGrid2;
		dispGrid2File.close();
	}*/

	/*computeDisparity(p_support,tri_1,dispGrid1,imLS,imRS,false,D1);
	computeDisparity(p_support,tri_2,dispGrid2,imLS,imRS,true,D2);*/

	computePlaneDisparity(p_support,tri_1,dispGrid1,imLS,imRS,false,D1);
	computePlaneDisparity(p_support,tri_2,dispGrid2,imLS,imRS,true,D2);

	leftRightConsistencyCheck(D1,D2);

	//将区域内一致（similarity）的像素点的数量过少的区域的视差设置为invalid（-10）
	/*removeSmallSegments(D1);
	if (!param.postprocess_only_left)
	removeSmallSegments(D2);*/

	//对invalid的区域做填补
	/*gapInterpolation(D1);
	if (!param.postprocess_only_left)
		gapInterpolation(D2);

	if (param.filter_adaptive_mean) {
		adaptiveMean(D1);
		if (!param.postprocess_only_left)
			adaptiveMean(D2);
	}

	if (param.filter_median) {
		median(D1);
		if (!param.postprocess_only_left)
			median(D2);
	}*/

	invalidD(D1);
	invalidD(D2);

	std::ofstream d1File ("temp/d1Matrix.txt");
	if (d1File.is_open()) 
	{
		d1File << D1;
		d1File.close();
	}

	std::ofstream d2File ("temp/d2Matrix.txt");
	if (d2File.is_open()) 
	{
		d2File << D2;
		d2File.close();
	}
}


int ElasStereoMatch(const std::string& imL, const std::string& imR)
{
	Image imgL(imL);
	Image imgR(imR);

	// check for correct size
	if (imgL.width <= 0 || imgL.height <= 0 || imgR.width <= 0 || imgR.height <= 0 ||
		imgL.width != imgR.width || imgL.height != imgR.height) {
		cout << "ERROR: Images must be of same size, but" << endl;
		cout << "       I1: " << imgL.width << " x " << imgL.height <<
			", I2: " << imgR.width << " x " << imgR.height << endl;
		return -1;
	}

	// get image width and height
	int width = imgL.width;
	int height = imgL.height;

	//disparity images
	Matrix<float> D1_data, D2_data;

	// process
	Elas::parameters param(Elas::MIDDLEBURY);
	param.width = width;
	param.height = height;
	param.postprocess_only_left = true;
	param.speckle_size = 50;
	param.grid_size = 3;
	param.candidate_stepsize = 1;
	param.lr_threshold = 2;
	param.gamma = 15;
	param.disp_max = 59;
	/*param.filter_adaptive_mean = true;
	param.filter_median = true;*/

	if (param.subsampling)
	{
		D1_data = Matrix<float>(height / 2, width / 2);
		D2_data = Matrix<float>(height / 2, width / 2);
		param.width = width / 2;
		param.height = height / 2;
	}
	else
	{
		D1_data = Matrix<float>(height, width);
		D2_data = Matrix<float>(height, width);
	}

	Elas elas(param);
	elas.process(imgL, imgR, D1_data, D2_data);

	MeshIO<float>::exportPLY(D1_data, imgL, 0.f, "temp/elasL.ply");
	MeshIO<float>::exportPLY(D2_data, imgR, 0.f, "temp/elasR.ply");

	// find maximum disparity for scaling output disparity images to [0..255]
	//D1_data.norm(0.f, 255.f);
	//D2_data.norm(0.f, 255.f);

	Image dispL = ImageExpandTools::Mat2Image(D1_data,4);
	Image dispR = ImageExpandTools::Mat2Image(D2_data,4);

	// save disparity images

	dispL.save("temp/cones_dispL.png");
	dispR.save("temp/cones_dispR.png");
	return 1;
}