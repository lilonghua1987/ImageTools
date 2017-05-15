#pragma once

#pragma warning(disable:4996)

extern "C"
{
	#include "OList/general.h"
    #include "OList/error.h"
	#include "OList/olist.h"
	#include "OList/chronos.h"
	#include "graphics.h"
	#include "incode.h"
}

#include <vector>

#include "../CoreStruct.h"
#include "../Matrix.h"

class Delaunay3D
{
public:
	Delaunay3D(const std::vector<Point3D<int>>& pointList);
	~Delaunay3D(void);

private:
	Tetra *MakeTetra(Face *f,Point3 *v,int n);

	Tetra *FirstTetra(Point3 *v, int n);

	List InCoDe(Point3 *v, int n);

private:
	boolean safeFaceFlag;	/* Whether Testing AFL inserted Face to not*/
	/* twice inserted. This situation happens  */
	/* for numerical errors and cause loop the */
	/* whole algorithm.			   */

	boolean safeTetraFlag;	/* Analogous to SafeFaceFlag but check the */
	/* tetrahedra instead of faces (Faster!)   */
public:
	Matrix<int> tetraMat;
};

