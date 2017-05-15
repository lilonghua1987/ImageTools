#include "Delaunay3D.h"

#include "OList/listobj.c"
#include "OList/listhash.c"
#include "OList/listscan.c"
#include "OList/list.c"
#include "OList/chronos.c"
#include "OList/error.c"

#include "geometry.c"
#include "ggveclib.c"
#include "stat.c"
#include "unifgrid.c"
#include "file.c"

StatInfo SI;			/* Statistic Infomations */

/************** Program Flags **************/

boolean CheckFlag	= ON;	/* Whether checking each built tetrahedron */
/* is a Delaunay one.			   */

boolean StatFlag	= OFF;	/* Whether printing Statistic Infomations. */

boolean NumStatFlag	= OFF;	/* Whether printing only numerical values  */
/* of Statistic Infomations.		   */

boolean UGSizeFlag	= OFF;	/* Whether UG size is user defined	   */
int	UGSize		= 1;	/* The UG size user proposed		   */


Delaunay3D::Delaunay3D(const std::vector<Point3D<int>>& pointList)
	:safeTetraFlag(ON)
	,safeFaceFlag(ON)
{
	List   T=NULL_LIST;
	int n = pointList.size();
	Point3 *v = (Point3*)malloc(n*sizeof(Point3));

	for (int i = 0; i < n; i++)
	{
		v[i].x = pointList.at(i).x;
		v[i].y = pointList.at(i).y;
		v[i].z = pointList.at(i).z;
	}
	

	double sec;

	SetProgramName("InCoDe");

	SI.Point=n;

	ResetChronos(USER_CHRONOS);
	StartChronos(USER_CHRONOS);
	T=InCoDe(v,n);
	free(v);
	StopChronos(USER_CHRONOS);

	SI.Tetra=CountList(T);
	tetraMat = Matrix<int>(T->nobject,T->objectsize);

	sec=ReadChronos(USER_CHRONOS);

	SI.Secs=sec;

	if(StatFlag && !NumStatFlag) PrintStat();
	if(StatFlag && NumStatFlag) PrintNumStat();
	if(!StatFlag)
		printf("Points:%7i Secs:%6.2f Tetras:%7i\n",SI.Point,SI.Secs,SI.Tetra);

	ShortTetra *t;
	int i = 0;
	while(ExtractList(&t,T))
	{
		for (int j = 0; j < (int)tetraMat.channel; j++)
		{
			tetraMat.at(i++,j) = t->v[j];
		}
	}

	//free(T);
}


Delaunay3D::~Delaunay3D(void)
{
}


/***************************************************************************
*									   *
* MakeTetra								   *
*									   *
* Given a face find the dd-nearest point to it and joining it to the face  *
* build a new Delaunay tetrahedron.					   *
*									   *
***************************************************************************/
Tetra* Delaunay3D::MakeTetra(Face *f,Point3 *v,int n)
{
	Plane p,Mp;
	boolean found=FALSE;
	double Radius=BIGNUMBER,  rad;
	Line Lc;
	int pind=0;
	Tetra *t;
	Point3 Center, c;
	int i;

	if(!CalcPlane(&(v[f->v[0]]),&(v[f->v[1]]),&(v[f->v[2]]),&p))
		Error("MakeTetra, Face with collinar vertices!\n",EXIT);

	CalcLineofCenter(&(v[f->v[0]]),&(v[f->v[1]]),&(v[f->v[2]]),&Lc);

	for(i=0;i<n;i++)
	{
		if((i!=f->v[0]) &&
			(i!=f->v[1]) &&
			(i!=f->v[2]) &&
			RightSide(&p,&(v[i])) )
		{
			CalcMiddlePlane(&(v[i]),&(v[f->v[0]]),&Mp);
			if(CalcLinePlaneInter(&Lc,&Mp,&c))
			{
				rad=V3SquaredDistanceBetween2Points(&c, &(v[i]));

				if(!RightSide(&p,&c)) rad=-rad;
				if(rad==Radius) Error("MakeTetra, Five cocircular points!\n",EXIT);
				if(rad<Radius)
				{
					found=TRUE;
					Radius=rad;
					Center=c;
					pind=i;
				}
			}
		}
	}

	if(!found) return NULL;

	t=BuildTetra(f,pind);

	if(CheckFlag)  CheckTetra(t,v,n);

	return t;
}


/***************************************************************************
*									   *
* FirstTetra								   *
*									   *
* Build the first Delaunay tetrahedron of the triangulation.		   *
*									   *
***************************************************************************/
Tetra* Delaunay3D::FirstTetra(Point3 *v, int n)
{
	int i, MinIndex=0;
	double Radius, MinRadius=BIGNUMBER;

	Tetra *t;
	Plane p[3], Mp;
	Line Lc;
	Face f;
	Point3 c;

	boolean found=FALSE;

	f.v[0]=0;		/* The first point of the face is the 0 */
	/* (any point is equally good).		*/

	for(i=0;i<n;i++)	/* The 2nd point of the face is the	*/
		if(i!=f.v[0]) 	/* euclidean nearest to first point	*/
		{
			Radius=V3SquaredDistanceBetween2Points(&(v[0]), &(v[i]));
			if(Radius<MinRadius)
			{
				MinRadius=Radius;
				MinIndex=i;
			}
		}

		f.v[1]=MinIndex;
		/* The 3rd point is that with previous	*/
		/* ones builds the smallest circle.	*/

		CalcMiddlePlane(&(v[f.v[0]]),&(v[f.v[1]]), &(p[0]));

		MinRadius=BIGNUMBER;

		for(i=0;i<n;i++)
			if(i!=f.v[0] && i!=f.v[1])
			{
				CalcMiddlePlane(&(v[f.v[0]]), &(v[i]),&(p[1]));
				if(CalcPlane(&(v[f.v[0]]),&(v[f.v[1]]),&(v[i]),&(p[2])))
					if(CalcPlaneInter(&(p[0]), &(p[1]), &(p[2]), &c))
					{
						Radius=V3DistanceBetween2Points(&c, &(v[0]));
						if(Radius<MinRadius)
						{
							MinRadius=Radius;
							MinIndex=i;
						}
					}
			}

			f.v[2]=MinIndex;


			/* The first tetrahedron construction is analogous to normal */
			/* MakeTetra, only we dont limit search to an halfspace.     */

			MinRadius=BIGNUMBER;

			CalcPlane(&(v[f.v[0]]),&(v[f.v[1]]),&(v[f.v[2]]),&p[0]);
			CalcLineofCenter(&(v[f.v[0]]),&(v[f.v[1]]),&(v[f.v[2]]),&Lc);

			for(i=0;i<n;i++)
				if(i!=f.v[0] && i!=f.v[1] && i!=f.v[2] )
				{
					CalcMiddlePlane(&(v[i]),&(v[f.v[0]]),&Mp);
					if(CalcLinePlaneInter(&Lc,&Mp,&c))
					{
						Radius=V3SquaredDistanceBetween2Points(&c, &(v[i]));

						if(MinRadius==Radius) Error("FirstTetra, Five cocircular points!\n",EXIT);
						if(Radius<MinRadius)
						{
							found=TRUE;
							MinRadius=Radius;
							MinIndex=i;
						}
					}
				}

				if(!found) Error("FirstTetra, Planar dataset, unable to build first tetrahedron.\n",EXIT);

				if(!RightSide(&p[0],&(v[MinIndex])))	ReverseFace(&f);

				t=BuildTetra(&f, MinIndex);

				CheckTetra(t,v,n);

				ReverseFace(t->f[0]); /* First Face in first Tetra   */
				/* must be outward oriented    */
				return t;
}



/***************************************************************************
*									   *
* InCoDe								   *
*									   *
* Given a vector v of n Point3 this function returns the tetrahedra list   *
* of the Delaunay triangulation of points. This Functions uses the incre-  *
* mental construction algoritm InCoDe [Cignoni 92].			   *
*									   *
* We build the triangulation adding to an initial Delaunay tetrahedron new *
* Delaunay tetrahedra built from its faces (using MakeTetra function).	   *
*									   *
* This algorithm make use of Speed up techniques suggeste in the paper	   *
* yelding an average linear performance (against a teorethical cubic worst *
* case.									   *
*									   *
* [Cignoni 92]								   *
* P. Cignoni, C. Montani, R. Scopigno		  d			   *
*"A Merge-First Divide and Conquer Algorithm for E  Delaunay Triangulation"*
* CNUCE Internal Report C92/16 Oct 1992					   *
*									   *
***************************************************************************/
List Delaunay3D::InCoDe(Point3 *v, int n)
{
	List Q=NULL_LIST;
	List T=NULL_LIST;
	List OldFace=NULL_LIST;

	Tetra *t;
	ShortTetra *st;
	Face  *f;
	UG g;

	Q=NewList(FIFO,sizeof(Face));			/* Initialize Active Face */
	ChangeEqualObjectList(EqualFace,Q);		/* List Q.		  */
	HashList(n/4,HashFace,Q);

	if(safeFaceFlag)
	{						/* Initialize list for	  */
		OldFace=NewList(FIFO,sizeof(Face)); 	/* preventing double face */
		ChangeEqualObjectList(EqualFace,OldFace);	/* looping on numerical   */
		HashList(n/4,HashFace,OldFace);		/* errors		  */
	}

	T=NewList(FIFO,sizeof(ShortTetra));		/* Initialize Built Tetra-*/
	ChangeEqualObjectList(EqualTetra,T);		/* hedra List T.	  */
	if(safeTetraFlag) HashList(n/4,HashTetra,T);

	if(UGSizeFlag) BuildUG(v,n,UGSize,&g); 	/* Initialize Uniform Grid */
	else BuildUG(v,n,n,&g);

	t=FirstTetra(v,n);

	for(int i = 0; i < 4; i++)
	{
		InsertList(t->f[i],Q);
		if(safeFaceFlag)  InsertList(t->f[i],OldFace);
		for(int j=0;j<3;j++)
			if(g.UsedPoint[t->f[i]->v[j]]==-1) 
				g.UsedPoint[t->f[i]->v[j]]=1;
			else 
				g.UsedPoint[t->f[i]->v[j]]++; 
	}

	st=Tetra2ShortTetra(t);

	free(t);

	InsertList(st,T);

	while(ExtractList(&f,Q))
	{
		t = FastMakeTetra(f,v,n,&g);

		if(t==NULL) SI.CHFace++;
		else
		{
			st=Tetra2ShortTetra(t);

			if(safeTetraFlag) if(MemberList(st, T))
				Error("Cyclic Tetrahedra Creation\n",EXIT);
			InsertList(st,T);

			for(int i = 1; i < 4; i++)
				if(MemberList(t->f[i],Q))
				{
					DeleteCurrList(Q);
					//free(t->f[i]);

					for(int j = 0; j < 3; j++)
						g.UsedPoint[t->f[i]->v[j]]--;
				}
				else
				{
					InsertList(t->f[i],Q);

					for(int j = 0; j < 3; j++)
						if(g.UsedPoint[t->f[i]->v[j]]==-1) 
							g.UsedPoint[t->f[i]->v[j]]=1;
						else g.UsedPoint[t->f[i]->v[j]]++;

						if(safeFaceFlag)
						{
							if(MemberList(t->f[i], OldFace))
								Error("Cyclic insertion in Active Face List\n",EXIT);
							InsertList(t->f[i],OldFace);
						}
				}
				free(t->f[0]);
				free(t);
		}
		for(int i = 0; i < 3; i++)                   
		{
			g.UsedPoint[f->v[i]]--;
			/*  if(g.UsedPoint[f->v[i]]==0) printf("Point %i completed\n",f->v[i]);*/
		}
		if(!safeFaceFlag) free(f);
	}


	return T;
}