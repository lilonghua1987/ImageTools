// SLIC.cpp: implementation of the SLIC class.
//===========================================================================
// This code implements the zero parameter superpixel segmentation technique
// described in:
//
//
//
// "SLIC Superpixels Compared to State-of-the-art Superpixel Methods"
//
// Radhakrishna Achanta, Appu Shaji, Kevin Smith, Aurelien Lucchi, Pascal Fua,
// and Sabine Susstrunk,
//
// IEEE TPAMI, Volume 34, Issue 11, Pages 2274-2282, November 2012.
//
//
//===========================================================================
// Copyright (c) 2013 Radhakrishna Achanta.
//
// For commercial use please contact the author:
//
// Email: firstname.lastname@epfl.ch
//===========================================================================

#include "SLIC.h"
#include <cfloat>
#include <cmath>
#include <iostream>
#include <fstream>

// For superpixels
const int dx4[4] = {-1,  0,  1,  0};
const int dy4[4] = { 0, -1,  0,  1};
const int dx8[8] = {-1, -1,  0,  1, 1, 1, 0, -1};
const int dy8[8] = { 0, -1, -1, -1, 0, 1, 1,  1};

// For supervoxels
const int dx10[10] = {-1,  0,  1,  0, -1,  1,  1, -1,  0, 0};
const int dy10[10] = { 0, -1,  0,  1, -1, -1,  1,  1,  0, 0};
const int dz10[10] = { 0,  0,  0,  0,  0,  0,  0,  0, -1, 1};

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

SLIC::SLIC()
{
}

SLIC::~SLIC()
{
}

//==============================================================================
///	RGB2XYZ
///

//=================================================================================
/// DrawContoursAroundSegments
///
/// Internal contour drawing option exists. One only needs to comment the if
/// statement inside the loop that looks at neighbourhood.
//=================================================================================
void SLIC::DrawContoursAroundSegments(Image& img, const std::vector<int>& labels)
{
	const int dx8[8] = {-1, -1,  0,  1, 1, 1, 0, -1};
	const int dy8[8] = { 0, -1, -1, -1, 0, 1, 1,  1};

	int sz = labels.size();

	vector<bool> istaken(sz, false);

	int mainindex(0);
	for (ulong j = 0; j < img.height; j++)
	{
		for (ulong k = 0; k < img.width; k++)
		{
			int np(0);
			for( int i = 0; i < 8; i++ )
			{
				int x = k + dx8[i];
				int y = j + dy8[i];

				if ((x >= 0 && x < (int)img.width) && (y >= 0 && y < (int)img.height))
				{
					int index = y*(int)img.width + x;

					if( false == istaken[index] )//comment this to obtain internal contours
					{
						if( labels[mainindex] != labels[index] ) np++;
					}
				}
			}
			if( np > 1 )//change to 2 or 3 for thinner lines
			{
				img.setPixel(j, k, Pixel<uchar>());
				istaken[mainindex] = true;
			}
			mainindex++;
		}
	}
}

//=================================================================================
/// DrawContoursAroundSegmentsTwoColors
///
/// Internal contour drawing option exists. One only needs to comment the if
/// statement inside the loop that looks at neighbourhood.
//=================================================================================
void SLIC::DrawContoursAroundSegmentsTwoColors(Image& img,const std::vector<int>& labels)
{
	const int dx[8] = {-1, -1,  0,  1, 1, 1, 0, -1};
	const int dy[8] = { 0, -1, -1, -1, 0, 1, 1,  1};

	int sz = labels.size();

	vector<bool> istaken(sz, false);

	vector<int> contourx(sz);
	vector<int> contoury(sz);
	int mainindex(0);
	int cind(0);
	for( ulong j = 0; j < img.height; j++ )
	{
		for( ulong k = 0; k < img.width; k++ )
		{
			int np(0);
			for( ulong i = 0; i < 8; i++ )
			{
				int x = k + dx8[i];
				int y = j + dy8[i];

				if ((x >= 0 && x < (int)img.width) && (y >= 0 && y < (int)img.height))
				{
					int index = y*(int)img.width + x;

					//if( false == istaken[index] )//comment this to obtain internal contours
					{
						if( labels[mainindex] != labels[index] ) np++;
					}
				}
			}
			if( np > 1 )
			{
				contourx[cind] = k;
				contoury[cind] = j;
				istaken[mainindex] = true;
				//img[mainindex] = color;
				cind++;
			}
			mainindex++;
		}
	}

	int numboundpix = cind;//int(contourx.size());

	for( int j = 0; j < numboundpix; j++ )
	{
		int ii = contoury[j]*img.width + contourx[j];
		img.setPixel(contoury[j],contourx[j],Pixel<uchar>(255,255,255));
		//----------------------------------
		// Uncomment this for thicker lines
		//----------------------------------
		for( int n = 0; n < 8; n++ )
		{
			ulong x = contourx[j] + dx[n];
			ulong y = contoury[j] + dy[n];
			if( (x >= 0 && x < img.width) && (y >= 0 && y < img.height) )
			{
				ulong ind = y*img.width + x;
				if(!istaken[ind]) img.setPixel(y,x,Pixel<uchar>());
			}
		}
	}
}


//==============================================================================
///	DetectLabEdges
//==============================================================================
void SLIC::DetectLabEdges(const Matrix<float>& lab,vector<double>& edges)
{
	for( int j = 1; j < (int)lab.row-1; j++ )
	{
		for( int k = 1; k < (int)lab.column-1; k++ )
		{
			double dx = 0,dy = 0;
			for (ulong c = 0; c < lab.channel; c++)
			{
				dx += pow((lab.at(j,k-1,c) - lab.at(j,k+1,c)),2.0f);
				dy += pow((lab.at(j-1,k,c) - lab.at(j+1,k,c)),2.0f);
			}

			//edges[i] = (sqrt(dx) + sqrt(dy));
			edges[j*lab.column+k] = (dx + dy);
		}
	}
}

//===========================================================================
///	PerturbSeeds
//===========================================================================
void SLIC::PerturbSeeds(
	vector<double>&				kseedsl,
	vector<double>&				kseedsa,
	vector<double>&				kseedsb,
	vector<double>&				kseedsx,
	vector<double>&				kseedsy,
	const vector<double>&		edges)
{	
	int numseeds = kseedsl.size();

	for( int n = 0; n < numseeds; n++ )
	{
		int ox = static_cast<int>(kseedsx[n]);//original x
		int oy = static_cast<int>(kseedsy[n]);//original y
		int oind = oy*(int)lab.column + ox;

		int storeind = oind;
		for( int i = 0; i < 10; i++ )
		{
			int nx = ox+dx10[i];//new x
			int ny = oy+dy10[i];//new y

			if( nx >= 0 && nx < (int)lab.column && ny >= 0 && ny < (int)lab.row)
			{
				int nind = ny*(int)lab.column + nx;
				if( edges[nind] < edges[storeind])
				{
					storeind = nind;
				}
			}
		}
		if(storeind != oind)
		{
			kseedsx[n] = storeind%lab.column;
			kseedsy[n] = storeind/lab.column;
			kseedsl[n] = lab.at((ulong)kseedsy[n],(ulong)kseedsx[n],0);
			kseedsa[n] = lab.at((ulong)kseedsy[n],(ulong)kseedsx[n],1);
			kseedsb[n] = lab.at((ulong)kseedsy[n],(ulong)kseedsx[n],2);
		}
	}
}


//===========================================================================
///	GetLABXYSeeds_ForGivenStepSize
///
/// The k seed values are taken as uniform spatial pixel samples.
//===========================================================================
void SLIC::GetLABXYSeeds_ForGivenStepSize(
	vector<double>&				kseedsl,
	vector<double>&				kseedsa,
	vector<double>&				kseedsb,
	vector<double>&				kseedsx,
	vector<double>&				kseedsy,
	const int&					STEP,
	const bool&					perturbseeds,
	const vector<double>&		edgemag)
{
	int numseeds(0);
	int n(0);

	//int xstrips = m_width/STEP;
	//int ystrips = m_height/STEP;
	int xstrips = static_cast<int>(0.5+double(m_width)/double(STEP));
	int ystrips = static_cast<int>(0.5+double(m_height)/double(STEP));

	int xerr = m_width  - STEP*xstrips;
	int yerr = m_height - STEP*ystrips;

	double xerrperstrip = double(xerr)/double(xstrips);
	double yerrperstrip = double(yerr)/double(ystrips);

	int xoff = STEP/2;
	int yoff = STEP/2;
	//-------------------------
	numseeds = xstrips*ystrips;
	//-------------------------
	kseedsl.resize(numseeds);
	kseedsa.resize(numseeds);
	kseedsb.resize(numseeds);
	kseedsx.resize(numseeds);
	kseedsy.resize(numseeds);

	for( int y = 0; y < ystrips; y++ )
	{
		int ye = static_cast<int>(y*yerrperstrip);
		for( int x = 0; x < xstrips; x++ )
		{
			int xe = static_cast<int>(x*xerrperstrip);
			int i = (y*STEP+yoff+ye)*m_width + (x*STEP+xoff+xe);
			
			kseedsl.push_back(lab.at((y*STEP+yoff+ye),(x*STEP+xoff+xe),0));
			kseedsa.push_back(lab.at((y*STEP+yoff+ye),(x*STEP+xoff+xe),1));
			kseedsb.push_back(lab.at((y*STEP+yoff+ye),(x*STEP+xoff+xe),2));
			kseedsx[n] = (x*STEP+xoff+xe);
			kseedsy[n] = (y*STEP+yoff+ye);
			n++;
		}
	}

	
	if(perturbseeds)
	{
		PerturbSeeds(kseedsl, kseedsa, kseedsb, kseedsx, kseedsy, edgemag);
	}
}

//===========================================================================
///	GetLABXYSeeds_ForGivenK
///
/// The k seed values are taken as uniform spatial pixel samples.
//===========================================================================
void SLIC::GetLABXYSeeds_ForGivenK(
	vector<double>&				kseedsl,
	vector<double>&				kseedsa,
	vector<double>&				kseedsb,
	vector<double>&				kseedsx,
	vector<double>&				kseedsy,
	const int&					K,
	const bool&					perturbseeds,
	const vector<double>&		edgemag)
{
	int sz = (int)(lab.row*lab.column);
	double step = sqrt(double(sz)/double(K));
	int xoff = static_cast<int>(step/2);
	int yoff = static_cast<int>(step/2);
	
	int n(0);int r(0);
	for( int y = 0; y < m_height; y++ )
	{
		int Y = static_cast<int>(y*step + yoff);
		if( Y > m_height-1 ) break;

		for( int x = 0; x < m_width; x++ )
		{
			//int X = x*step + xoff;//square grid
			int X = static_cast<int>(x*step + (xoff<<(r&0x1)));//hex grid
			if(X > m_width-1) break;
			
			kseedsl.push_back(lab.at(Y,X,0));
			kseedsa.push_back(lab.at(Y,X,1));
			kseedsb.push_back(lab.at(Y,X,2));
			kseedsx.push_back(X);
			kseedsy.push_back(Y);
			n++;
		}
		r++;
	}

	if(perturbseeds)
	{
		PerturbSeeds(kseedsl, kseedsa, kseedsb, kseedsx, kseedsy, edgemag);
	}
}


//===========================================================================
///	PerformSuperpixelSegmentation_VariableSandM
///
///	Magic SLIC - no parameters
///
///	Performs k mean segmentation. It is fast because it looks locally, not
/// over the entire image.
/// This function picks the maximum value of color distance as compact factor
/// M and maximum pixel distance as grid step size S from each cluster (13 April 2011).
/// So no need to input a constant value of M and S. There are two clear
/// advantages:
///
/// [1] The algorithm now better handles both textured and non-textured regions
/// [2] There is not need to set any parameters!!!
///
/// SLICO (or SLIC Zero) dynamically varies only the compactness factor S,
/// not the step size S.
//===========================================================================
void SLIC::PerformSuperpixelSegmentation_VariableSandM(
	vector<double>&				kseedsl,
	vector<double>&				kseedsa,
	vector<double>&				kseedsb,
	vector<double>&				kseedsx,
	vector<double>&				kseedsy,
	std::vector<int>&		    klabels,
	const int&					STEP,
	const int&                  m,
	const int&					NUMITR)
{
	int sz = m_width*m_height;
	const int numk = kseedsl.size();
	//double cumerr(99999.9);
	int numitr(0);

	//----------------
	int offset = STEP;
	if(STEP < 10) offset = static_cast<int>(STEP*1.5);
	//----------------

	vector<double> sigmal(numk, 0);
	vector<double> sigmaa(numk, 0);
	vector<double> sigmab(numk, 0);
	vector<double> sigmax(numk, 0);
	vector<double> sigmay(numk, 0);
	vector<int> clustersize(numk, 0);
	vector<double> inv(numk, 0);//to store 1/clustersize[k] values
	vector<double> distxy(sz, DBL_MAX);
	vector<double> distlab(sz, DBL_MAX);
	vector<double> distvec(sz, DBL_MAX);
	vector<double> maxlab(numk, m*m);//THIS IS THE VARIABLE VALUE OF M, just start with 10
	vector<double> maxxy(numk, STEP*STEP);//THIS IS THE VARIABLE VALUE OF M, just start with 10

	double invxywt = 1.0/(STEP*STEP);//NOTE: this is different from how usual SLIC/LKM works

	while( numitr < NUMITR )
	{
		numitr++;

		distvec.assign(sz, DBL_MAX);
		for( int n = 0; n < numk; n++ )
		{
			int y1 = static_cast<int>(max(0.0,			kseedsy[n]-offset));
			int y2 = static_cast<int>(min(m_height*1.0,	kseedsy[n]+offset));
			int x1 = static_cast<int>(max(0.0,			kseedsx[n]-offset));
			int x2 = static_cast<int>(min(m_width*1.0,	kseedsx[n]+offset));

			for( int y = y1; y < y2; y++ )
			{
				for( int x = x1; x < x2; x++ )
				{
					int i = y*m_width + x;
					assert( y < m_height && x < m_width && y >= 0 && x >= 0 );

					double l = lab.at(y,x,0);
					double a = lab.at(y,x,1);
					double b = lab.at(y,x,2);

					distlab[i] =	(l - kseedsl[n])*(l - kseedsl[n]) +
									(a - kseedsa[n])*(a - kseedsa[n]) +
									(b - kseedsb[n])*(b - kseedsb[n]);

					distxy[i] =		(x - kseedsx[n])*(x - kseedsx[n]) +
									(y - kseedsy[n])*(y - kseedsy[n]);

					//------------------------------------------------------------------------
					double dist = distlab[i]/maxlab[n] + distxy[i]*invxywt;//only varying m, prettier superpixels
					//double dist = distlab[i]/maxlab[n] + distxy[i]/maxxy[n];//varying both m and S
					//------------------------------------------------------------------------
					
					if( dist < distvec[i] )
					{
						distvec[i] = dist;
						klabels[i]  = n;
					}
				}
			}
		}
		//-----------------------------------------------------------------
		// Assign the max color distance for a cluster
		//-----------------------------------------------------------------
		if(0 == numitr)
		{
			maxlab.assign(numk,1);
			maxxy.assign(numk,1);
		}
		for( int i = 0; i < sz; i++ )
		{
			if(maxlab[klabels[i]] < distlab[i]) maxlab[klabels[i]] = distlab[i];
			if(maxxy[klabels[i]] < distxy[i]) maxxy[klabels[i]] = distxy[i];
		}
		//-----------------------------------------------------------------
		// Recalculate the centroid and store in the seed values
		//-----------------------------------------------------------------
		sigmal.assign(numk, 0);
		sigmaa.assign(numk, 0);
		sigmab.assign(numk, 0);
		sigmax.assign(numk, 0);
		sigmay.assign(numk, 0);
		clustersize.assign(numk, 0);

		for( int j = 0; j < sz; j++ )
		{
			int temp = klabels[j];
			assert(klabels[j] >= 0);
			sigmal[klabels[j]] += lab.at(j/lab.column,j%lab.column,0);
			sigmaa[klabels[j]] += lab.at(j/lab.column,j%lab.column,1);
			sigmab[klabels[j]] += lab.at(j/lab.column,j%lab.column,2);
			sigmax[klabels[j]] += (j%lab.column);
			sigmay[klabels[j]] += (j/lab.column);

			clustersize[klabels[j]]++;
		}

		for(int k = 0; k < numk; k++ )
		{
			//_ASSERT(clustersize[k] > 0);
			if( clustersize[k] <= 0 ) clustersize[k] = 1;
			inv[k] = 1.0/double(clustersize[k]);//computing inverse now to multiply, than divide later
		}
		
		for(int k = 0; k < numk; k++ )
		{
			kseedsl[k] = sigmal[k]*inv[k];
			kseedsa[k] = sigmaa[k]*inv[k];
			kseedsb[k] = sigmab[k]*inv[k];
			kseedsx[k] = sigmax[k]*inv[k];
			kseedsy[k] = sigmay[k]*inv[k];
		}
	}
}


//===========================================================================
///	EnforceLabelConnectivity
///
///		1. finding an adjacent label for each new component at the start
///		2. if a certain component is too small, assigning the previously found
///		    adjacent label to this component, and not incrementing the label.
//===========================================================================
void SLIC::EnforceLabelConnectivity(
	const std::vector<int>&		labels,//input labels that need to be corrected to remove stray labels
	const int&					width,
	const int&					height,
	std::vector<int>&			nlabels,//new labels
	int&						numlabels,//the number of labels changes in the end if segments are removed
	const int&					K) //the number of superpixels desired by the user
{

	const int sz = width*height;
	const int SUPSZ = sz/K;
	nlabels.resize(sz, -1);
	int label(0);
	std::vector<int> xvec(sz);
	std::vector<int> yvec(sz);
	int oindex(0);
	int adjlabel(0);//adjacent label
	for( int j = 0; j < height; j++ )
	{
		for( int k = 0; k < width; k++ )
		{
			if( 0 > nlabels[oindex] )
			{
				nlabels[oindex] = label;
				//--------------------
				// Start a new segment
				//--------------------
				xvec[0] = k;
				yvec[0] = j;
				//-------------------------------------------------------
				// Quickly find an adjacent label for use later if needed
				//-------------------------------------------------------
				for( int n = 0; n < 4; n++ )
				{
					int x = xvec[0] + dx4[n];
					int y = yvec[0] + dy4[n];
					if( (x >= 0 && x < width) && (y >= 0 && y < height) )
					{
						int nindex = y*width + x;
						if(nlabels[nindex] >= 0) adjlabel = nlabels[nindex];
					}
				}

				int count(1);
				for( int c = 0; c < count; c++ )
				{
					for( int n = 0; n < 4; n++ )
					{
						int x = xvec[c] + dx4[n];
						int y = yvec[c] + dy4[n];

						if( (x >= 0 && x < width) && (y >= 0 && y < height) )
						{
							int nindex = y*width + x;

							if( 0 > nlabels[nindex] && labels[oindex] == labels[nindex] )
							{
								xvec[count] = x;
								yvec[count] = y;
								nlabels[nindex] = label;
								count++;
							}
						}

					}
				}
				//-------------------------------------------------------
				// If segment size is less then a limit, assign an
				// adjacent label found before, and decrement label count.
				//-------------------------------------------------------
				if(count <= SUPSZ >> 2)
				{
					for( int c = 0; c < count; c++ )
					{
						int ind = yvec[c]*width+xvec[c];
						nlabels[ind] = adjlabel;
					}
					label--;
				}
				label++;
			}
			oindex++;
		}
	}
	numlabels = label;
}

//===========================================================================
///	PerformSLICO_ForGivenStepSize
///
/// There is option to save the labels if needed.
//===========================================================================
void SLIC::PerformSLICO_ForGivenStepSize(
	const Image&			    img,
	std::vector<int>&			klabels,
	int&						numlabels,
	const int&					STEP,
	const double&				m)
{
	vector<double> kseedsl(0);
	vector<double> kseedsa(0);
	vector<double> kseedsb(0);
	vector<double> kseedsx(0);
	vector<double> kseedsy(0);

	//--------------------------------------------------
	m_width  = (int)img.width;
	m_height = (int)img.height;
	int sz = m_width*m_height;
	
	//--------------------------------------------------
	klabels.resize( sz, -1 );
	//--------------------------------------------------
	lab = Matrix<float>(img.height,img.width,img.channel);
	ColorConversion::RGBToLab(img,lab);
	//--------------------------------------------------

	bool perturbseeds(true);
	std::vector<double> edgemag(0);
	if(perturbseeds) DetectLabEdges(lab, edgemag);
	GetLABXYSeeds_ForGivenStepSize(kseedsl, kseedsa, kseedsb, kseedsx, kseedsy, STEP, perturbseeds, edgemag);

	PerformSuperpixelSegmentation_VariableSandM(kseedsl,kseedsa,kseedsb,kseedsx,kseedsy,klabels,STEP,10);
	numlabels = kseedsl.size();

	std::vector<int> nlabels;
	EnforceLabelConnectivity(klabels, m_width, m_height, nlabels, numlabels, static_cast<int>(double(sz)/double(STEP*STEP)));
	klabels.swap(nlabels);
}

//===========================================================================
///	PerformSLICO_ForGivenK
///
/// Zero parameter SLIC algorithm for a given number K of superpixels.
//===========================================================================
void SLIC::PerformSLICO_ForGivenK(
	const Image&			img,
	std::vector<int>&			klabels,
	int&						numlabels,
	const int&					K,//required number of superpixels
	const double&				m)//weight given to spatial distance
{
	vector<double> kseedsl(0);
	vector<double> kseedsa(0);
	vector<double> kseedsb(0);
	vector<double> kseedsx(0);
	vector<double> kseedsy(0);

	//--------------------------------------------------
	m_width  = img.width;
	m_height = img.height;
	int sz = m_width*m_height;
	//--------------------------------------------------
	klabels.resize( sz, -1 );
	//--------------------------------------------------

	lab = Matrix<float>(img.height,img.width,img.channel);
	ColorConversion::RGBToLab(img,lab);
	//--------------------------------------------------

	bool perturbseeds(true);
	std::vector<double> edgemag(img.height*img.width,0);
	if(perturbseeds) DetectLabEdges(lab, edgemag);
	GetLABXYSeeds_ForGivenK(kseedsl, kseedsa, kseedsb, kseedsx, kseedsy, K, perturbseeds, edgemag);

	int STEP = static_cast<int>(sqrt(double(sz)/double(K)) + 2);//adding a small value in the even the STEP size is too small.
	//PerformSuperpixelSLIC(kseedsl, kseedsa, kseedsb, kseedsx, kseedsy, klabels, STEP, edgemag, m);
	PerformSuperpixelSegmentation_VariableSandM(kseedsl,kseedsa,kseedsb,kseedsx,kseedsy,klabels,STEP,(int)m,15);
	numlabels = kseedsl.size();

	std::vector<int> nlabels;
	EnforceLabelConnectivity(klabels, m_width, m_height, nlabels, numlabels, K);
	klabels.swap(nlabels);
}