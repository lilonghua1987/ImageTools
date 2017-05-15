// SLIC.h: interface for the SLIC class.
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
// For Example:
//  int numlabels(0);
//	int m_spcount = 200,m_compactness = 10.0;//Compactness factor. use a value ranging from 10 to 40 depending on your needs. Default is 10
//	Image img(imgLeft.c_str());
//	int* labels = new int[img.width*img.height];
//	SLIC slic;
//	slic.PerformSLICO_ForGivenK(img, labels, numlabels, m_spcount, m_compactness);//for a given number K of superpixels
//	//slic.PerformSLICO_ForGivenStepSize(img, width, height, labels, numlabels, m_stepsize, m_compactness);//for a given grid step size
//	//slic.DrawContoursAroundSegments(img, labels, width, height, 0);//for black contours around superpixels
//	slic.DrawContoursAroundSegmentsTwoColors(img, labels);//for black-and-white contours around superpixels
//	//slic.SaveSuperpixelLabels(labels,width,height,picvec[k],saveLocation);
//	if(labels) delete [] labels;
//	ImageTools::SaveImage(img,"temp//slic.jpg");
//===========================================================================
#pragma once

#if !defined(_SLIC_H_INCLUDED_)
#define _SLIC_H_INCLUDED_

#define _CRT_SECURE_NO_WARNINGS
#include <vector>
#include <string>
#include <algorithm>
using namespace std;
#include "ColorConversion.h"


class SLIC  
{
public:
	SLIC();
	virtual ~SLIC();
	//============================================================================
	// Superpixel segmentation for a given step size (superpixel size ~= step*step)
	//============================================================================
	void PerformSLICO_ForGivenStepSize(
		const Image&			    img,//Each 32 bit unsigned int contains ARGB pixel values.
		std::vector<int>&			klabels,
		int&						numlabels,
		const int&					STEP,
		const double&				m);
	//============================================================================
	// Superpixel segmentation for a given number of superpixels
	//============================================================================
	void PerformSLICO_ForGivenK(
		const Image&			    img,
		std::vector<int>&			klabels,
		int&						numlabels,
		const int&					K,
		const double&				m);

	//============================================================================
	// Function to draw boundaries around superpixels of a given 'color'.
	// Can also be used to draw boundaries around supervoxels, i.e layer by layer.
	//============================================================================
	void DrawContoursAroundSegments(Image& img, const std::vector<int>& labels);

	void DrawContoursAroundSegmentsTwoColors(Image&	img,const std::vector<int>& labels);

private:

	//============================================================================
	// Magic SLIC. No need to set M (compactness factor) and S (step size).
	// SLICO (SLIC Zero) varies only M dynamicaly, not S.
	//============================================================================
	void PerformSuperpixelSegmentation_VariableSandM(
		vector<double>&				kseedsl,
		vector<double>&				kseedsa,
		vector<double>&				kseedsb,
		vector<double>&				kseedsx,
		vector<double>&				kseedsy,
		std::vector<int>&		    klabels,
		const int&					STEP,
		const int&                  m = 10,
		const int&					NUMITR = 10);
	//============================================================================
	// Pick seeds for superpixels when step size of superpixels is given.
	//============================================================================
	void GetLABXYSeeds_ForGivenStepSize(
		vector<double>&				kseedsl,
		vector<double>&				kseedsa,
		vector<double>&				kseedsb,
		vector<double>&				kseedsx,
		vector<double>&				kseedsy,
		const int&					STEP,
		const bool&					perturbseeds,
		const vector<double>&		edgemag);
	//============================================================================
	// Pick seeds for superpixels when number of superpixels is input.
	//============================================================================
	void GetLABXYSeeds_ForGivenK(
		vector<double>&				kseedsl,
		vector<double>&				kseedsa,
		vector<double>&				kseedsb,
		vector<double>&				kseedsx,
		vector<double>&				kseedsy,
		const int&					STEP,
		const bool&					perturbseeds,
		const vector<double>&		edges);

	//============================================================================
	// Move the seeds to low gradient positions to avoid putting seeds at region boundaries.
	//============================================================================
	void PerturbSeeds(
		vector<double>&				kseedsl,
		vector<double>&				kseedsa,
		vector<double>&				kseedsb,
		vector<double>&				kseedsx,
		vector<double>&				kseedsy,
		const vector<double>&		edges);
	//============================================================================
	// Detect color edges, to help PerturbSeeds()
	//============================================================================
	void DetectLabEdges(const Matrix<float>& lab,vector<double>& edges);
	//============================================================================
	// Post-processing of SLIC segmentation, to avoid stray labels.
	//============================================================================
	void EnforceLabelConnectivity(
		const std::vector<int>&		labels,
		const int&					width,
		const int&					height,
		std::vector<int>&			nlabels,//input labels that need to be corrected to remove stray labels
		int&						numlabels,//the number of labels changes in the end if segments are removed
		const int&					K); //the number of superpixels desired by the user


public:
	int										m_width;
	int										m_height;

    Matrix<float> lab;
};

#endif // !defined(_SLIC_H_INCLUDED_)
