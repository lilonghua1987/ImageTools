#pragma once

#include "ImageExpandTools.h"
#include "ImageProcess.h"
#include "universe.h"
#include "tools.h"

// threshold function
#define THRESHOLD(size, c) (c/size)

struct edge
{
	float w;
    int a, b;

	edge()
		:w(0),a(0),b(0)
	{};

	bool operator<(const edge e)
	{
		return w < e.w;
	};
};



class ImageSegment
{
public:
	ImageSegment(const Image& src, float sigma, float c, int min_size, int& num_ccs);
	~ImageSegment(void);

	universe* segmentGraph(int num_vertices, int num_edges, edge *edges, float c) ;

private:
	void Segment(const Image& src, float sigma, float c, int min_size, int& num_ccs);
};

