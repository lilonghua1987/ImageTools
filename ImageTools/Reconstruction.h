#pragma once

#include "ImageExpandTools.h"
#include "tools.h"
#include "Log.h"

class Reconstruction
{
public:
	Reconstruction(void);
	~Reconstruction(void);

	static void shading(const Image& img, int iter, float Ps, float Qs);
};