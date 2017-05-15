#include "StereoMatch.h"

inline int gcd(int a, int b)
{
	if (a < b) return gcd(b, a);
	if (b == 0) return a;
	int r = a % b;
	if (r == 0) return b;
	return gcd(b, r);
}

inline void generatePermutation(std::vector<int>& buf)
{
	for (size_t i = 0; i < buf.size(); i++)
	{
		buf.at(i) = i;
	}

	int j = 0;
	for (size_t i = 0; i < buf.size() - 1; i++)
	{
		j = i + (int) (((double)rand()/(RAND_MAX+1.0))*(buf.size() - i));
		int tmp = buf.at(i); 
		buf.at(i) = buf.at(j); 
		buf.at(j) = tmp;
	}
}


StereoMatch::StereoMatch(const Image& imgL,const Image& imgR)
	:imgL(imgL)
	,imgR(imgR)
	,E(0.0)
	,numLabels(0)
	,mSpcount(1000)
	,mCompactness(20)
	,truncate(40.0)
	,size(3)
{
	if(!imgL.data || !imgR.data) return;
	initParams();
	//RawCosts();
}


StereoMatch::~StereoMatch(void)
{
}


void StereoMatch::initParams()
{
	disp_base = 0; disp_max = 0; disp_size = 1;
	imgSize.x = (long)imgL.width;
	imgSize.y = (long)imgL.height;

	leftX = Matrix<long>(imgSize.y,imgSize.x);
	rightX = Matrix<long>(imgSize.y,imgSize.x);

	if (!leftX.data || !rightX.data)
	{ 
		fprintf(stderr, "Not enough memory!\n"); exit(1); 
	}

	for (long y = 0; y < imgSize.y; y++)
	{
		for (long x = 0; x < imgSize.x; x++)
		{
			leftX.at(y,x) = rightX.at(y,x) = OCCLUDED;
		}
	}

	uniqueFlag = true;

	ptrImL =  Matrix<Energy::Var>(imgSize.y,imgSize.x);
	ptrImR =  Matrix<Energy::Var>(imgSize.y,imgSize.x);

	/*segmL = imgL;
	segmR = imgR;*/
}


void StereoMatch::InitSubPixel()
{
	if (params.sub_pixel && imgL.data && !imgLMin.data)
	{
		imgLMin  = imgL;
		imgLMax  = imgL;
		imgRMin = imgR;
		imgRMax = imgR;

		setSubPixel(imgL,  imgLMin,  imgLMax);
		setSubPixel(imgR, imgRMin, imgRMax);
	}
}


void StereoMatch::setSubPixel(const Image& Im, Image& ImMin, Image& ImMax)
{
	Point2D<long> p;
	//int I, I1, I2, I3, I4, I_min, I_max;

	Pixel<BYTE> pixelM,pixelMin,pixelMax;
	const int neighborNum = 4;
	Pixel<BYTE> neighbors[neighborNum];
	const Point2D<long> neighborhoods[] = {Point2D<long>(-1,0),Point2D<long>(1,0),Point2D<long>(0,-1),Point2D<long>(0,1)};
	Pixel<int> meanI[neighborNum];

	for (p.y=0; p.y < imgSize.y; p.y++)
	{
		for (p.x=0; p.x < imgSize.x; p.x++)
		{
			pixelM = pixelMin = pixelMax = Im.get<BYTE>(p);

			for (int i = 0; i < neighborNum; i++)
			{
				neighbors[i] = Im.get<BYTE>(tools::bound((p+neighborhoods[i]),Point2D<long>(0,0),Point2D<long>(imgSize.x-1,imgSize.y-1)));
				meanI[i] = (Pixel<int>(pixelM) + Pixel<int>(neighbors[i]))/2;

				/* red component */
				if(pixelMin.red > meanI[i].red) pixelMin.red = meanI[i].red;
				if(pixelMax.red < meanI[i].red) pixelMax.red = meanI[i].red;

				if (imgL.channel == 3)
				{
					/* green component */
					if(pixelMin.green > meanI[i].green) pixelMin.green = meanI[i].green;
					if(pixelMax.green < meanI[i].green) pixelMax.green = meanI[i].green;

					/* blue component */
					if(pixelMin.blue > meanI[i].blue) pixelMin.blue = meanI[i].blue;
					if(pixelMax.blue < meanI[i].blue) pixelMax.blue = meanI[i].blue;
				}
			}

			ImMin.set(p,pixelMin);
			ImMax.set(p,pixelMax);
		}
	}
}


void StereoMatch::SetDispRange(int _disp_base, int _disp_max)
{
	disp_base = _disp_base;
	disp_max = _disp_max;
	disp_size = disp_max - disp_base + 1;
	if (! (disp_base <= disp_max) ) { fprintf(stderr, "Error: wrong disparity range!\n"); exit(1); }
	printf("Disparity range: x_min = %d, x_max = %d\n", disp_base, disp_max);
	//RawCosts();
}


void StereoMatch::swapImage()
{
	long c_tmp;
	Image r_tmp;
	Matrix<long> l_tmp;

	printf("SWAP_IMAGES\n\n");

	c_tmp = disp_base;
	disp_base = -disp_max;
	disp_max = -c_tmp;

	r_tmp = imgL;
	imgL = imgR;
	imgR = r_tmp;

	r_tmp = imgLMin;
	imgLMin = imgRMin;
	imgRMin = r_tmp;

	r_tmp = imgLMax;
	imgLMax = imgRMax;
	imgRMax = r_tmp;

	l_tmp = leftX;
	leftX = rightX;
	rightX = l_tmp;
}


void StereoMatch::makeUnique()
{
	Point2D<long> p, d, pd, d2;

	printf("MAKE_UNIQUE\n\n");

	for (p.y = 0; p.y < imgSize.y; p.y++)
	{
		for (p.x=0; p.x<imgSize.x; p.x++)
		{
			rightX.at(p.y,p.x) = OCCLUDED;
		}
	}

	for (p.y=0; p.y<imgSize.y; p.y++)
	{
		for (p.x=0; p.x<imgSize.x; p.x++)
		{
			d.x = leftX.at(p.y,p.x); if (d.x == OCCLUDED) continue;

			pd = p + d;
			if (pd>=Point2D<long>(0,0) && pd<imgSize)
			{
				d2.x = rightX.at(pd.y,pd.x);

				if (d2.x != OCCLUDED)
				{
					leftX.at((pd+d2).y,(pd+d2).x) = OCCLUDED;
				}

				rightX.at(pd.y,pd.x) = -d.x;
			}
			else leftX.at(p.y,p.x) = OCCLUDED;
		}
	}

	uniqueFlag = true;
}


void StereoMatch::crossCheck()
{
	Point2D<long> p, pd;

	printf("CROSS_CHECK\n\n");

	for (p.y=0; p.y<imgSize.y; p.y++)
	{
		for (p.x=0; p.x<imgSize.x; p.x++)
		{
			long d = leftX.at(p.y,p.x); if (d == OCCLUDED) continue;

			pd.x = p.x + d; pd.y = p.y;
			if (pd>=Point2D<long>(0,0) && pd<imgSize)
			{
				long d2 = rightX.at(pd.y,pd.x);
				if (-d == d2) continue;
			}
			leftX.at(p.y,p.x) = OCCLUDED;
		}
	}

	for (p.y=0; p.y<imgSize.y; p.y++)
	{
		for (p.x=0; p.x<imgSize.x; p.x++)
		{
			long d = rightX.at(p.y,p.x); if (d == OCCLUDED) continue;

			pd.x = p.x + d; pd.y = p.y;
			if (pd>=Point2D<long>(0,0) && pd<imgSize)
			{
				long d2 = leftX.at(pd.y,pd.x);
				if (-d == d2) continue;
			}
			rightX.at(p.y,p.x) = OCCLUDED;
		}
	}

	uniqueFlag = true;
}


void StereoMatch::fillOcclusion()
{
	Point2D<long> p, q;

	printf("FILL_OCCLUSIONS\n\n");

	for (p.y=0; p.y<imgSize.y; p.y++)
	{
		long d = 0;
		for (p.x=0; p.x<imgSize.x; p.x++)
		{
			d = leftX.at(p.y,p.x);
			if (d != OCCLUDED) break;
		}
		if (p.x == imgSize.x) continue;
		for (q.x=0, q.y=p.y; q.x<p.x; q.x++)
		{
			leftX.at(q.y,q.x) = d;
		}

		for (; p.x<imgSize.x; p.x++)
		{
			if (leftX.at(p.y,p.x) == OCCLUDED)
			{
				leftX.at(p.y,p.x) = d;
			}
			else
			{
				d = leftX.at(p.y,p.x);
			}
		}
	}

	uniqueFlag = false;
}


double StereoMatch::GetK()
{
	Point2D<long> p_base;
	long d = 0;
	Point2D<long> p;
	double delta, sum = 0, num = 0;
	std::vector<double> array;

	int k = (disp_size + 2) / 4; /* 0.25 times the number of disparities */
	if (k < 3) k = 3; if (k > disp_size) k = disp_size;

	array = std::vector<double>(k);

	p_base.x = tools::Max(-disp_base, 0);
	p_base.y = 0;

	for (p.y = p_base.y; p.y < imgSize.y; p.y++)
	{
		for (p.x = p_base.x; p.x < imgSize.x - disp_max; p.x++)
		{
			/* compute k'th smallest value among dataPenalty(p, p+d) for all d */
			int i = 0;

			for (d = disp_base; d <= disp_max; d++)
			{
				if (params.sub_pixel)
				{
					delta = dataPenaltySubpixel(p, Point2D<long>(p.x + d,p.y));
				}
				else
				{
					delta = dataPenalty(p, Point2D<long>(p.x + d, p.y));
					//delta = dataCost(p, p+d);
				}
				if (i < k) array[i++] = delta;
				else
				{
					for (int i = 0; i< k; i++)
						if (delta < array[i])
						{
							double tmp = delta;
							delta = array[i];
							array[i] = tmp;
						}
				}
			}

			delta = array[0];
			for (int i = 1; i < k; i++) if (delta < array[i]) delta = array[i];

			sum += delta;
			num ++;
		}
	}

	if (num == 0) { fprintf(stderr, "GetK: Not enough samples!\n"); exit(1); }
	if (sum == 0) { fprintf(stderr, "GetK failed: K is 0!\n"); exit(1); }

	double K = sum/num;
	printf("Computing statistics: dataPenalty noise is %f\n", K);
	return K;
}


void StereoMatch::setParameters(const Parameters& params)
{
	if (labelsL.size() == 0 && labelsR.size() == 0)
	{
		SLIC slic;
		slic.PerformSLICO_ForGivenK(imgL, labelsL, numLabels, mSpcount, mCompactness);
		Image imL(imgL);
		slic.DrawContoursAroundSegments(imL, labelsL);
		imL.save(params.saveDir + "slic_seg.png");
		labL = slic.lab;
		slic.PerformSLICO_ForGivenK(imgR, labelsR, numLabels, mSpcount, mCompactness);
		labR = slic.lab;
		slic.~SLIC();
		std::cout << "The image segment is over !" << std::endl;
	}
	this->params = params;
	InitSubPixel();
}


void StereoMatch::SaveXLeft(const std::string& fileName, bool flage)
{
	Point2D<ulong> p;

	printf("Saving left x-disparity map as %s\n\n", fileName);
	Image img(imgSize.x,imgSize.y);

	for (p.y=0; p.y<(ulong)imgSize.y; p.y++)
	{
		for (p.x=0; p.x<(ulong)imgSize.x; p.x++)
		{
			int d = leftX.at(p.y,p.x), c;
			if (d==OCCLUDED) c = 0;
			else
			{
				if (flage) c = d - disp_base;
				else      c = disp_max - d;
			}

			img.at<BYTE>(p.y, p.x) = (uchar)(c * params.scale);
		}
	}

	//ImageTools::NormalImage(img,0,255);
	img.save(fileName);
}


void StereoMatch::SaveScaledXLeft(const std::string& fileName, bool flag /* = false */)
{
	Point2D<long> p;

	printf("Saving scaled left x-disparity map as %s\n\n", fileName);
	Image img(imgSize.x,imgSize.y,FIT_BITMAP,3);
	Pixel<BYTE> point;

	for (p.y=0; p.y<imgSize.y; p.y++)
	{
		for (p.x=0; p.x<imgSize.x; p.x++)
		{
			int d = leftX.at(p.y,p.x), c;
			if (d==OCCLUDED) { point.red = 255; point.green = point.blue = 0; }
			else
			{
				if (disp_size == 0) c = 255;
				else if (flag) c = 255 - (255-64)*(disp_max - d)/disp_size;
				else           c = 255 - (255-64)*(d - disp_base)/disp_size;
				point.red = point.green = point.blue = c;
			}
			img.set(p.y,p.x,point);
		}
	}

	img.save(fileName);
}


void StereoMatch::SaveXRight(const std::string& fileName, bool flage)
{
	Point2D<ulong> p;

	printf("Saving right x-disparity map as %s\n\n", fileName);
	Image img(imgSize.x,imgSize.y);

	for (p.y=0; p.y<(ulong)imgSize.y; p.y++)
	{
		for (p.x=0; p.x<(ulong)imgSize.x; p.x++)
		{
			if (rightX.at(p.y, p.x) == OCCLUDED) img.at<BYTE>(p.y, p.x) = 0;
			else img.at<BYTE>(p.y, p.x) = (uchar)(rightX.at(p.y, p.x) * params.scale);
		}
	}

	//ImageTools::NormalImage(img,0,255);
	img.save(fileName);
}


void StereoMatch::KZ1()
{
	if ( params.K < 0 ||
		params.I_threshold < 0 ||
		params.lambda1 < 0 ||
		params.lambda2 < 0 ||
		params.denominator < 1 )
	{
		fprintf(stderr, "Error in KZ1: wrong parameter!\n");
		exit(1);
	}

	/* printing parameters */
	if (params.denominator == 1)
	{
		printf("KZ1:  K = %d\n", params.K);
		printf("      I_threshold = %d, lambda1 = %d, lambda2 = %d\n",
			params.I_threshold, params.lambda1, params.lambda2);
	}
	else
	{
		printf("KZ1:  K = %d/%d\n", params.K, params.denominator);
		printf("      I_threshold = %d, lambda1 = %d/%d, lambda2 = %d/%d\n",
			params.I_threshold, params.lambda1, params.denominator,
			params.lambda2, params.denominator);
	}
	printf("      sub_pixel = %s, data_cost = L%d\n",
		params.sub_pixel ? "true" : "false", params.data_cost==Parameters::L1 ? 1 : 2);

	if (disp_max <= 0)
		KZ1_visibility = true;
	else
	{
		KZ1_visibility = false;
		printf("Visibility constraint is not enforced! (not a stereo case)\n");
	}

	Run_KZ_BVZ(METHOD_KZ1);

	uniqueFlag = false;
}


void StereoMatch::KZ2()
{
	if ( params.K < 0 ||
		params.I_threshold2 < 0 ||
		params.lambda1 < 0 ||
		params.lambda2 < 0 ||
		params.denominator < 1 )
	{
		fprintf(stderr, "Error in KZ2: wrong parameter!\n");
		exit(1);
	}
	if (!uniqueFlag) makeUnique();

	/* printing parameters */
	if (params.denominator == 1)
	{
		printf("KZ2:  K = %d\n", params.K);
		printf("      I_threshold2 = %d, lambda1 = %d, lambda2 = %d\n",
			params.I_threshold2, params.lambda1, params.lambda2);
	}
	else
	{
		printf("KZ2:  K = %d/%d\n", params.K, params.denominator);
		printf("      I_threshold2 = %d, lambda1 = %d/%d, lambda2 = %d/%d\n",
			params.I_threshold2, params.lambda1, params.denominator,
			params.lambda2, params.denominator);
	}
	printf("      sub_pixel = %s, data_cost = L%d\n",
		params.sub_pixel ? "true" : "false", params.data_cost==Parameters::L1 ? 1 : 2);

	Run_KZ_BVZ(METHOD_KZ2);
}


void StereoMatch::BVZ()
{
	Point2D<long> p;

	if ( params.occlusion_penalty < 0 ||
		params.I_threshold < 0 ||
		params.lambda1 < 0 ||
		params.lambda2 < 0 ||
		params.denominator < 1 )
	{
		fprintf(stderr, "Error in BVZ: wrong parameter!\n");
		exit(1);
	}

	/* printing parameters */
	printf("BVZ:  occlusion_penalty = %d\n", params.occlusion_penalty);
	if (params.denominator == 1)
	{
		printf("      I_threshold = %d, lambda1 = %d, lambda2 = %d\n",
			params.I_threshold, params.lambda1, params.lambda2);
	}
	else
	{
		printf("      I_threshold = %d, lambda1 = %d/%d, lambda2 = %d/%d\n",
			params.I_threshold, params.lambda1, params.denominator,
			params.lambda2, params.denominator);
	}
	printf("      sub_pixel = %s, data_cost = L%d\n",
		params.sub_pixel ? "true" : "false", params.data_cost==Parameters::L1 ? 1 : 2);

	Run_KZ_BVZ(METHOD_BVZ);

	uniqueFlag = false;
}


void StereoMatch::Run_KZ_BVZ(Method method)
{
	unsigned int seed = (unsigned int)time(NULL);
	printf("Random seed = %d\n", seed);
	srand(seed);

	int label_num = disp_size;//最大视差范围
	if ( method == METHOD_BVZ && params.occlusion_penalty < infinityP ) label_num ++;

	switch (method)
	{
		case METHOD_KZ1: KZ1_ComputeEnergy(); break;
		case METHOD_KZ2: KZ2_ComputeEnergy(); break;
		case METHOD_BVZ: BVZ_ComputeEnergy(); break;
	}
	printf("E = %0.4f\n", E);

	/* starting the algorithm */
	std::vector<int> permutation(label_num, 0);/* contains random permutation of 0, 1, ..., label_num-1 */
	std::vector<bool> buf(label_num, false); /* if buf[l] is true then expansion of label corresponding to l cannot decrease the energy */ 
	int buf_num = label_num;/* number of 'false' entries in buf */
	int step = 0;
	double E_old = 0;
	for (int iter = 0; iter < params.iter_max && buf_num > 0; iter++)
	{
		//第一次或者每次迭代，随机产生标签（视差）顺序
		if (iter==0 || params.randomize_every_iteration)
			generatePermutation(permutation);

		for (int index = 0; index < label_num; index++)
		{
			int label = permutation[index];
			if (buf[label]) continue;

			long a = disp_base + label;
			if (a > disp_max) a = OCCLUDED;

			E_old = E;

			switch (method)
			{
			    case METHOD_KZ1: KZ1_Expand(a); break;
				case METHOD_KZ2: KZ2_Expand(a); break;
				case METHOD_BVZ: BVZ_Expand(a); break;
			}

			if(!RUNTIMECHECK) //release do that
			{
				double E_tmp = E;
				switch (method)
				{
				    case METHOD_KZ1: KZ1_ComputeEnergy(); break;
					case METHOD_KZ2: KZ2_ComputeEnergy(); break;
					case METHOD_BVZ: BVZ_ComputeEnergy(); break;
				}
				/*if (E_tmp != E)
				{
					fprintf(stderr, "E and E_tmp are different! (E = %0.4f, E_tmp = %0.4f)\n", E, E_tmp);
					exit(1);
				}*/
			}

			step ++;
			if (E_old == E) printf("-");
			else printf("*");
			fflush(stdout);

			if (E_old == E)
			{
				if (!buf[label]) 
				{ 
					buf[label] = true;
					buf_num --;
				}
			}
			else
			{
				for (int i=0; i<label_num; i++) buf[i] = false;
				buf[label] = true;
				buf_num = label_num - 1;
			}
		}
		printf(" E = %0.4f\n", E); fflush(stdout);
	}

	printf("%.1f iterations\n", ((float)step)/label_num);
}


void StereoMatch::RawCosts()
{
	timer.start();
	costRaw = Matrix<float>(imgL.height,imgL.width,abs(disp_base)+1);
	timer.timeDisplay("RawCosts");
}



void runStereoMatch(const Image& imgL, const Image& imgR, const std::string& dir)
{
	Parameters params = {
		false,Parameters::L2, 1, /* subpixel, data_cost, denominator */
		5, 8, 1, -1, -1,       /* I_threshold, I_threshold2, interaction_radius, lambda1, lambda2 */
		-1,                 /* K */
		StereoMatch::infinityP,           /* occlusion_penalty */
		5, false,    /* iter_max, randomize_every_iteration */
		5,                  /* w_size */
		4,
		dir
	};

	Method method = METHOD_KZ1;

	int lambda = -1, denom = 1, tmp_denom = 1;
	if(method == METHOD_BVZ) lambda = 10;
	lambda *= params.denominator;
	params.lambda1 *= tmp_denom;
	params.lambda2 *= tmp_denom;
	params.K *= tmp_denom;
	params.denominator *= tmp_denom;

	params.iter_max = 9;

	StereoMatch match(imgL,imgR);
	match.SetDispRange(-59,0);

	if (method == METHOD_BVZ)
	{
		if (params.lambda1<0) params.lambda1 = 2*lambda;
		if (params.lambda2<0) params.lambda2 = lambda;
		tmp_denom = gcd(params.lambda1, gcd(params.lambda2, params.denominator));

		if (tmp_denom)
		{
			params.lambda1 /= tmp_denom;
			params.lambda2 /= tmp_denom;
			params.denominator /= tmp_denom;
		}

		match.setParameters(params);
		match.BVZ();

		match.SaveXLeft(params.saveDir+"bvz.png");
		match.SaveScaledXLeft(params.saveDir + "bvz_scaled.png");
	}else
	{
		if (lambda<0 && (params.K<0 || params.lambda1 || params.lambda2<0))
		{
			match.setParameters(params);
			double v = match.GetK()/5; tmp_denom = 1;
			while (v < 3) { v *= 2; tmp_denom *= 2; }
			lambda = int(v + (float)0.5);
			lambda *= params.denominator;
			params.lambda1 *= tmp_denom;
			params.lambda2 *= tmp_denom;
			params.K *= tmp_denom;
			params.denominator *= tmp_denom;
		}
		if (params.K<0) params.K = 5*lambda;
		if (params.lambda1<0) params.lambda1 = 3*lambda;
		if (params.lambda2<0) params.lambda2 = lambda;
		tmp_denom = gcd(params.K, gcd(params.lambda1, gcd(params.lambda2, params.denominator)));
		if (tmp_denom)
		{
			params.K /= tmp_denom;
			params.lambda1 /= tmp_denom;
			params.lambda2 /= tmp_denom;
			params.denominator /= tmp_denom;
		}

		match.setParameters(params);
		
		if(method == METHOD_KZ1) match.KZ1();
		if(method == METHOD_KZ2) match.KZ2();

		match.SaveXLeft(params.saveDir + "KZ1_KZ2.png");
		match.SaveXRight(params.saveDir + "KZ1_KZ2_r.png");
		//match.SaveYLeft("temp/test_y.jpg");
		match.SaveScaledXLeft(params.saveDir + "KZ1_KZ2_scaled.png");
		/*match.CORR();
		match.swapImage();
		match.CORR();
		match.swapImage();*/
		match.crossCheck();
		match.SaveXLeft(params.saveDir + "KZ1_KZ2_cross.png");
		match.SaveScaledXLeft(params.saveDir + "KZ1_KZ2_cross_scaled.png");
		match.fillOcclusion();
		match.SaveXLeft(params.saveDir + "KZ1_KZ2_fill.png");
		match.SaveScaledXLeft(params.saveDir + "KZ1_KZ2_fill_scaled.png");
	}
}