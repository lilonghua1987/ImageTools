#include "MrfStereo.h"


MrfStereo::MrfStereo(const Image& imgL,const Image& imgR,const mrf::Parameters& paramet)
	:imgL(imgL)
	,imgR(imgR)
	,paramet(paramet)
	,mrf(nullptr)
	,dcost(nullptr)
	,scost(nullptr)
	,truncate(40)
	,size(2)
{
	costRaw = Matrix<float>(imgL.height,imgL.width,paramet.nD);
	labL = Matrix<float>(imgL.height,imgL.width,imgL.channel);
	ColorConversion::RGBToLab(imgL,labL);
	labR = Matrix<float>(imgR.height,imgR.width,imgR.channel);
	ColorConversion::RGBToLab(imgR,labR);
}


MrfStereo::~MrfStereo(void)
{
	if(mrf != nullptr) delete mrf;
}


void MrfStereo::computeDSI(const Image& imL, const Image& imR, MRF::CostVal *&dsi, int nD, int birchfield, int squaredDiffs, int truncDiffs)
{
	int width = (int)imL.width, height = (int)imL.height, nB = imL.channel;
	dsi = new MRF::CostVal[width * height * nD];

	int nColors = __min(3, nB);

	// worst value for sumdiff below 
	int worst_match = nColors * (squaredDiffs ? 255 * 255 : 255);
	// truncation threshold - NOTE: if squared, don't multiply by nColors (Eucl. dist.)
	int maxsumdiff = squaredDiffs ? truncDiffs * truncDiffs : nColors * abs(truncDiffs);
	// value for out-of-bounds matches
	int badcost = __min(worst_match, maxsumdiff);

	int dsiIndex = 0;
	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
			for (int d = 0; d < nD; d++) {
				int x2 = x-d;
				int dsiValue;
				if (x2 >= 0 && d < nD) { // in bounds
					int sumdiff = 0;
					for (int b = 0; b < nColors; b++) {
						int diff = 0;
						if (birchfield) {
							// Birchfield/Tomasi cost
							int im1c = imL.atVal<BYTE>(y,x,b);
							int im1l = x == 0?   im1c : (im1c + imL.atVal<BYTE>(y,x-1,b)) / 2;
							int im1r = x == width-1? im1c : (im1c + imL.atVal<BYTE>(y,x+1,b)) / 2;
							int im2c = imR.atVal<BYTE>(y,x2,b);
							int im2l = x2 == 0?   im2c : (im2c + imR.atVal<BYTE>(y,x2-1,b)) / 2;
							int im2r = x2 == width-1? im2c : (im2c + imR.atVal<BYTE>(y,x2+1,b)) / 2;
							int min1 = __min(im1c, __min(im1l, im1r));
							int max1 = __max(im1c, __max(im1l, im1r));
							int min2 = __min(im2c, __min(im2l, im2r));
							int max2 = __max(im2c, __max(im2l, im2r));
							int di1 = __max(0, __max(im1c - max2, min2 - im1c));
							int di2 = __max(0, __max(im2c - max1, min1 - im2c));
							diff = __min(di1, di2);
						} else {
							// simple absolute difference
							int di = imL.atVal<BYTE>(y,x,b) - imR.atVal<BYTE>(y,x,b);
							diff = abs(di);
						}
						// square diffs if requested (Birchfield too...)
						sumdiff += (squaredDiffs ? diff * diff : diff);
					}
					// truncate diffs
					dsiValue = __min(sumdiff, maxsumdiff);
				} else { // out of bounds: use maximum truncated cost
					dsiValue = badcost;
				}
				//int x0=-140, y0=-150;
				//if (x==x0 && y==y0)
				//    printf("dsi(%d,%d,%2d)=%3d\n", x, y, d, dsiValue); 

				// The cost of pixel p and label l is stored at dsi[p*nLabels+l]
				dsi[dsiIndex++] = (MRF::CostVal)dsiValue;
			}
		}
	}
}


void MrfStereo::computeCues(const Image& img, MRF::CostVal *&hCue, MRF::CostVal *&vCue,int gradThresh, int gradPenalty)
{
	int width = (int)img.width, height = (int)img.height, channel = img.channel;
	hCue = new MRF::CostVal[width * height];
	vCue = new MRF::CostVal[width * height];

	// below we compute sum of squared colordiffs, so need to adjust threshold accordingly (results in RMS)
	gradThresh *= channel * gradThresh;

	int n = 0;
	for (int y = 0; y < height; y++) 
	{
		for (int x = 0; x < width; x++) 
		{
			int sx = 0, sy = 0;
			for (int b = 0; b < channel; b++) 
			{
				int dx = img.atVal<BYTE>(y,x,b) - img.atVal<BYTE>(y,x + (x < width-1),b);
				int dy = img.atVal<BYTE>(y,x,b) - img.atVal<BYTE>(y + (y < height-1),x,b);
				sx += dx * dx;
				sy += dy * dy;
			}
			hCue[n] = (MRF::CostVal)(sx < gradThresh ? gradPenalty : 1);
			vCue[n] = (MRF::CostVal)(sy < gradThresh ? gradPenalty : 1);
			//hc.Pixel(x, y, 0) = 100*hCue[n];
			//vc.Pixel(x, y, 0) = 100*vCue[n];
			n++;
		}
	}
}



void MrfStereo::WTA(MRF::CostVal *dsi, int nD, Image &disp)
{
	int n = 0;
	for (int y = 0; y < (int)disp.height; y++)
	{
		for (int x = 0; x < (int)disp.width; x++)
		{
			int minval = (int)dsi[n++]; // dsi(x,y,0)
			int mind = 0;
			for (int d = 1; d < nD; d++)
			{
				int val = (int)dsi[n++]; // dsi(x,y,d)
				if (val < minval) 
				{
					minval = val;
					mind = d;
				}
			}
			disp.at<BYTE>(y,x) = mind;
		}
	}
}



void MrfStereo::getDisparities(MRF *mrf, Image &disp)
{
	int n = 0;
	for (int y = 0; y < (int)disp.height; y++)
	{
		for (int x = 0; x < (int)disp.width; x++)
		{
			disp.at<BYTE>(y,x) = mrf->getLabel(n++);
		}
	}
}


void MrfStereo::setDisparities(const Image& disp, MRF *mrf)
{

	int n = 0;
	for (int y = 0; y < (int)disp.height; y++)
	{
		for (int x = 0; x < (int)disp.width; x++) 
		{
			mrf->setLabel(n++, disp.at<BYTE>(y,x));
		}
	}
}


void MrfStereo::dataCost()
{
	Point2D<long> pL,pR;
	for (long i = 0; i < (long)labL.row; i++)
	{
		for (long j = 0; j < (long)labL.column; j++)
		{
			pL = Point2D<long>(j,i);
			for (int d = 0; d < abs(paramet.nD); d++)
			{
				pR = Point2D<long>(j-d,i);

				if (pR < Point2D<long>(imgL.width,imgL.height) && pR >= Point2D<long>(0,0))
				{
					double sumW = 0;
					double sumAW = 0;

					for (auto i = -size; i <= size; i++)
					{
						for (auto j = -size; j <= size; j++)
						{
							Point2D<long> pwL(pL.x+j,pL.y+i);
							Point2D<long> pwR(pR.x+j,pR.y+i);
							if (pwL < Point2D<long>(imgL.width,imgL.height) && pwL >= Point2D<long>(0,0) && pwR < Point2D<long>(imgR.width,imgR.height) && pwR >= Point2D<long>(0,0))
							{
								double wL = weight(labL,pL,pwL);
								double wR = weight(labR,pR,pwR);
								sumW += wL*wR;
								sumAW += wL*wR*rawCost(pwL,pwR);
							}
						}
					}
					costRaw.at(pL.y,pL.x,abs(d)) = (float)((sumW == 0)?0:sumAW/sumW);
				}
			}
		}
	}
}


void MrfStereoMatch(const Image& imL, const Image& imR)
{
	const char *algs[] = {"ICM", "Expansion", "Swap", "TRW-S", "BP-S", "BP-M"};
	static const int MAXITER = 500;
	mrf::Parameters paramet;
	paramet.nD = 60;
	paramet.MRFalg = aExpansion;
	paramet.gradThresh = -1;
	paramet.smoothexp = 2;
	paramet.smoothmax = 5;
	//paramet.lambda = 20;
	MrfStereo match(imL,imR,paramet);
	Image disp(imL.width,imL.height);
	MRF::CostVal *dsi = NULL;
	match.computeDSI(imL, imR, dsi, paramet.nD, paramet.birchfield, paramet.squaredDiffs, paramet.truncDiffs);
	/*match.dataCost();
	dsi = match.costRaw.date;*/
	match.dcost = new DataCost(dsi);

	MRF::CostVal *hCue = NULL, *vCue = NULL;

	if (paramet.gradThresh > 1) 
	{
		match.computeCues(imL, hCue, vCue, paramet.gradThresh, paramet.gradPenalty);
		match.scost = new SmoothnessCost(paramet.smoothexp, paramet.smoothmax, paramet.lambda, hCue, vCue);
	} else 
	{
		match.scost = new SmoothnessCost(paramet.smoothexp, paramet.smoothmax, paramet.lambda);
	}
	EnergyFunction *energy = new EnergyFunction(match.dcost, match.scost);


	int outerIter, innerIter;
	int noChange = 0;
	for (int numAlg = aICM; numAlg <= aBPM; numAlg++) 
	{
		outerIter = MAXITER;
		innerIter = 3;
		if (paramet.MRFalg < 9 && numAlg != paramet.MRFalg) continue;
		switch (numAlg) 
		{
		case aICM:       match.mrf = new ICM(imL.width, imL.height, paramet.nD, energy); innerIter = 5; break;
		case aExpansion: match.mrf = new Expansion(imL.width, imL.height, paramet.nD, energy); break;
		case aSwap:      match.mrf = new Swap(imL.width, imL.height, paramet.nD, energy); break;
		case aTRWS:      match.mrf = new TRWS(imL.width, imL.height, paramet.nD, energy); break;
		case aBPS:       match.mrf = new BPS(imL.width, imL.height, paramet.nD, energy); break;
		case aBPM:       match.mrf = new MaxProdBP(imL.width, imL.height, paramet.nD, energy); outerIter = MAXITER/2; break;
		default: break;
		}
		if (RUNTIMECHECK)
			printf("******* Running %s for up to %d x %d iterations\n",algs[numAlg], outerIter, innerIter);

		match.mrf->initialize();

		bool initializeToWTA = (numAlg == aICM); 
		if (initializeToWTA) 
		{
			if (RUNTIMECHECK)
				printf("performing WTA\n");

			match.WTA(dsi,paramet.nD, disp);
			match.setDisparities(disp, match.mrf);
		} else
		{
			match.mrf->clearAnswer();
		}


		MRF::EnergyVal E, Ed, Es, Eold;
		Ed = match.mrf->dataEnergy();
		Es = match.mrf->smoothnessEnergy();
		E = Ed + Es; // mrf->totalEnergy();
		Eold = E;
		if (RUNTIMECHECK)
		{
			printf("Energy = %d (Ed=%d, Es=%d) at start\n", E, Ed, Es);
			printf("sec,%s,%s\n", algs[numAlg], numAlg == aTRWS ? "lowerBound," : "");
			printf("0,%d,%s\n", E, numAlg == aTRWS ? "," : "");
		}

		float t, tot_t = 0;
		double lowerBound = 0;
		int iter;
		for (iter = 0; iter < outerIter; iter++) 
		{
			match.mrf->optimize(innerIter, t);
			tot_t += t ;

			Ed = match.mrf->dataEnergy();
			Es = match.mrf->smoothnessEnergy();
			E = Ed + Es; // mrf->totalEnergy();
			if (numAlg == aTRWS)
				lowerBound = match.mrf->lowerBound();
			if (RUNTIMECHECK)
				printf("Energy = %d (Ed=%d, Es=%d)", E, Ed, Es);
			if (numAlg == aTRWS && RUNTIMECHECK)
				printf(", lower bound = %.1f", lowerBound);
			if (RUNTIMECHECK)
				printf(", %.1f secs\n", tot_t);

			if (RUNTIMECHECK) {
				if (numAlg == aTRWS)
					printf("%.1f,%d,%.1f,\n", tot_t, E, lowerBound);
				else
					printf("%.1f,%d,\n", tot_t, E);
			}

			if (E == Eold) {
				if (numAlg <= aSwap) // ICM, Expansion and Swap converge
					break;
				noChange++;
				if (noChange >= 10) // if energy hasn't changed for 10 iterations, it's save to quit
					break;
			} else
				noChange = 0;

			// TODO: print warning if GC energy increases, which can happen for some energies

			Eold = E;
		}

		match.getDisparities(match.mrf,disp);

		for (ulong i = 0; i < disp.height; i++)
		{
			for (ulong j = 0; j < disp.width; j++)
			{
				disp.at<BYTE>(i, j) = (BYTE)(disp.at<BYTE>(i, j) * 4);
			}
		}

		//ImageTools::NormalImage(disp,0,255);

		std::string fileName("temp/");
		fileName.append(algs[numAlg]);
		fileName.append(".png");
		disp.save(fileName);
	}
}