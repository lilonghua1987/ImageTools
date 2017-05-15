#include "StereoMatch.h"

double StereoMatch::dataCost(const Point2D<long>& pL, const Point2D<long>& pR)
{
	double costT = CUTOFF;

	if (pR < Point2D<long>(imgL.width,imgL.height) && pR >= Point2D<long>(0,0))
	{
		double sumW = 0;
		double sumAW = 0;
		int d = pR.x - pL.x;
		costT = costRaw.at(pL.y,pL.x,abs(d));

		if (costT == 0)
		{
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
			costRaw.at(pL.y,pL.x,abs(d)) = float((sumW == 0)?0:sumAW/sumW);
			costT = (sumW == 0) ? 0 : sumAW / sumW;
		} 
	}
	
	return costT;
}


double StereoMatch::dataBTCost(const Point3D<long>& l, const Point3D<long>& r)
{
	// Birchfield/Tomasi cost
	int im1c = imgL.atVal<BYTE>(l.y,l.x,(BYTE)l.z);
	int im1l = l.x == 0?   im1c : (im1c + imgL.atVal<BYTE>(l.y,l.x-1,(BYTE)l.z)) / 2;
	int im1r = l.x == imgL.width-1? im1c : (im1c + imgL.atVal<BYTE>(l.y,l.x+1,(BYTE)l.z)) / 2;
	int im2c = imgR.atVal<BYTE>(r.y,r.x,(BYTE)r.z);
	int im2l = r.x == 0?   im2c : (im2c + imgR.atVal<BYTE>(r.y,r.x-1,(BYTE)r.z)) / 2;
	int im2r = r.x == imgR.width-1? im2c : (im2c + imgR.atVal<BYTE>(r.y,r.x+1,(BYTE)r.z)) / 2;
	int min1 = tools::Min(im1c, tools::Min(im1l, im1r));
	int max1 = tools::Max(im1c, tools::Max(im1l, im1r));
	int min2 = tools::Min(im2c, tools::Min(im2l, im2r));
	int max2 = tools::Max(im2c, tools::Max(im2l, im2r));
	int di1 = tools::Max(0, tools::Max(im1c - max2, min2 - im1c));
	int di2 = tools::Max(0, tools::Max(im2c - max1, min1 - im2c));

	return tools::Min(di1, di2);
}


double StereoMatch::dataPenalty(const Point2D<long>& l, const Point2D<long>& r)
{

	double d_sum = 0.0;
	for (auto i = 0; i < imgL.channel; i++)
	{
		double d = 0;

		switch (params.data_cost)
		{
		case Parameters::L1:
			{
				d = imgL.atVal<BYTE>(l.y,l.x,i) - imgR.atVal<BYTE>(r.y,r.x,i);
				if (d<0) d = -d;
				d_sum += tools::Min(d,double(CUTOFF));
			}
			break;
		case Parameters::L2:
			{
				d = imgL.atVal<BYTE>(l.y,l.x,i) - imgR.atVal<BYTE>(r.y,r.x,i);
				d = d*d;
				d_sum += tools::Min(d,double(CUTOFF));
			}
			break;
		case Parameters::BT:
			{
				d = dataBTCost(Point3D<long>(l.x,l.y,i),Point3D<long>(r.x,r.y,i));
				d_sum += d*d;
			}
		default:
			{
				d = imgL.atVal<BYTE>(l.y,l.x,i) - imgR.atVal<BYTE>(r.y,r.x,i);
				if (d<0) d = -d;
				d_sum += tools::Min(d,double(CUTOFF));
			}
			break;
		}
	}

	if(params.data_cost == Parameters::BT)
		return tools::Min(d_sum / imgL.channel, 255.0 * 255);
	else
		return d_sum / imgL.channel;
}


double StereoMatch::dataPenaltySubpixel(const Point2D<long>& l, const Point2D<long>& r)
{
	double dl, dr, d_sum = 0;
	double Il, Il_min, Il_max, Ir, Ir_min, Ir_max;

	for (auto i = 0; i < imgL.channel; i++)
	{
		Il     = imgL.atVal<BYTE>(l.y,l.x,i); Ir     = imgR.atVal<BYTE>(r.y,r.x,i);
		Il_min = imgLMin.atVal<BYTE>(l.y,l.x,i); Ir_min = imgRMin.atVal<BYTE>(r.y,r.x,i);
		Il_max = imgLMax.atVal<BYTE>(l.y,l.x,i); Ir_max = imgRMax.atVal<BYTE>(r.y,r.x,i);

		if      (Il < Ir_min) dl = Ir_min - Il;
		else if (Il > Ir_max) dl = Il - Ir_max;
		else dl = 0;

		if      (Ir < Il_min) dr = Il_min - Ir;
		else if (Ir > Il_max) dr = Ir - Il_max;
		else dr = 0;

		double d = tools::Min(dl, dr); 
		if (params.data_cost==Parameters::L2) d = d*d;
		if (d>CUTOFF) d = CUTOFF;
		d_sum += d;
	}

	return d_sum/imgL.channel;
}


double StereoMatch::smoothPenalty(const Matrix<float>& img, const std::vector<int>& labels, const Point2D<long>& p, const Point2D<long>& np, long disp, long ndisp)
{
	double d_max = 0, R = 0;

	//if (disp == ndisp) return 0;
	if (disp == OCCLUDED || ndisp == OCCLUDED) R = 0;
	else
	{
		R = abs(disp - ndisp);
		//if (R > params.interaction_radius) R = params.interaction_radius;
	}


	for (ulong i = 0; i< img.channel; i++)
	{
		/*double d = abs(img.at(p.y,p.x,i) - img.at(np.y,np.x,i));
		if (d_max < d) d_max = d;*/
		d_max += abs(img.at(p.y,p.x,i) - img.at(np.y,np.x,i));
	}

	int dS = inSegment(labels, p, np, (long)img.column);
	//if (R > 2) std::cout << "R=" << R << " dS=" << dS << std::endl;
	if (dS > 2)
	{
		return R*params.lambda1;
	}else
	{
		if (R > params.interaction_radius * 2) R = params.interaction_radius;
		if (d_max < params.I_threshold) 
			return R*params.lambda1;
		else                          
			return R*params.lambda2;
	}
}