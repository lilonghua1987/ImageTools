#include "AdaptiveWeight.h"
#include "ImageProcess.h"
#include "tools.h"
#include <thread>

AdaptiveWeight::AdaptiveWeight(Image imgL,Image imgR)
{
	long height = imgL.height;
	long width = imgL.width;

	this->imgL = imgL;
	this->imgR = imgR;

	T = 40;
	size = 35;
	rc = 5.0;
	rp = (float)(size/2.0);
	dMin = 0;
	dMax = 19;

	Matrix<float> rawCost(height,width,dMax - dMin + 1);
	RawCost(imgL,imgR,rawCost,dMax-dMin);
	std::cout<<"rawCost is end!"<<std::endl;

	Matrix<float> pSW_ref(height,width,size*size);
	Matrix<float> pSW_tar(height,width,size*size);
	AdaptiveSupportWeightComputation(imgL,pSW_ref,(int)rp);
	AdaptiveSupportWeightComputation(imgR,pSW_tar,(int)rp);
	std::cout<<"AdaptiveSupportWeightComputation is end!"<<std::endl;

	greCost = Matrix<float>(height,width,dMax - dMin + 1);
	SupportAggregation(pSW_ref,pSW_tar,rawCost,dMax,(int)rp);
	std::cout<<"SupportAggregation is end!"<<std::endl;

	
}


AdaptiveWeight::~AdaptiveWeight(void)
{
}


inline float AdaptiveWeight::getWeight(Matrix<Pixel<float> > labImg,int i,int j,int m,int n)
{
	float dC = 0;
	float dP = 0;
	float w = 0;

	if((i+m)<0 || (ulong)(i+m)>= labImg.row || (j+n)<0 || (ulong)(j+n)>= labImg.column)
	{
		dC = 0;
		dP = 0;
	}else
	{
		
		Pixel<float> p = labImg.at(i,j);
		Pixel<float> q = labImg.at(i+m,j+n);

		dC = (float)ImageExpandTools::EuclideanDistance(p,q);
		dP = (float)sqrt(m*m+n*n);
		dC/=rc;
		dP/=rp;
		w = exp(-(dC+dP));
	}

	return w;
}

inline float AdaptiveWeight::getTAD(const Image &imgL,const Image &imgR,int i,int j,int d)
{
	float minT = (float)T;
	Pixel<int> p = imgL.get<BYTE>(tools::bound(i,0,imgR.height-1),tools::bound(j,0,imgR.width-1));
	float Lp = (float)p.blue;
	float ap = (float)p.green;
	float bp = (float)p.red;
	Pixel<int> q = imgR.get<BYTE>(tools::bound(i+d,0,imgR.height-1),tools::bound(j+d,0,imgR.width-1));
	float Lq = (float)q.blue;
	float aq = (float)q.green;
	float bq = (float)q.red;

	float Ib = fabsf(Lp-Lq);
	float Ig = fabsf(ap-aq);
	float Ir = fabsf(bp-bq);

	minT = minT<(Ib+Ig+Ir)?minT:(Ib+Ig+Ir);

	return minT;
}


//
//int AdaptiveWeight::computeByThread(int bX,int bY,int eX,int eY)
//{
//	cout<<thread::hardware_concurrency<<endl;
//	cout<<"bX:"<<bX<<" bY:"<<bY<<" eX:"<<eX<<" eY:"<<eY<<endl;
//	int center = floor((size-1)/2.0);
//
//	for(int i=bY;i<std::min(eY,imgL.height);i++)
//	{
//		for(int j=bX;j<std::min(eX,imgL.width);j++)
//		{
//			for(int d=dMin;d<=dMax;d++)
//			{
//				double sumW = 0;
//				double sumAW = 0;
//				for(int m=-center;m<center;m++)
//				{
//					for(int n=-center;n<center;n++)
//					{						
//						double wL = getWeight(labImgL,i,j,m,n);
//						double wR = getWeight(labImgR,i+d,j+d,m,n);
//						sumW += wL*wR;
//						sumAW = wL*wR*getTAD(imgL,imgR,i,j,d);
//					}
//				}
//
//				greCost.at(i,j,d-dMin) = (sumW == 0)?0:sumAW/sumW;
//			}
//		}
//	}
//
//	return rand();
//}


void AdaptiveWeight::optimization()
{
	long height = imgL.height;
	long width = imgL.width;
	Image dispartyMat(width, height, imgL.imgType);

	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width; j++)
		{
			int minDIndex = 0;
			double minE = greCost.at(i,j,0);
			for (int d = dMin; d <= dMax; d++)
			{				
				if (minE > greCost.at(i,j,d-dMin))
				{
					minE = greCost.at(i,j,d-dMin);
					minDIndex = d;
				}
			}
			ImageTools::setPixel(dispartyMat,i,j,minDIndex);
		}
	}

	string imgName = tools::fileName(PATH,"dispartyImg_",width,height,".jpg");
	ImageTools::SaveImage(dispartyMat, imgName.c_str());
}

void AdaptiveWeight::RawCost( const Image &imgL, const Image &imgR, Matrix<float> &rawCost, int range )
{
	assert(imgL.width == imgR.width);
	assert(imgL.height == imgR.height);
	assert(imgL.channel == 3);
	int width = imgL.width,height = imgL.height;

	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width; j++)
		{
			for (int d = 0; d < abs(range); d++)
			{
				rawCost.at(i,j,d-dMin) = getTAD(imgL,imgR,i,j,d-dMin);
			}
		}
	}
}


void AdaptiveWeight::AdaptiveSupportWeightComputation( const Image &img,Matrix<float> &pSW,int mask )
{
	// color space conversion
	Matrix<Pixel<float> > labImg = ColorConversion::ColorSpaceConversionFromRGB2Lab(img);

	// computation
	for(ulong i=0; i<img.height; i++)
	{
		for(ulong j=0; j<img.width; j++)
		{			
			for(int m=-mask,k=0;m<mask;m++)
			{
				if ((i+m) <0 || (ulong)(i+m)>= img.height)
				{
					for (int n=-mask;n<mask;n++)
					{
						int wIndex = (m+mask)*size+(n+mask);
						pSW.at(i,j,wIndex) = 0;
					}					
				}else
				{
					for(int n=-mask;n<mask;n++)
					{	
						int wIndex = (m+mask)*size+(n+mask);
						if(pSW.at(i,j,wIndex)>0)
							continue;
						if ((j+n) <0 || (ulong)(j+n)>= img.width)
						{
							pSW.at(i,j,wIndex) = 0;
						}else
						{
							pSW.at(i,j,wIndex) = getWeight(labImg,i,j,m,n);
							if ((i+mask) < img.height && (j+mask) < img.width)
							{
								pSW.at(i+mask,j,mask) = pSW.at(i,j,wIndex);
								pSW.at(i,j+mask,(mask-1)*size) = pSW.at(i,j,wIndex);
								pSW.at(i+mask,j+mask,0) = pSW.at(i,j,wIndex);
							}
						}						
					}
				}				
			}
		}
	}
}


	
void AdaptiveWeight::SupportAggregation(Matrix<float> pSW_ref,Matrix<float> pSW_tar, Matrix<float> pRawCost,int range,int mask )
{
	for(ulong i=0;i<pSW_ref.row;i++)
	{
		for(ulong j=0;j<pSW_ref.column;j++)
		{
			for(int d=dMin;d<=dMax;d++)
			{
				float sumW = 0;
				float sumAW = 0;
				for(int m=-mask;m<mask;m++)
				{
					for(int n=-mask;n<mask;n++)
					{	
						int wIndex = (m+mask)*size+(n+mask);
						float wR = pSW_tar.at(i,j,wIndex);
						float wL = pSW_ref.at(i,j+d,wIndex);
						sumW += wL*wR;
						sumAW += wL*wR*pRawCost.at(i,j,d-dMin);
					}
				}

			    greCost.at(i,j,d-dMin) = (sumW == 0)?0:sumAW/sumW;
			}
		}
	}
}
