#include "sharpen.h"
#include "tools.h"

sharpen::sharpen(void)
{
}


sharpen::~sharpen(void)
{
}


double* sharpen::laplacianGrad(const Image &img){

	const int img_width = img.width;  //ͼ����
	const int img_height = img.height;  //ͼ��߶�

	double *grade =new double[img_width*img_height]; 
	// ��ʼ��
	for(int y=0; y<img_height ; y++ )
		for(int x=0 ; x<img_width ; x++ )
		{
			*(grade+y*img_width+x)=0;
		}
		// ����ģ��ϵ��
		static int nWeight[3][3] ;
		nWeight[0][0] = -1 ;
		nWeight[0][1] = -1 ;
		nWeight[0][2] = -1 ;
		nWeight[1][0] = -1 ;
		nWeight[1][1] = 8 ;
		nWeight[1][2] = -1 ;
		nWeight[2][0] = -1 ;
		nWeight[2][1] = -1 ;
		nWeight[2][2] = -1 ;
		// �������������ʾLaplacian ��������ֵ
		int nTmp[3][3];
		// ��ʱ����
		double dGrad;
		// ģ��ѭ�����Ʊ���
		int yy ;
		int xx ;
		// ���濪ʼ����Laplacian ���ӽ��м��㣬Ϊ�˱�֤��������Ҫ��
		// ������λ��ͼ�����ݵ��ڲ������������ѭ����������
		// y<img_height-1 ������y<img_height����Ӧ��x ����Ҳ��x<img_width-1
		// ������x<nWidth
		for(int y=1; y<img_height-1 ; y++ ){
			for(int x=1 ; x<img_width-1 ; x++ )
			{
				dGrad = 0 ;
				// Laplacian ������Ҫ�ĸ�������ֵ
				// ģ���һ��
				nTmp[0][0] = ImageTools::getPixel(img,y-1,x-1);
				nTmp[0][1] = ImageTools::getPixel(img,y-1,x);
				nTmp[0][2] = ImageTools::getPixel(img,y-1,x+1) ;
				// ģ��ڶ���
				nTmp[1][0] = ImageTools::getPixel(img,y,x-1);
				nTmp[1][1] = ImageTools::getPixel(img,y,x);
				nTmp[1][2] = ImageTools::getPixel(img,y,x+1);
				// ģ�������
				nTmp[2][0] = ImageTools::getPixel(img,y+1,x-1);
				nTmp[2][1] = ImageTools::getPixel(img,y+1,x);
				nTmp[2][2] = ImageTools::getPixel(img,y+1,x+1);
				// �����ݶ�
				for(yy=0; yy<3; yy++)
					for(xx=0; xx<3; xx++)
					{
						dGrad += nTmp[yy][xx] * nWeight[yy][xx] ;
					}
					// �ݶ�ֵд���ڴ�
					*(grade+y*img_width+x)=dGrad;
			}
		}

	return grade;

}


void sharpen::printLaplacianGrad(double *grad,const Image &templateImg){

	int img_height = templateImg.height;

	int img_width = templateImg.width;

	Image laplacianImg = ImageTools::ImageCopy(templateImg,0,0,img_width,img_height);

	for(int y=0; y<img_height ; y++ ){
		for(int x=0 ; x<img_width ; x++ )
		{
			if(*(grad+y*img_width+x)>50)
				ImageTools::setPixel(laplacianImg,y,x,0);
			else
				ImageTools::setPixel(laplacianImg,y,x,255);
		}
	}

	string imgName = tools::fileNameFromTime(PATH,"laplacianImg",".jpg");
	ImageTools::SaveImage(laplacianImg, imgName.c_str());
}



Image sharpen::laplacian(const Image &img,int k){

	assert(k==4 || k==8);

	const int img_width = img.width;  //ͼ����
	const int img_height = img.height;  //ͼ��߶�

	Image laplacianImg = ImageTools::ImageCopy(img,0,0,img_width,img_height);

	// ����ģ��ϵ��
	int nWeight[3][3] ;

	nWeight[0][0] = 0 ;
	nWeight[0][1] = 1 ;
	nWeight[0][2] = 0 ;
	nWeight[1][0] = 1 ;
	nWeight[1][1] = -k;
	nWeight[1][2] = 1 ;
	nWeight[2][0] = 0 ;
	nWeight[2][1] = 1 ;
	nWeight[2][2] = 0 ;

	if (k == 8)
	{
		nWeight[0][0] = 1 ;
		nWeight[0][2] = 1 ;
		nWeight[2][0] = 1 ;
		nWeight[2][2] = 1 ;
	}

	for(int y=0; y<img_height ; y++ ){
		for(int x=0 ; x<img_width ; x++ )
		{
			Pixel<int> dGrad ;
			// �����ݶ�
			for(int yy=-1; yy<2; yy++)
			{
				for(int xx=-1; xx<2; xx++)
				{
					dGrad += ImageTools::getColorPixel(img,tools::bound(y+yy,0,img_height-1),tools::bound(x+xx,0,img_width-1)) * nWeight[yy+1][xx+1] ;
				}
			}
			if(dGrad<0)
				dGrad = Pixel<int>(0,0,0,0);
			if(dGrad>255)
				dGrad = Pixel<int>(255,255,255,255);
			ImageTools::setColorPixel(laplacianImg,y,x,dGrad);
		}
	}

	/*string imgName = tools::fileNameFromTime(PATH,"laplacian",".jpg");
	ImageTools::SaveImage(laplacianImg, imgName.c_str());*/

	return laplacianImg;

}


//  һά��˹�ֲ�����������ƽ�����������ɵĸ�˹�˲�ϵ��
void sharpen::CreatGauss(double sigma, double **pdKernel, int *pnWidowSize)
{

 LONG i;

 //�������ĵ�
 int nCenter;

 //������һ�㵽���ĵ����
 double dDis;

 //�м����
 double dValue;
 double dSum;
 dSum = 0;

 // [-3*sigma,3*sigma] �������ݣ��Ḳ�Ǿ��󲿷��˲�ϵ��
 *pnWidowSize = static_cast<int>(1+ 2*ceil(3*sigma));

 nCenter = (*pnWidowSize)/2;

 *pdKernel = new double[*pnWidowSize];

 //���ɸ�˹����
 for(i=0;i<(*pnWidowSize);i++)
 {
  dDis = double(i - nCenter);
  dValue = exp(-(1/2)*dDis*dDis/(sigma*sigma))/(sqrt(2*PI)*sigma);
  (*pdKernel)[i] = dValue;
  dSum+=dValue;

 }
 //��һ��
 for(i=0;i<(*pnWidowSize);i++)
 {
  (*pdKernel)[i]/=dSum;
 }

}

//�ø�˹�˲���ƽ��ԭͼ��
Image sharpen::GaussianSmooth(const Image &img,double sigma){

	const int img_width = img.width;  //ͼ����
	const int img_height = img.height;  //ͼ��߶�

	Image resultImg = ImageTools::ImageClone(img);

	LONG x, y;
	LONG i;

	//��˹�˲�������
	int nWindowSize;

	//���ڳ���
	int nLen;

	//һά��˹�˲���
	double *pdKernel;

	//��˹ϵ����ͼ�����ݵĵ��
	double dDotMul;

	//�˲�ϵ���ܺ�
	double dWeightSum;

	double *pdTemp;
	pdTemp = new double[img_width*img_height];

	//����һά��˹����
	CreatGauss(sigma, &pdKernel, &nWindowSize);

	nLen = nWindowSize/2;

	//x�����˲�
	for(y=0;y<img_height;y++)
	{
		for(x=0;x<img_width;x++)
		{
			dDotMul = 0;
			dWeightSum = 0;
			for(i=(-nLen);i<=nLen;i++)
			{
				//�ж��Ƿ���ͼ���ڲ�
				if((i+x)>=0 && (i+x)<img_width)
				{
					dDotMul+=(double)ImageTools::getPixel(img,y,(i+x)) * pdKernel[nLen+i];
					dWeightSum += pdKernel[nLen+i];
				}
			}
			pdTemp[y*img_width+x] = dDotMul/dWeightSum;
		}
	}

	//y�����˲�
	for(x=0; x<img_width;x++)
	{
		for(y=0; y<img_height; y++)
		{
			dDotMul = 0;
			dWeightSum = 0;
			for(i=(-nLen);i<=nLen;i++)
			{
				if((i+y)>=0 && (i+y)< img_height)
				{
					dDotMul += (double)pdTemp[(y+i)*img_width+x]*pdKernel[nLen+i];
					dWeightSum += pdKernel[nLen+i];
				}
			}
			ImageTools::setPixel(resultImg,y,x,static_cast<int>(dDotMul/dWeightSum));
		}
	}

	delete []pdKernel;
	pdKernel = NULL;

	delete []pdTemp;
	pdTemp = NULL;

	string imgName = tools::fileNameFromTime(PATH,"cannyGaussImg",".jpg");
	ImageTools::SaveImage(resultImg, imgName.c_str());

	return resultImg;

}

// ������,���ݶ�
void sharpen::Grad(const Image &img,int *pGradX, int *pGradY, int *pMag){

	const int img_width = img.width;  //ͼ����
	const int img_height = img.height;  //ͼ��߶�
	LONG y,x;

	//x����ķ�����
	for(y=1;y<img_height-1;y++)
	{
		for(x=1;x<img_width-1;x++)
		{
			pGradX[y*img_width +x] =ImageTools::getPixel(img,y,x+1)-ImageTools::getPixel(img,y,x-1);
		}
	}

	//y��������
	for(x=1;x<img_width-1;x++)
	{
		for(y=1;y<img_height-1;y++)
		{
			pGradY[y*img_width +x] = ImageTools::getPixel(img,y+1,x) - ImageTools::getPixel(img,y-1,x);
		}
	}

	//���ݶ�

	//�м����
	double dSqt1;
	double dSqt2;

	for(y=0; y<img_height; y++)
	{
		for(x=0; x<img_width; x++)
		{
			//���׷������ݶ�
			dSqt1 = pGradX[y*img_width + x]*pGradX[y*img_width + x];
			dSqt2 = pGradY[y*img_width + x]*pGradY[y*img_width + x];
			pMag[y*img_width+x] = (int)(sqrt(dSqt1+dSqt2)+0.5);
		}
	}
}



//���������
void sharpen::NonmaxSuppress(int *pMag, int *pGradX, int *pGradY, Image &img){

	const int img_width = img.width;  //ͼ����
	const int img_height = img.height;  //ͼ��߶�

	LONG y,x;
	int nPos;

	//�ݶȷ���
	int gx;
	int gy;

	//�м����
	int g1,g2,g3,g4;
	double weight;
	double dTmp,dTmp1,dTmp2;

	//����ͼ���ԵΪ�����ܵķֽ��
	for(x=0;x<img_width;x++)
	{
		ImageTools::setPixel(img,0,x,0);
		ImageTools::setPixel(img,img_height-1,x,0);

	}
	for(y=0;y<img_height;y++)
	{
		ImageTools::setPixel(img,y,0,0);
		ImageTools::setPixel(img,y,img_width-1,0);
	}

	for(y=1;y<img_height-1;y++)
	{
		for(x=1;x<img_width-1;x++)
		{
			//��ǰ��
			nPos = y*img_width + x;

			//�����ǰ�����ݶȷ���Ϊ0�����Ǳ߽��
			if(pMag[nPos] == 0)
			{
				ImageTools::setPixel(img,y,x,0);
			}
			else
			{
				//��ǰ����ݶȷ���
				dTmp = pMag[nPos];

				//x,y������
				gx = pGradX[nPos];
				gy = pGradY[nPos];

				//���������y������x������˵����������������y����
				if(abs(gy) > abs(gx))
				{
					//�����ֵ����
					weight = abs(gx)/abs(gy);

					g2 = pMag[nPos-img_width];
					g4 = pMag[nPos+img_width];

					//���x,y�����������ķ�����ͬ
					//C Ϊ��ǰ���أ���g1-g4 ��λ�ù�ϵΪ��
					//g1 g2
					//      C
					//       g4 g3
					if(gx*gy>0)
					{
						g1 = pMag[nPos-img_width-1];
						g3 = pMag[nPos+img_width+1];
					}

					//���x,y��������ķ����������෴
					//C�ǵ�ǰ���أ���g1-g4�Ĺ�ϵΪ��
					//       g2 g1
					//        C
					//    g3 g4
					else
					{
						g1 = pMag[nPos-img_width+1];
						g3 = pMag[nPos+img_width-1];
					}
				}

				//���������x������y������˵�������ķ���������x����
				else
				{
					//��ֵ����
					weight = abs(gy)/abs(gx);

					g2 = pMag[nPos+1];
					g4 = pMag[nPos-1];

					//���x,y��������ķ�����������ͬ
					//��ǰ����C�� g1-g4�Ĺ�ϵΪ
					//  g3
					//  g4 C g2
					//       g1
					if(gx * gy > 0)
					{
						g1 = pMag[nPos+img_width+1];
						g3 = pMag[nPos-img_width-1];
					}

					//���x,y�����������ķ����෴
					// C��g1-g4�Ĺ�ϵΪ
					//   g1
					//    g4 C g2
					//     g3
					else
					{
						g1 = pMag[nPos-img_width+1];
						g3 = pMag[nPos+img_width-1];
					}
				}

				//���� g1-g4 ���ݶȽ��в�ֵ
				{
					dTmp1 = weight*g1 + (1-weight)*g2;
					dTmp2 = weight*g3 + (1-weight)*g4;

					//��ǰ���ص��ݶ��Ǿֲ������ֵ
					//�õ�����Ǳ߽��
					if(dTmp>=dTmp1 && dTmp>=dTmp2)
					{
						ImageTools::setPixel(img,y,x,128);
					}
					else
					{
						//�������Ǳ߽��
						ImageTools::setPixel(img,y,x,0);
					}
				}
			}
		}
	}

	string imgName = tools::fileNameFromTime(PATH,"cannyNonmaxSuppressImg",".jpg");
	ImageTools::SaveImage(img, imgName.c_str());
}



// ͳ��pMag��ֱ��ͼ���ж���ֵ
void sharpen::EstimateThreshold(const Image &img,int *pMag,int *pThrHigh, int *pThrLow,double dRatHigh, double dRatLow){

	const int img_width = img.width;  //ͼ����
	const int img_height = img.height;  //ͼ��߶�
	LONG y,x,k;

	//������Ĵ�С���ݶ�ֵ�ķ�Χ�йأ�������ñ�������㷨
	//��ô�ݶȵķ�Χ���ᳬ��pow(2,10)
	int nHist[256];

	//���ܱ߽���
	int nEdgeNum;

	//����ݶ���
	int nMaxMag;

	int nHighCount;

	nMaxMag = 0;

	//��ʼ��
	for(k=0;k<256;k++)
	{
		nHist[k] = 0;
	}
	//ͳ��ֱ��ͼ,����ֱ��ͼ������ֵ
	for(y=0;y<img_height;y++)
	{
		for(x=0;x<img_width;x++)
		{
			if(ImageTools::getPixel(img,y,x)==128)
			{
				nHist[pMag[y*img_width+x]]++;
			}
		}
	}

	nEdgeNum = nHist[0];
	nMaxMag = 0;

	//ͳ�ƾ����������ֵ���ơ����ж�������
	for(k=1;k<256;k++)
	{
		if(nHist[k] != 0)
		{
			nMaxMag = k;
		}

		//�ݶ�Ϊ0�ĵ��ǲ�����Ϊ�߽���
		//����non-maximum suppression���ж�������
		nEdgeNum += nHist[k];

	}

	//�ݶȱȸ���ֵ*pThrHigh С�����ص�����Ŀ
	nHighCount = (int)(dRatHigh * nEdgeNum + 0.5);

	k=1;
	nEdgeNum = nHist[1];

	//�������ֵ
	while((k<(nMaxMag-1)) && (nEdgeNum < nHighCount))
	{
		k++;
		nEdgeNum += nHist[k];
	}

	*pThrHigh = k;

	//����ֵ
	*pThrLow = (int)((*pThrHigh) * dRatLow + 0.5);

}

//���ú���Ѱ�ұ߽����
void sharpen::Hysteresis(Image &img, int *pMag,double dRatLow, double dRatHigh){

	const int img_width = img.width;  //ͼ����
	const int img_height = img.height;  //ͼ��߶�

	LONG y,x;

	int nThrHigh,nThrLow;

	int nPos;
	//����TraceEdge ������Ҫ�ĵ���ֵ���Լ�Hysteresis����ʹ�õĸ���ֵ
	EstimateThreshold(img,pMag,&nThrHigh,&nThrLow,dRatHigh,dRatLow);

	//Ѱ�Ҵ���dThrHigh�ĵ㣬��Щ�����������߽�㣬
	//Ȼ����TraceEdge�������ٸõ��Ӧ�ı߽�
	for(y=0;y<img_height;y++)
	{
		for(x=0;x<img_width;x++)
		{
			nPos = y*img_width + x;

			//����������ǿ��ܵı߽�㣬�����ݶȴ��ڸ���ֵ��
			//��������Ϊһ���߽�����
			if((ImageTools::getPixel(img,y,x)==128) && (pMag[nPos] >= nThrHigh))
			{
				//���øõ�Ϊ�߽��
				ImageTools::setPixel(img,y,x,255);
				TraceEdge(img,y,x,nThrLow,pMag);
			}

		}
	}

	string imgName = tools::fileNameFromTime(PATH,"cannyTraceEdgeImg",".jpg");
	ImageTools::SaveImage(img, imgName.c_str());

	//�������Ѿ�������Ϊ�߽��
	for(y=0;y<img_height;y++)
	{
		for(x=0;x<img_width;x++)
		{
			nPos = y*img_width + x;

			if(ImageTools::getPixel(img,y,x) != 255)
			{
				ImageTools::setPixel(img,y,x,0);
			}
		}
	}

	string imgName2 = tools::fileNameFromTime(PATH,"cannyHysteresisImg",".jpg");
	ImageTools::SaveImage(img, imgName2.c_str());
}

//����Hysteresis ִ�еĽ������һ�����ص㿪ʼ�����������Ը����ص�Ϊ�߽�����һ���߽��
//һ���߽�����б߽�㣬���������˵ݹ��㷨
//       �ӣ�x,y)������������б߽��ĸ��٣�����ֻ����pResult��û�д����ҿ����Ǳ߽�
//  ������أ�=128��������ֵΪ0�����õ㲻�����Ǳ߽�㣬����ֵΪ255�����õ��Ѿ��Ǳ߽��

void sharpen::TraceEdge(Image &img, int y, int x, int nThrLow,int *pMag){

	const int img_width = img.width;  //ͼ����
	const int img_height = img.height;  //ͼ��߶�

	//��8�������ؽ��в�ѯ
	int xNum[8] = {1,1,0,-1,-1,-1,0,1};
	int yNum[8] = {0,1,1,1,0,-1,-1,-1};

	LONG yy,xx,k;

	for(k=0;k<8;k++)
	{
		yy = y+yNum[k];
		xx = x+xNum[k];

		if(ImageTools::getPixel(img,yy,xx)==128 && pMag[yy*img_width+xx]>=nThrLow )
		{
			//�õ���Ϊ�߽��
			ImageTools::setPixel(img,yy,xx,255);

			//�Ըõ�Ϊ�����ٽ��и���
			TraceEdge(img,yy,xx,nThrLow,pMag);
		}
	}

}


// Canny����
void sharpen::Canny(const Image &img, double sigma, double dRatLow, double dRatHigh){

	const int img_width = img.width;  //ͼ����
	const int img_height = img.height;  //ͼ��߶�

	//x��������ָ��
	int *pGradX;
	pGradX = new int[img_width*img_height];

	//y����
	int *pGradY;
	pGradY = new int[img_width*img_height];

	//�ݶȵķ���
	int *pGradMag;
	pGradMag = new int[img_width*img_height];

	//��ԭͼ��˹�˲�
	Image pGaussSmooth = GaussianSmooth(img,sigma);

	//���㷽�������ݶȵķ���
	Grad(pGaussSmooth,pGradX,pGradY,pGradMag);

	//Ӧ�÷��������
	NonmaxSuppress(pGradMag,pGradX,pGradY,pGaussSmooth);

	//Ӧ��Hysteresis���ҵ����б߽�
	Hysteresis(pGaussSmooth,pGradMag,dRatLow,dRatHigh);

	delete[] pGradX;
	pGradX = NULL;
	delete[] pGradY;
	pGradY = NULL;
	delete[] pGradMag;
	pGradMag = NULL;

	string imgName = tools::fileNameFromTime(PATH,"cannyImg",".jpg");
	ImageTools::SaveImage(pGaussSmooth, imgName.c_str());
}