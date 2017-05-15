#include "sharpen.h"
#include "tools.h"

sharpen::sharpen(void)
{
}


sharpen::~sharpen(void)
{
}


double* sharpen::laplacianGrad(const Image &img){

	const int img_width = img.width;  //图像宽度
	const int img_height = img.height;  //图像高度

	double *grade =new double[img_width*img_height]; 
	// 初始化
	for(int y=0; y<img_height ; y++ )
		for(int x=0 ; x<img_width ; x++ )
		{
			*(grade+y*img_width+x)=0;
		}
		// 设置模板系数
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
		// 这个变量用来表示Laplacian 算子像素值
		int nTmp[3][3];
		// 临时变量
		double dGrad;
		// 模板循环控制变量
		int yy ;
		int xx ;
		// 下面开始利用Laplacian 算子进行计算，为了保证计算所需要的
		// 的数据位于图像数据的内部，下面的两重循环的条件是
		// y<img_height-1 而不是y<img_height，相应的x 方向也是x<img_width-1
		// 而不是x<nWidth
		for(int y=1; y<img_height-1 ; y++ ){
			for(int x=1 ; x<img_width-1 ; x++ )
			{
				dGrad = 0 ;
				// Laplacian 算子需要的各点像素值
				// 模板第一行
				nTmp[0][0] = ImageTools::getPixel(img,y-1,x-1);
				nTmp[0][1] = ImageTools::getPixel(img,y-1,x);
				nTmp[0][2] = ImageTools::getPixel(img,y-1,x+1) ;
				// 模板第二行
				nTmp[1][0] = ImageTools::getPixel(img,y,x-1);
				nTmp[1][1] = ImageTools::getPixel(img,y,x);
				nTmp[1][2] = ImageTools::getPixel(img,y,x+1);
				// 模板第三行
				nTmp[2][0] = ImageTools::getPixel(img,y+1,x-1);
				nTmp[2][1] = ImageTools::getPixel(img,y+1,x);
				nTmp[2][2] = ImageTools::getPixel(img,y+1,x+1);
				// 计算梯度
				for(yy=0; yy<3; yy++)
					for(xx=0; xx<3; xx++)
					{
						dGrad += nTmp[yy][xx] * nWeight[yy][xx] ;
					}
					// 梯度值写入内存
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

	const int img_width = img.width;  //图像宽度
	const int img_height = img.height;  //图像高度

	Image laplacianImg = ImageTools::ImageCopy(img,0,0,img_width,img_height);

	// 设置模板系数
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
			// 计算梯度
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


//  一维高斯分布函数，用于平滑函数中生成的高斯滤波系数
void sharpen::CreatGauss(double sigma, double **pdKernel, int *pnWidowSize)
{

 LONG i;

 //数组中心点
 int nCenter;

 //数组中一点到中心点距离
 double dDis;

 //中间变量
 double dValue;
 double dSum;
 dSum = 0;

 // [-3*sigma,3*sigma] 以内数据，会覆盖绝大部分滤波系数
 *pnWidowSize = static_cast<int>(1+ 2*ceil(3*sigma));

 nCenter = (*pnWidowSize)/2;

 *pdKernel = new double[*pnWidowSize];

 //生成高斯数据
 for(i=0;i<(*pnWidowSize);i++)
 {
  dDis = double(i - nCenter);
  dValue = exp(-(1/2)*dDis*dDis/(sigma*sigma))/(sqrt(2*PI)*sigma);
  (*pdKernel)[i] = dValue;
  dSum+=dValue;

 }
 //归一化
 for(i=0;i<(*pnWidowSize);i++)
 {
  (*pdKernel)[i]/=dSum;
 }

}

//用高斯滤波器平滑原图像
Image sharpen::GaussianSmooth(const Image &img,double sigma){

	const int img_width = img.width;  //图像宽度
	const int img_height = img.height;  //图像高度

	Image resultImg = ImageTools::ImageClone(img);

	LONG x, y;
	LONG i;

	//高斯滤波器长度
	int nWindowSize;

	//窗口长度
	int nLen;

	//一维高斯滤波器
	double *pdKernel;

	//高斯系数与图像数据的点乘
	double dDotMul;

	//滤波系数总和
	double dWeightSum;

	double *pdTemp;
	pdTemp = new double[img_width*img_height];

	//产生一维高斯数据
	CreatGauss(sigma, &pdKernel, &nWindowSize);

	nLen = nWindowSize/2;

	//x方向滤波
	for(y=0;y<img_height;y++)
	{
		for(x=0;x<img_width;x++)
		{
			dDotMul = 0;
			dWeightSum = 0;
			for(i=(-nLen);i<=nLen;i++)
			{
				//判断是否在图像内部
				if((i+x)>=0 && (i+x)<img_width)
				{
					dDotMul+=(double)ImageTools::getPixel(img,y,(i+x)) * pdKernel[nLen+i];
					dWeightSum += pdKernel[nLen+i];
				}
			}
			pdTemp[y*img_width+x] = dDotMul/dWeightSum;
		}
	}

	//y方向滤波
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

// 方向导数,求梯度
void sharpen::Grad(const Image &img,int *pGradX, int *pGradY, int *pMag){

	const int img_width = img.width;  //图像宽度
	const int img_height = img.height;  //图像高度
	LONG y,x;

	//x方向的方向导数
	for(y=1;y<img_height-1;y++)
	{
		for(x=1;x<img_width-1;x++)
		{
			pGradX[y*img_width +x] =ImageTools::getPixel(img,y,x+1)-ImageTools::getPixel(img,y,x-1);
		}
	}

	//y方向方向导数
	for(x=1;x<img_width-1;x++)
	{
		for(y=1;y<img_height-1;y++)
		{
			pGradY[y*img_width +x] = ImageTools::getPixel(img,y+1,x) - ImageTools::getPixel(img,y-1,x);
		}
	}

	//求梯度

	//中间变量
	double dSqt1;
	double dSqt2;

	for(y=0; y<img_height; y++)
	{
		for(x=0; x<img_width; x++)
		{
			//二阶范数求梯度
			dSqt1 = pGradX[y*img_width + x]*pGradX[y*img_width + x];
			dSqt2 = pGradY[y*img_width + x]*pGradY[y*img_width + x];
			pMag[y*img_width+x] = (int)(sqrt(dSqt1+dSqt2)+0.5);
		}
	}
}



//非最大抑制
void sharpen::NonmaxSuppress(int *pMag, int *pGradX, int *pGradY, Image &img){

	const int img_width = img.width;  //图像宽度
	const int img_height = img.height;  //图像高度

	LONG y,x;
	int nPos;

	//梯度分量
	int gx;
	int gy;

	//中间变量
	int g1,g2,g3,g4;
	double weight;
	double dTmp,dTmp1,dTmp2;

	//设置图像边缘为不可能的分界点
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
			//当前点
			nPos = y*img_width + x;

			//如果当前像素梯度幅度为0，则不是边界点
			if(pMag[nPos] == 0)
			{
				ImageTools::setPixel(img,y,x,0);
			}
			else
			{
				//当前点的梯度幅度
				dTmp = pMag[nPos];

				//x,y方向导数
				gx = pGradX[nPos];
				gy = pGradY[nPos];

				//如果方向导数y分量比x分量大，说明导数方向趋向于y分量
				if(abs(gy) > abs(gx))
				{
					//计算插值比例
					weight = abs(gx)/abs(gy);

					g2 = pMag[nPos-img_width];
					g4 = pMag[nPos+img_width];

					//如果x,y两个方向导数的符号相同
					//C 为当前像素，与g1-g4 的位置关系为：
					//g1 g2
					//      C
					//       g4 g3
					if(gx*gy>0)
					{
						g1 = pMag[nPos-img_width-1];
						g3 = pMag[nPos+img_width+1];
					}

					//如果x,y两个方向的方向导数方向相反
					//C是当前像素，与g1-g4的关系为：
					//       g2 g1
					//        C
					//    g3 g4
					else
					{
						g1 = pMag[nPos-img_width+1];
						g3 = pMag[nPos+img_width-1];
					}
				}

				//如果方向导数x分量比y分量大，说明导数的方向趋向于x分量
				else
				{
					//插值比例
					weight = abs(gy)/abs(gx);

					g2 = pMag[nPos+1];
					g4 = pMag[nPos-1];

					//如果x,y两个方向的方向导数符号相同
					//当前像素C与 g1-g4的关系为
					//  g3
					//  g4 C g2
					//       g1
					if(gx * gy > 0)
					{
						g1 = pMag[nPos+img_width+1];
						g3 = pMag[nPos-img_width-1];
					}

					//如果x,y两个方向导数的方向相反
					// C与g1-g4的关系为
					//   g1
					//    g4 C g2
					//     g3
					else
					{
						g1 = pMag[nPos-img_width+1];
						g3 = pMag[nPos+img_width-1];
					}
				}

				//利用 g1-g4 对梯度进行插值
				{
					dTmp1 = weight*g1 + (1-weight)*g2;
					dTmp2 = weight*g3 + (1-weight)*g4;

					//当前像素的梯度是局部的最大值
					//该点可能是边界点
					if(dTmp>=dTmp1 && dTmp>=dTmp2)
					{
						ImageTools::setPixel(img,y,x,128);
					}
					else
					{
						//不可能是边界点
						ImageTools::setPixel(img,y,x,0);
					}
				}
			}
		}
	}

	string imgName = tools::fileNameFromTime(PATH,"cannyNonmaxSuppressImg",".jpg");
	ImageTools::SaveImage(img, imgName.c_str());
}



// 统计pMag的直方图，判定阈值
void sharpen::EstimateThreshold(const Image &img,int *pMag,int *pThrHigh, int *pThrLow,double dRatHigh, double dRatLow){

	const int img_width = img.width;  //图像宽度
	const int img_height = img.height;  //图像高度
	LONG y,x,k;

	//该数组的大小和梯度值的范围有关，如果采用本程序的算法
	//那么梯度的范围不会超过pow(2,10)
	int nHist[256];

	//可能边界数
	int nEdgeNum;

	//最大梯度数
	int nMaxMag;

	int nHighCount;

	nMaxMag = 0;

	//初始化
	for(k=0;k<256;k++)
	{
		nHist[k] = 0;
	}
	//统计直方图,利用直方图计算阈值
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

	//统计经过“非最大值抑制”后有多少像素
	for(k=1;k<256;k++)
	{
		if(nHist[k] != 0)
		{
			nMaxMag = k;
		}

		//梯度为0的点是不可能为边界点的
		//经过non-maximum suppression后有多少像素
		nEdgeNum += nHist[k];

	}

	//梯度比高阈值*pThrHigh 小的像素点总书目
	nHighCount = (int)(dRatHigh * nEdgeNum + 0.5);

	k=1;
	nEdgeNum = nHist[1];

	//计算高阈值
	while((k<(nMaxMag-1)) && (nEdgeNum < nHighCount))
	{
		k++;
		nEdgeNum += nHist[k];
	}

	*pThrHigh = k;

	//低阈值
	*pThrLow = (int)((*pThrHigh) * dRatLow + 0.5);

}

//利用函数寻找边界起点
void sharpen::Hysteresis(Image &img, int *pMag,double dRatLow, double dRatHigh){

	const int img_width = img.width;  //图像宽度
	const int img_height = img.height;  //图像高度

	LONG y,x;

	int nThrHigh,nThrLow;

	int nPos;
	//估计TraceEdge 函数需要的低阈值，以及Hysteresis函数使用的高阈值
	EstimateThreshold(img,pMag,&nThrHigh,&nThrLow,dRatHigh,dRatLow);

	//寻找大于dThrHigh的点，这些点用来当作边界点，
	//然后用TraceEdge函数跟踪该点对应的边界
	for(y=0;y<img_height;y++)
	{
		for(x=0;x<img_width;x++)
		{
			nPos = y*img_width + x;

			//如果该像素是可能的边界点，并且梯度大于高阈值，
			//该像素作为一个边界的起点
			if((ImageTools::getPixel(img,y,x)==128) && (pMag[nPos] >= nThrHigh))
			{
				//设置该点为边界点
				ImageTools::setPixel(img,y,x,255);
				TraceEdge(img,y,x,nThrLow,pMag);
			}

		}
	}

	string imgName = tools::fileNameFromTime(PATH,"cannyTraceEdgeImg",".jpg");
	ImageTools::SaveImage(img, imgName.c_str());

	//其他点已经不可能为边界点
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

//根据Hysteresis 执行的结果，从一个像素点开始搜索，搜索以该像素点为边界起点的一条边界的
//一条边界的所有边界点，函数采用了递归算法
//       从（x,y)坐标出发，进行边界点的跟踪，跟踪只考虑pResult中没有处理并且可能是边界
//  点的像素（=128），像素值为0表明该点不可能是边界点，像素值为255表明该点已经是边界点

void sharpen::TraceEdge(Image &img, int y, int x, int nThrLow,int *pMag){

	const int img_width = img.width;  //图像宽度
	const int img_height = img.height;  //图像高度

	//对8邻域像素进行查询
	int xNum[8] = {1,1,0,-1,-1,-1,0,1};
	int yNum[8] = {0,1,1,1,0,-1,-1,-1};

	LONG yy,xx,k;

	for(k=0;k<8;k++)
	{
		yy = y+yNum[k];
		xx = x+xNum[k];

		if(ImageTools::getPixel(img,yy,xx)==128 && pMag[yy*img_width+xx]>=nThrLow )
		{
			//该点设为边界点
			ImageTools::setPixel(img,yy,xx,255);

			//以该点为中心再进行跟踪
			TraceEdge(img,yy,xx,nThrLow,pMag);
		}
	}

}


// Canny算子
void sharpen::Canny(const Image &img, double sigma, double dRatLow, double dRatHigh){

	const int img_width = img.width;  //图像宽度
	const int img_height = img.height;  //图像高度

	//x方向导数的指针
	int *pGradX;
	pGradX = new int[img_width*img_height];

	//y方向
	int *pGradY;
	pGradY = new int[img_width*img_height];

	//梯度的幅度
	int *pGradMag;
	pGradMag = new int[img_width*img_height];

	//对原图高斯滤波
	Image pGaussSmooth = GaussianSmooth(img,sigma);

	//计算方向导数和梯度的幅度
	Grad(pGaussSmooth,pGradX,pGradY,pGradMag);

	//应用非最大抑制
	NonmaxSuppress(pGradMag,pGradX,pGradY,pGaussSmooth);

	//应用Hysteresis，找到所有边界
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