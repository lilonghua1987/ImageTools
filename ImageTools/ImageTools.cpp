#include "ImageTools.h"


ImageTools::ImageTools(void)
{
}


ImageTools::~ImageTools(void)
{
}


bool ImageTools::SaveImage(const Image &src, const char imgName[])
{
	IMAGE_FORMAT imgFormat = src.imgFormat;
	BOOL bSuccess = FALSE;

	if(src.data) 
	{
		if(imgFormat != FIF_UNKNOWN ) 
		{
			// check that the plugin has sufficient writing and export capabilities ...
			WORD bpp = FreeImage_GetBPP(src.data);
			if(FreeImage_FIFSupportsWriting(imgFormat) && FreeImage_FIFSupportsExportBPP(imgFormat, bpp)) 
			{
				// ok, we can save the file
				bSuccess = FreeImage_Save(imgFormat,src.data,imgName);
				// unless an abnormal bug, we are done !
			}
		}
	}
	if(!bSuccess) std::cout<<"Save image fail!"<<std::endl;
	return bSuccess?true:false;
	/*FreeImage_Save(src.imgFormat,src.data,imgName);*/
}



Image ImageTools::ImageCopy(const Image &src,long left, long top, long right, long bottom)
{
	Image img;

	img.data = FreeImage_Copy(src.data,left,top,right,bottom);
	img.channel = src.channel;
	img.imgFormat = src.imgFormat;
	img.height = bottom - top;
	img.imgExt = src.imgExt;
	img.imgType = src.imgType;
	img.width = right - left;
	img.bitsSize = img.width*img.height*img.channel;

	return img;
}


Image ImageTools::ImageCopy(const Image &src,Point2D<int> top, Point2D<int> bottom)
{
	return ImageCopy(src,top.x,top.y,bottom.x,bottom.y);
}


void ImageTools::NormalImage(Image &src, float minValue, float maxValue)
{
	if (minValue >= maxValue) return;

	Pixel<float> minPoint,maxPoint;
	minPoint = maxPoint = src.getPixel<float>(0,0);

	for (ulong i = 0; i < src.height; i++)
	{
		for (ulong j = 0; j < src.width; j++)
		{
			Pixel<float> p =src.getPixel<float>(i, j);

			if (p.red < minPoint.red)
			{
				minPoint.red = p.red;
			}

			if (p.red > maxPoint.red)
			{
				maxPoint.red = p.red;
			}

			if (p.green < minPoint.green)
			{
				minPoint.green = p.green;
			}

			if (p.green > maxPoint.green)
			{
				maxPoint.green = p.green;
			}

			if (p.blue < minPoint.blue)
			{
				minPoint.blue = p.blue;
			}

			if (p.blue > maxPoint.blue)
			{
				maxPoint.blue = p.blue;
			}

			if (p.alpha < minPoint.alpha)
			{
				minPoint.alpha = p.alpha;
			}

			if (p.alpha > maxPoint.alpha)
			{
				maxPoint.alpha = p.alpha;
			}
		}
	}

	//
	float boundR = maxPoint.red - minPoint.red;
	float boundG = maxPoint.green - minPoint.green;
	float boundB = maxPoint.blue - minPoint.blue;
	float boundA = maxPoint.alpha - minPoint.alpha;
	float bound = maxValue - minValue;

	Pixel<float> pRate;

	if (src.channel == 1)
	{
		if ((maxPoint.red == maxValue) && (minPoint.red == minValue))
		{
			return;
		}
		if (boundR <= 0) return;
		pRate.red = bound/boundR;
	}

	if (src.channel == 3)
	{
		/*if ((maxPoint == Pixel<int>(maxValue,maxValue,maxValue)) && (minPoint == Pixel<int>(minValue,minValue,minValue)))
		{
		return;
		}*/
		if (boundR <= 0 || boundG <= 0 || boundB <= 0) return;
		pRate.red = bound/boundR;
		pRate.green = bound/boundG;
		pRate.blue = bound/boundB;
	}
	if (src.channel == 4)
	{
		if (boundR <= 0 || boundG <= 0 || boundB <= 0 || boundA <= 0) return;
		pRate.red = bound/boundR;
		pRate.green = bound/boundG;
		pRate.blue = bound/boundB;
		pRate.alpha = bound/boundA;
	}

	/*Log log("image.txt");*/
	for (ulong i = 0; i < src.height; i++)
	{
		for (ulong j = 0; j < src.width; j++)
		{
			Pixel<float> pixel = (src.getPixel<float>(i,j)-minPoint)*pRate+Pixel<float>(minValue,minValue,minValue,minValue);
			src.setPixel(i,j,pixel);
		}
		/*log.systemLog("\n");*/
	}
}


Image ImageTools::ImageToGrey(const Image &src)
{
	Image img;
	img.data = FreeImage_ConvertToGreyscale(src.data);
	img.channel = 1;
	img.imgType = src.imgType;
	img.width = src.width;
	img.height = src.height;
	img.imgFormat = src.imgFormat;
	img.imgExt = src.imgExt;
	img.bitsSize = img.width*img.height*img.channel;
	return img;
}


Image ImageTools::ImageClone(const Image &src)
{
	Image img;
	img.data = FreeImage_Clone(src.data);
	img.channel = src.channel;
	img.imgType = src.imgType;
	img.width = src.width;
	img.height = src.height;
	img.imgFormat = src.imgFormat;
	img.imgExt = src.imgExt;
	img.bitsSize = img.width*img.height*img.channel;
	return img;
}



/*
该函数用三种剪切(shear)来旋转一个1位、8位灰度位图或一个24位、32位
彩色图象，旋转角度由angle参数以度为单位指定。旋转是围绕图象中心进行的。
被旋转的图象保持源图象的大小和长宽比(通常目标图象比较大)，所以应该在
将一幅图象旋转90°、180°或270°时使用该函数。
*/
Image ImageTools::ImageRotateClassic(const Image& src, double angle)
{
	Image img;
	img.data = FreeImage_RotateClassic(src.data,angle);
	img.channel = src.channel;
	img.imgType = src.imgType;
	img.width = FreeImage_GetWidth(img.data);
	img.height = FreeImage_GetHeight(img.data);
	img.imgFormat = src.imgFormat;
	img.imgExt = src.imgExt;
	img.bitsSize = img.width*img.height*img.channel;
	return img;
}



/*
旋转角度由angle参数以度为单位指定。水平和垂直平移(以像素为单
位)由xshift参数和yshift参数指定。旋转围绕着由xorigin和yorigin(同样以
像素为单位)指定的中心进行。当usemask被设为TRUE时，与图象无关的部分
被设为黑色，否则用反射技术来填充无关的像素。
*/
Image ImageTools::ImageRotateEx(const Image &src, double angle,Point2D<double> shiftXY, Point2D<double> rotateCenter, BOOL mask)
{
	Image img;
	img.data = FreeImage_RotateEx(src.data,angle,shiftXY.x,shiftXY.y,rotateCenter.x,rotateCenter.y,mask);
	img.channel = src.channel;
	img.imgType = src.imgType;
	img.width = FreeImage_GetWidth(img.data);
	img.height = FreeImage_GetHeight(img.data);
	img.imgFormat = src.imgFormat;
	img.imgExt = src.imgExt;
	img.bitsSize = img.width*img.height*img.channel;
	return img;
}



/*
沿垂直轴将输入img水平翻转。
*/
bool ImageTools::ImageFlipHorizontal(Image &src)
{
	return FreeImage_FlipHorizontal(src.data) ? true:false;
}


/*
沿水平轴将输入img垂直翻转
*/
bool ImageTools::ImageFlipVertical(Image &src)
{
	return FreeImage_FlipVertical(src.data) ? true:false;
}



/*
RBOX 箱形(Box),脉冲,傅立叶窗,1阶(常量) B样条
BILINEAR 双线性(Bilinear)滤镜
BSPLINE 4阶(立方)B样条
BICUBIC Mitchell-Netravali双参数立方(two-param cubic)滤
CATMULLROM Catmull-Rom样条,Overhauser样条
LANCZOS3 Lanczos-windowed sinc滤镜
*/
Image ImageTools::ImageRescale(const Image &src, long dstwidth, long dstheight,IMAGEFILTER filterType)
{
	FREE_IMAGE_FILTER filter;
	switch (filterType)
	{
	case BOX:
		filter = FILTER_BOX;
		break;
	case BILINEAR:
		filter = FILTER_BILINEAR;
		break;
	case BSPLINE:
		filter = FILTER_BSPLINE;
		break;
	case BICUBIC:
		filter = FILTER_BICUBIC;
		break;
	case CATMULLROM:
		filter = FILTER_CATMULLROM;
		break;
	case LANCZOS3:
		filter = FILTER_LANCZOS3;
		break;
	default:
		filter = FILTER_BOX;
		break;
	}

	Image img;
	img.data = FreeImage_Rescale(src.data,dstwidth,dstheight,filter);
	img.channel = src.channel;
	img.imgType = src.imgType;
	img.width = dstwidth;
	img.height = dstheight;
	img.imgFormat = src.imgFormat;
	img.imgExt = src.imgExt;
	img.bitsSize = img.width*img.height*img.channel;
	return img;
}



/*
将一幅8位、24位或32位图象的亮度作一定数量的调整，该数量由percentage参
数给定，其值为[-100..100]之间的一个数。0值意味着无变化，小于0使图象变深
而大于0使图象变浅。
如果调用成功，函数返回TRUE，否则返回FALSE(例如当无法处理源data的
位深度时)。
*/
bool ImageTools::ImageAdjustBrightness(Image &src, double percentage)
{
	return FreeImage_AdjustBrightness(src.data,percentage) ? true:false;
}


/*
将一幅8位、 24位或32位图象的对比度作一定数量的调整， 该数量由percentage参
数给定， 其值为[-100..100]之间的一个数。 0值意味着无变化， 小于0使图象对比
度减小而大于0使图象对比度增加
*/
bool ImageTools::ImageAdjustContrast(Image &src, double percentage)
{
	return FreeImage_AdjustContrast(src.data,percentage) ? true : false;
}



/*
反转(Invert)每个象素数据。
*/
bool ImageTools::ImageInvert(Image &src)
{
	return FreeImage_Invert(src.data) ? true:false;
}


Image ImageTools::ImageGetChannel(const Image& src, BYTE channel)
{
	assert(src.depth == 24 || src.depth == 32);
	assert(channel >= 0 && channel <= 4);
	Image img;
	img.channel = src.channel;
	img.imgType = src.imgType;
	img.width = src.width;
	img.height = src.height;
	img.depth = src.depth;
	img.imgFormat = src.imgFormat;
	img.imgExt = src.imgExt;
	img.bitsSize = img.width*img.height*img.channel;

	switch (channel)
	{
	case FICC_RGB:
		{
			img.data = FreeImage_GetChannel(src.data,FICC_RGB);
			break;
		}
	case FICC_RED:
		{
			img.data = FreeImage_GetChannel(src.data,FICC_RED);
			img.depth = 8;
			img.channel = 1;
			break;
		}
	case FICC_GREEN:
		{
			img.data = FreeImage_GetChannel(src.data,FICC_GREEN);
			img.depth = 8;
			img.channel = 1;
			break;
		}
	case FICC_BLUE:
		{
			img.data = FreeImage_GetChannel(src.data,FICC_BLUE);
			img.depth = 8;
			img.channel = 1;
			break;
		}
	case FICC_ALPHA:
		{
			img.data = FreeImage_GetChannel(src.data,FICC_ALPHA);
			img.depth = 8;
			img.channel = 1;
			break;
		}
	default:
		break;
	}

	return img;
}


bool ImageTools::getHistoGram(const Image& img, unsigned long hist[], IMAGE_CHANNEL channel /* = FICC_BLACK */)
{
	BOOL flage = FreeImage_GetHistogram(img.data,hist,channel);
	return flage ? true : false;
}


Image ImageTools::Threshold(const Image& src,BYTE t)
{
	assert(t>=0 && t <=255);

	Image img;
	img.channel = 1;
	img.imgType = src.imgType;
	img.width = src.width;
	img.height = src.height;
	img.depth = 8;
	img.imgFormat = src.imgFormat;
	img.imgExt = src.imgExt;
	img.data = FreeImage_ConvertTo8Bits(FreeImage_Threshold(src.data,t));
	img.bitsSize = img.width*img.height*img.channel;

	return img;
}


Image ImageTools::ImageTo32Bits(const Image& src)
{
	Image img;
	img.channel = 4;
	img.imgType = src.imgType;
	img.width = src.width;
	img.height = src.height;
	img.depth = 32;
	img.imgFormat = src.imgFormat;
	img.imgExt = src.imgExt;
	img.data = FreeImage_ConvertTo32Bits(src.data);	
	img.bitsSize = img.width*img.height*img.channel;
	img.redMask = FreeImage_GetRedMask(img.data);
	img.greenMask = FreeImage_GetGreenMask(img.data);
	img.blueMask = FreeImage_GetBlueMask(img.data);
	img.alphaMask=FreeImage_GetTransparencyCount(img.data);

	return img;
}