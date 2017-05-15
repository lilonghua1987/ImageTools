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
�ú��������ּ���(shear)����תһ��1λ��8λ�Ҷ�λͼ��һ��24λ��32λ
��ɫͼ����ת�Ƕ���angle�����Զ�Ϊ��λָ������ת��Χ��ͼ�����Ľ��еġ�
����ת��ͼ�󱣳�Դͼ��Ĵ�С�ͳ����(ͨ��Ŀ��ͼ��Ƚϴ�)������Ӧ����
��һ��ͼ����ת90�㡢180���270��ʱʹ�øú�����
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
��ת�Ƕ���angle�����Զ�Ϊ��λָ����ˮƽ�ʹ�ֱƽ��(������Ϊ��
λ)��xshift������yshift����ָ������תΧ������xorigin��yorigin(ͬ����
����Ϊ��λ)ָ�������Ľ��С���usemask����ΪTRUEʱ����ͼ���޹صĲ���
����Ϊ��ɫ�������÷��似��������޹ص����ء�
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
�ش�ֱ�Ὣ����imgˮƽ��ת��
*/
bool ImageTools::ImageFlipHorizontal(Image &src)
{
	return FreeImage_FlipHorizontal(src.data) ? true:false;
}


/*
��ˮƽ�Ὣ����img��ֱ��ת
*/
bool ImageTools::ImageFlipVertical(Image &src)
{
	return FreeImage_FlipVertical(src.data) ? true:false;
}



/*
RBOX ����(Box),����,����Ҷ��,1��(����) B����
BILINEAR ˫����(Bilinear)�˾�
BSPLINE 4��(����)B����
BICUBIC Mitchell-Netravali˫��������(two-param cubic)��
CATMULLROM Catmull-Rom����,Overhauser����
LANCZOS3 Lanczos-windowed sinc�˾�
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
��һ��8λ��24λ��32λͼ���������һ�������ĵ�������������percentage��
����������ֵΪ[-100..100]֮���һ������0ֵ��ζ���ޱ仯��С��0ʹͼ�����
������0ʹͼ���ǳ��
������óɹ�����������TRUE�����򷵻�FALSE(���統�޷�����Դdata��
λ���ʱ)��
*/
bool ImageTools::ImageAdjustBrightness(Image &src, double percentage)
{
	return FreeImage_AdjustBrightness(src.data,percentage) ? true:false;
}


/*
��һ��8λ�� 24λ��32λͼ��ĶԱȶ���һ�������ĵ����� ��������percentage��
�������� ��ֵΪ[-100..100]֮���һ������ 0ֵ��ζ���ޱ仯�� С��0ʹͼ��Ա�
�ȼ�С������0ʹͼ��Աȶ�����
*/
bool ImageTools::ImageAdjustContrast(Image &src, double percentage)
{
	return FreeImage_AdjustContrast(src.data,percentage) ? true : false;
}



/*
��ת(Invert)ÿ���������ݡ�
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