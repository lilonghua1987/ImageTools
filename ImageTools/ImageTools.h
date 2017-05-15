#pragma once
#include "Image.h"

enum IMAGEFILTER
{
	BOX,BILINEAR,BSPLINE,BICUBIC,CATMULLROM,LANCZOS3
};

class ImageTools
{
public:
	ImageTools(void);
	~ImageTools(void);

	static Pixel<int> getColorPixel(const Image &img, unsigned long row, unsigned long column)
	{
		assert(img.width > column);
		assert(column >= 0);
		assert(img.height > row);
		assert(row >= 0);
		assert(img.channel == 3 || img.channel == 1);

		Pixel<int> value;

		if (img.channel == 1)
		{
			int tmp = getPixel(img,row,column);
			value = Pixel<int>(tmp,tmp,tmp);

			return value;
		}

		if(img.imgType == FIT_RGBF) {

			FIRGBF *bits = (FIRGBF *)FreeImage_GetScanLine(img.data, row);
			value.red = (int)bits[column].red;
			value.green = (int)bits[column].green;
			value.blue = (int)bits[column].blue;
		}else if(img.imgType == FIT_RGBA16)
		{
			FIRGBA16 *bits = (FIRGBA16 *)FreeImage_GetScanLine(img.data, row);
			value.red = (int)bits[column].red;
			value.green = (int)bits[column].green;
			value.blue = (int)bits[column].blue;
			value.alpha = (int)bits[column].alpha;
		}else if(img.imgType == FIT_RGBAF)
		{
			FIRGBAF *bits = (FIRGBAF *)FreeImage_GetScanLine(img.data, row);
			value.red = (int)bits[column].red;
			value.green = (int)bits[column].green;
			value.blue = (int)bits[column].blue;
			value.alpha = (int)bits[column].alpha;
		}

		if((img.imgType == FIT_BITMAP) && (img.channel == 3)) {
			BYTE *bits = (BYTE *)FreeImage_GetScanLine(img.data, row);
			//int bytespp = FreeImage_GetLine(img.data) / FreeImage_GetWidth(img.data);
			bits+=img.channel*column;
			value.red = (int)bits[FI_RGBA_RED];
			value.green = (int)bits[FI_RGBA_GREEN];
			value.blue = (int)bits[FI_RGBA_BLUE];
		}
		return value;
	}

	static int getPixel(const Image &img, unsigned long row, unsigned long column)
	{
		assert(img.width > column);
		assert(column >= 0);
		assert(img.height > row);
		assert(row >= 0);
		assert(img.channel == 1);

		BYTE *bits = (BYTE *)FreeImage_GetScanLine(img.data, row);
		return bits[column];
	};

	static void setPixel(Image &img, unsigned long row, unsigned long column, int value)
	{
		assert(img.width > column);
		assert(column >= 0);
		assert(img.height > row);
		assert(row >= 0);
		assert(img.channel == 1);

		BYTE *bits = (BYTE *)FreeImage_GetScanLine(img.data, row);
		bits[column] = (BYTE)value;
	};

	static void setColorPixel(Image &img, unsigned long row, unsigned long column, Pixel<int> point)
	{
		assert(img.width > column);
		assert(column >= 0);
		assert(img.height > row);
		assert(row >= 0);
		assert(img.channel == 3 || img.channel == 1);

		if (img.channel == 1)
		{
			setPixel(img,row,column,point.red);
			return;
		}

		if(img.imgType == FIT_RGBF) {

			FIRGBF *bits = (FIRGBF *)FreeImage_GetScanLine(img.data, row);
			bits[column].red = (float)point.red;
			bits[column].green = (float)point.green;
			bits[column].blue = (float)point.blue;
		}else if(img.imgType == FIT_RGBA16)
		{
			FIRGBA16 *bits = (FIRGBA16 *)FreeImage_GetScanLine(img.data, row);
			bits[column].red = (WORD)point.red;
			bits[column].green = (WORD)point.green;
			bits[column].blue = (WORD)point.blue;
			bits[column].alpha = (WORD)point.alpha;
		}else if(img.imgType == FIT_RGBAF)
		{
			FIRGBAF *bits = (FIRGBAF *)FreeImage_GetScanLine(img.data, row);
			bits[column].red = (float)point.red;
			bits[column].green = (float)point.green;
			bits[column].blue = (float)point.blue;
			bits[column].alpha = (float)point.alpha;
		}

		if((img.imgType == FIT_BITMAP) && (img.channel == 3)) {
			BYTE *bits = (BYTE *)FreeImage_GetScanLine(img.data, row);
			//int bytespp = FreeImage_GetLine(img.data) / FreeImage_GetWidth(img.data);
			bits+=img.channel*column;
			bits[FI_RGBA_RED] = (BYTE)point.red;
			bits[FI_RGBA_GREEN] = (BYTE)point.green;
			bits[FI_RGBA_BLUE] = (BYTE)point.blue;
		}
	}

	static bool SaveImage(const Image &src, const char imgName[]);

	static Image ImageCopy(const Image& src,long left, long top, long right, long bottom);

	static Image ImageCopy(const Image& src,Point2D<int> top, Point2D<int> bottom);

	static void NormalImage(Image &src, float minValue, float maxValue);

	static Image ImageToGrey(const Image &src);
	static Image ImageClone(const Image &src);

	static Image ImageRotateClassic(const Image &src, double angle);
	static Image ImageRotateEx(const Image &src, double angle,Point2D<double> shiftXY, Point2D<double> rotateCenter, BOOL mask=true);
	static bool ImageFlipHorizontal(Image &src);
	static bool ImageFlipVertical(Image &src);

	static Image ImageRescale(const Image &src, long dstwidth, long dstheight,IMAGEFILTER filterType=IMAGEFILTER::BILINEAR);

	static bool ImageAdjustBrightness(Image &src, double percentage);
	static bool ImageAdjustContrast(Image &src, double percentage);
	static bool ImageInvert(Image &src);

	static Image ImageGetChannel(const Image& src, BYTE channel);
	static bool getHistoGram(const Image& img, unsigned long hist[], IMAGE_CHANNEL channel = FICC_BLACK);

	/*
	*通过使用一个在[0 . . . 255]之间的阈值T来将控制。
	*函数首先将位图转换为8位灰度位图， 然后凡低于T的亮度级都被设置为0， 否则被设置为1
	*/
	static Image Threshold(const Image& src,BYTE t);

	static Image ImageTo32Bits(const Image& src);

};

