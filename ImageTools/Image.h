#pragma once

#include "FreeImage.h"
#include <iostream>
#include <assert.h>
#include "CoreStruct.h"

#pragma comment(lib, "lib/FreeImage.lib")

/*
图像类的封装
*/

typedef FREE_IMAGE_TYPE IMAGE_TYPE;
typedef FREE_IMAGE_FORMAT IMAGE_FORMAT;
typedef FREE_IMAGE_COLOR_CHANNEL IMAGE_CHANNEL;

class HeaderInfo
{
public:
	friend class Image;
	HeaderInfo()
		:make(std::string())
		, model(std::string())
		, focal(-1)
	{
	}
	~HeaderInfo(){};

private:
	std::string make;
	std::string model;
	float focal;

public:
	std::string getMake()const
	{
		return make;
	}

	std::string getModel()const
	{
		return model;
	}

	float getFocal()const
	{
		return focal;
	}

	bool hasInfo()const
	{
		return !make.empty() && !model.empty() && (focal > 0);
	}
};

class Image
{
public:

	Image(void)
		: data(NULL)
		, width(0)
		, height(0)
		, imgType(FIT_BITMAP)
		, imgFormat(FIF_BMP)
		, imgExt(std::string(".bmp"))
		, channel(0)
		, depth(0)
		, redMask(0)
		, greenMask(0)
		, blueMask(0)
		, alphaMask(0)
		, bitsSize(0)
	{
	}

	Image(unsigned long width, unsigned long height, IMAGE_TYPE imgType=FIT_BITMAP, unsigned char channel=1)
		: data(NULL)
		, width(width)
		, height(height)
		, imgExt(std::string(".bmp"))
		, imgType(imgType)
		, imgFormat(FIF_BMP)
		, channel(channel)
		, depth(0)
	{
		if (!createImage(width, height, imgType, channel))
		{
			return;
		} 
	}

	Image(const Point2D<unsigned long>& imSize, IMAGE_TYPE imgType=FIT_BITMAP, unsigned char channel=1)
		:data(NULL)
		,width(imSize.x)
		,height(imSize.y)
		, imgExt(std::string(".bmp"))
		, imgType(imgType)
		, imgFormat(FIF_BMP)
		, channel(channel)
		, depth(0)
	{
		if (!createImage(width, height, imgType, channel))
		{
			return;
		} 
	}

	Image(const Image &image)
	{
		/*此处data为动态分配的内存空间，拷贝时必须动态开辟空间
		否则就会出现两次释放同一内存空间
		*/
		this->data = FreeImage_Clone(image.data);
		this->channel = image.channel;
		this->depth = image.depth;
		this->imgType = image.imgType;
		this->width = image.width;
		this->height = image.height;
		this->imgFormat = image.imgFormat;
		this->imgExt = image.imgExt;
		fillInfo(data);
	}

	Image(const char srcPath[])
	{
		data = loadImage(srcPath);
	}

	Image(const std::string& srcPath)
	{
		data = loadImage(srcPath.c_str());
	}

	~Image(void)
	{
		if (data)
			FreeImage_Unload(data);
		//imgExt 不是new出来的 不需要释放
		//delete只能删除堆上内存; 栈上的内存, 是不需要手动来释放的
	}

public:
	bool createImage(unsigned long width, unsigned long height, IMAGE_TYPE imgType = FIT_BITMAP, unsigned char channel = 1)
	{
		this->channel = channel;
		if (!initCreat(imgType))
		{
			return false;
		}

		this->data = FreeImage_AllocateT(imgType, width, height,depth);
		this->width = width;
		this->height = height;
		fillInfo(data);

		return true;
	}

	Image& operator = (const Image &image)
	{
		if(this == &image)
			return *this;     //防止出现自赋值

		if (data)
			FreeImage_Unload(data);
		this->data = FreeImage_Clone(image.data);
		this->channel = image.channel;
		this->depth = image.depth;
		this->imgType = image.imgType;
		this->width = image.width;
		this->height = image.height;
		this->imgFormat = image.imgFormat;
		this->imgExt = image.imgExt;
		/*this->filePath = image.filePath;*/
		fillInfo(data);

		return *this;
	}

	template<typename T>
	Pixel<T> getPixel(unsigned long row, unsigned long column)const
	{
		assert(width > column);
		assert(column >= 0);
		assert(height > row);
		assert(row >= 0);

		Pixel<T> value;

		switch (imgType)
		{
		case FIT_UNKNOWN:
			break;
		case FIT_BITMAP:
			{
				BYTE *bits = (BYTE *)FreeImage_GetScanLine(data, row);

				if (channel == 1)
				{
					value.red = (T)bits[column];
				} 
				else
				{
					bits+=channel*column;
					value.red = (T)bits[FI_RGBA_RED];
					value.green = (T)bits[FI_RGBA_GREEN];
					value.blue = (T)bits[FI_RGBA_BLUE];
					if(channel == 4) value.alpha = (T)bits[FI_RGBA_ALPHA];
				}
			}
			break;
		case FIT_UINT16:
			{
				unsigned short *bits = (unsigned short *)FreeImage_GetScanLine(data, row);
				value.red = (T)bits[column];
			}
			break;
		case FIT_INT16:
			{
				short *bits = (short *)FreeImage_GetScanLine(data, row);
				value.red = (T)bits[column];
			}
			break;
		case FIT_UINT32:
			{
				unsigned long *bits = (unsigned long *)FreeImage_GetScanLine(data, row);
				value.red = (T)bits[column];
			}
			break;
		case FIT_INT32:
			{
				long *bits = (long *)FreeImage_GetScanLine(data, row);
				value.red = (T)bits[column];
			}
			break;
		case FIT_FLOAT:
			{
				float *bits = (float *)FreeImage_GetScanLine(data, row);
				value.red = (T)bits[column];
			}
			break;
		case FIT_DOUBLE:
			{
				double *bits = (double *)FreeImage_GetScanLine(data, row);
				value.red = (T)bits[column];
			}
			break;
		case FIT_COMPLEX:
			{
				FICOMPLEX *bits = (FICOMPLEX *)FreeImage_GetScanLine(data, row);
				value.red = (T)bits[column].r;
				value.green = (T)bits[column].i;
			}
			break;
		case FIT_RGB16:
			{
				FIRGB16 *bits = (FIRGB16 *)FreeImage_GetScanLine(data, row);
				value.red = (T)bits[column].red;
				value.green = (T)bits[column].green;
				value.blue = (T)bits[column].blue;
			}
			break;
		case FIT_RGBA16:
			{
				FIRGBA16 *bits = (FIRGBA16 *)FreeImage_GetScanLine(data, row);
				value.red = (T)bits[column].red;
				value.green = (T)bits[column].green;
				value.blue = (T)bits[column].blue;
				value.alpha = (T)bits[column].alpha;
			}
			break;
		case FIT_RGBF:
			{
				FIRGBF *bits = (FIRGBF *)FreeImage_GetScanLine(data, row);
				value.red = (T)bits[column].red;
				value.green = (T)bits[column].green;
				value.blue = (T)bits[column].blue;
			}
			break;
		case FIT_RGBAF:
			{
				FIRGBAF *bits = (FIRGBAF *)FreeImage_GetScanLine(data, row);
				value.red = (T)bits[column].red;
				value.green = (T)bits[column].green;
				value.blue = (T)bits[column].blue;
				value.alpha = (T)bits[column].alpha;
			}
			break;
		}
		return value;
	}


	template<typename T>
	Pixel<T> getPixel(const Point2D<unsigned long>& point)const
	{
		return getPixel(point.y,point.x);
	}

	template<typename T>
	void setPixel(unsigned long row, unsigned long column, const Pixel<T> &point)
	{
		assert(width > column);
		assert(column >= 0);
		assert(height > row);
		assert(row >= 0);
		assert(channel > 0);

		switch (imgType)
		{
		case FIT_UNKNOWN:
			break;
		case FIT_BITMAP:
			{
				BYTE *bits = (BYTE *)FreeImage_GetScanLine(data, row);

				if (channel == 1)
				{
					bits[column] = (BYTE)point.red;
				} 
				else
				{
					bits+=channel*column;
					bits[FI_RGBA_RED] = (BYTE)point.red;
					bits[FI_RGBA_GREEN] = (BYTE)point.green;
					bits[FI_RGBA_BLUE] = (BYTE)point.blue;
					if(channel==4) bits[FI_RGBA_ALPHA] = (BYTE)point.alpha;
				}
			}
			break;
		case FIT_UINT16:
			{
				unsigned short *bits = (unsigned short *)FreeImage_GetScanLine(data, row);
				bits[column] = (unsigned short)point.red;
			}
			break;
		case FIT_INT16:
			{
				short *bits = (short *)FreeImage_GetScanLine(data, row);
				bits[column] = (short)point.red;
			}
			break;
		case FIT_UINT32:
			{
				unsigned long *bits = (unsigned long *)FreeImage_GetScanLine(data, row);
				bits[column] = (unsigned long)point.red;
			}
			break;
		case FIT_INT32:
			{
				long *bits = (long *)FreeImage_GetScanLine(data, row);
				bits[column] = (long)point.red;
			}
			break;
		case FIT_FLOAT:
			{
				float *bits = (float *)FreeImage_GetScanLine(data, row);
				bits[column] = (float)point.red;
			}
			break;
		case FIT_DOUBLE:
			{
				double *bits = (double *)FreeImage_GetScanLine(data, row);
				bits[column] = (double)point.red;
			}
			break;
		case FIT_COMPLEX:
			{
				FICOMPLEX *bits = (FICOMPLEX *)FreeImage_GetScanLine(data, row);
				bits[column].r = (double)point.red;
				bits[column].i = (double)point.green;
			}
			break;
		case FIT_RGB16:
			{
				FIRGB16 *bits = (FIRGB16 *)FreeImage_GetScanLine(data, row);
				bits[column].red = (unsigned short)point.red;
				bits[column].green = (unsigned short)point.green;
				bits[column].blue = (unsigned short)point.blue;
			}
			break;
		case FIT_RGBA16:
			{
				FIRGBA16 *bits = (FIRGBA16 *)FreeImage_GetScanLine(data, row);
				bits[column].red = (unsigned short)point.red;
				bits[column].green = (unsigned short)point.green;
				bits[column].blue = (unsigned short)point.blue;
				bits[column].alpha = (unsigned short)point.alpha;
			}
			break;
		case FIT_RGBF:
			{
				FIRGBF *bits = (FIRGBF *)FreeImage_GetScanLine(data, row);
				bits[column].red = (float)point.red;
				bits[column].green = (float)point.green;
				bits[column].blue = (float)point.blue;
			}
			break;
		case FIT_RGBAF:
			{
				FIRGBAF *bits = (FIRGBAF *)FreeImage_GetScanLine(data, row);
				bits[column].red = (float)point.red;
				bits[column].green = (float)point.green;
				bits[column].blue = (float)point.blue;
				bits[column].alpha = (float)point.alpha;
			}
			break;
		}
	}

	template<typename T>
	void setPixel(const Point2D<unsigned long>& point, const Pixel<T> &pixel)
	{
		setPixel(point.y,point.x,pixel);
	}

	//返回一个指向位图数据位中给定扫描位首的引用
	template<typename T>
	T& at(unsigned long row, unsigned long column)const
	{
		assert(width > column);
		assert(column >= 0);
		assert(height > row);
		assert(row >= 0);
		//assert(channel == 1);

		T *bits = (T *)FreeImage_GetScanLine(data, row);
		return bits[column];
	};

	template<typename T>
	T& at(const Point2D<unsigned long>& point)const
	{
		return at(point.y,point.x);
	}

	template<typename T>
	T& atVal(unsigned long row, unsigned long column, BYTE channel)const
	{
		assert(width > column);
		assert(column >= 0);
		assert(height > row);
		assert(row >= 0);
		assert(channel < this->channel && channel >= 0);

		T *bits = (T *)FreeImage_GetScanLine(data, row);
		return bits[column*this->channel + channel];
	}

	//返回对应行的首地址
	template<typename T>
	T* operator [](unsigned long row)const
	{
		assert(column >= 0);
		assert(height > row);
		assert(row >= 0);

		T *bits = (T *)FreeImage_GetScanLine(data, row);
		return bits;
	};

	template<class T>
	T& operator ()(unsigned long row, unsigned long column)const
	{
		assert(width > column);
		assert(column >= 0);
		assert(height > row);
		assert(row >= 0);

		T *bits = (T *)FreeImage_GetScanLine(data, row);
		return bits[column];
	}

	//只适用于位图
	template<typename T>
	Pixel<T> get(unsigned long row, unsigned long column)const
	{
		assert(width > column);
		assert(column >= 0);
		assert(height > row);
		assert(row >= 0);
		assert(imgType == FIT_BITMAP);

		Pixel<T> value;

		BYTE *bits = (BYTE *)FreeImage_GetScanLine(data, row);

		if (channel == 1)
		{
			value.red = (T)bits[column];
		} 
		else
		{
			bits+=channel*column;
			value.red = (T)bits[FI_RGBA_RED];
			value.green = (T)bits[FI_RGBA_GREEN];
			value.blue = (T)bits[FI_RGBA_BLUE];
			if(channel == 4) value.alpha = (T)bits[FI_RGBA_ALPHA];
		}
		return value;
	}

	//只适用于位图
	template<typename T>
	Pixel<T> get(const Point2D<unsigned long>& point)const
	{
		return get<T>(point.y,point.x);
	}

	//只适用于位图
	template<typename T>
	void set(unsigned long row, unsigned long column, const Pixel<T> &point)
	{
		assert(width > column);
		assert(column >= 0);
		assert(height > row);
		assert(row >= 0);
		assert(channel > 0);
		assert(imgType == FIT_BITMAP);

		BYTE *bits = (BYTE *)FreeImage_GetScanLine(data, row);

		if (channel == 1)
		{
			bits[column] = (BYTE)point.red;
		} 
		else
		{
			bits+=channel*column;
			bits[FI_RGBA_RED] = (BYTE)point.red;
			bits[FI_RGBA_GREEN] = (BYTE)point.green;
			bits[FI_RGBA_BLUE] = (BYTE)point.blue;
			if(channel==4) bits[FI_RGBA_ALPHA] = (BYTE)point.alpha;
		}
	}

	//只适用于位图
	template<typename T>
	void set(const Point2D<unsigned long>& point, const Pixel<T> &pixel)
	{
		set(point.y,point.x,pixel);
	}

	//返回一个指向位图的数据位的指针
	template<typename T>
	T* getPtr()const
	{
		return (T*)FreeImage_GetBits(data);
	};

	bool bound(const Point2D<int>& p)const
	{
		if (p.x < 0 || p.y < 0) return false;

		if (p.x >= (int)width || p.y >= (int)height) return false;

		return true;
	}

	//
	template<typename T>
	T maxVal()const
	{
		assert(channel == 1);

		T maxV = (std::numeric_limits<T>::min)();
		for (ulong i = 0; i < height; i++)
		{
			for (ulong j = 0; j < width; j++)
			{
				if (at<T>(i,j) > maxV) maxV = at<T>(i,j);
			}
		}

		return maxV;
	};


	//
	template<typename T>
	T minVal()const
	{
		assert(channel == 1);

		T minV = (std::numeric_limits<T>::max)();
		for (ulong i = 0; i < height; i++)
		{
			for (ulong j = 0; j < width; j++)
			{
				if (at<T>(i,j) < minV) minV = at<T>(i,j);
			}
		}

		return minV;
	};


	//
	bool save(const std::string& fileName)
	{
		assert(!fileName.empty());
		assert(data != nullptr);
		std::string imgeName(fileName);
		if(imgeName.find_last_of('.') == std::string::npos) imgeName += imgExt;

		imgFormat = FreeImage_GetFileType(imgeName.c_str(), 0);

		if(imgFormat == FIF_UNKNOWN) imgFormat = FreeImage_GetFIFFromFilename(imgeName.c_str());

		BOOL bSuccess = false;

		if(data && (imgFormat != FIF_UNKNOWN )) 
		{
			// check that the plugin has sufficient writing and export capabilities ...
			WORD bpp = FreeImage_GetBPP(data);
			if(FreeImage_FIFSupportsWriting(imgFormat) && FreeImage_FIFSupportsExportBPP(imgFormat, bpp)) 
			{
				bSuccess = FreeImage_Save(imgFormat,data,imgeName.c_str());
			}
		}
		if(!bSuccess) std::cout<<"Save image fail!"<<std::endl;
		return bSuccess?true:false;
	}

	//geometry
	void drawLine(const Point2D<int>& begineP, const Point2D<int>& endP,const Pixel<BYTE>& color)
	{
		int x = begineP.x;
		int y = begineP.y;
		int dx = abs(endP.x - x);
		int dy = abs(endP.y - y);
		int s1 = endP.x > x ? 1 : -1;
		int s2 = endP.y > y ? 1 : -1;

		bool interchange = false;	// 默认不互换 dx、dy
		if (dy > dx)				// 当斜率大于 1 时，dx、dy 互换
		{
			int temp = dx;
			dx = dy;
			dy = temp;
			interchange = true;
		}

		int p = 2 * dy - dx;
		for(int i = 0; i < dx; i++)
		{
			setPixel(y, x, color);
			if (p >= 0)
			{
				if (!interchange)		// 当斜率 < 1 时，选取上下象素点
					y += s2;
				else					// 当斜率 > 1 时，选取左右象素点
					x += s1;
				p -= 2 * dx;
			}
			if (!interchange)
				x += s1;				// 当斜率 < 1 时，选取 x 为步长
			else
				y += s2;				// 当斜率 > 1 时，选取 y 为步长
			p += 2 * dy;
	   }
		setPixel(endP, color);
	}

private:
	FIBITMAP* loadImage(const char srcPath[])
	{
		data = NULL;

		imgType = FIT_BITMAP;

		/*filePath = std::string(srcPath);*/

		imgFormat = FIF_UNKNOWN;
		// check the file signature and deduce its format
		// (the second argument is currently not used by FreeImage)
		imgFormat = FreeImage_GetFileType(srcPath, 0);

		if(imgFormat == FIF_UNKNOWN) {
			// no signature ?
			// try to guess the file format from the file extension
			imgFormat = FreeImage_GetFIFFromFilename(srcPath);
		}

		// check that the plugin has reading capabilities ...
		if((imgFormat != FIF_UNKNOWN) && FreeImage_FIFSupportsReading(imgFormat)) {
			// ok, let's load the file
			data = FreeImage_Load(imgFormat, srcPath, 0);
			imgType =  FreeImage_GetImageType(data);  //图像类型
			width = FreeImage_GetWidth(data);
			height = FreeImage_GetHeight(data);
			// unless a bad file format, we are done !
		}

		if (data)
		{
			switch (imgFormat)
			{
			case FIF_BMP:
				imgExt = ".bmp";
				break;
			case FIF_ICO:
				imgExt = ".ico";
				break;
			case FIF_JPEG:
				imgExt = ".jpg";
				break;
			case FIF_PNG:
				imgExt = ".png";
				break;
			case FIF_PPM:
				imgExt = ".ppm";
				break;
			case FIF_PGM:
				imgExt = ".pgm";
				break;
			case FIF_PBM:
				imgExt = ".pbm";
				break;
			case FIF_TIFF:
				imgExt = ".tiff";
				break;
			case FIF_PSD:
				imgExt = ".psd";
				break;
			case FIF_GIF:
				imgExt = ".gif";
				break;
			case FIF_HDR:
				imgExt = ".hdr";
				break;
			case FIF_RAW:
				imgExt = ".raw";
				break;
			default:
				imgExt = ".jpg";
				break;
			}

			//
			channel = (unsigned char)(FreeImage_GetLine(data) / width);
			if (!initCreat(imgType))
			{
				//channel = FreeImage_GetLine(data) / width;
				depth = FreeImage_GetBPP(data);
			}

			fillInfo(data);
		}
		return data;
	}

	bool initCreat(IMAGE_TYPE imgType)
	{
		bool flage = false;
		//this->imgType = imgType;
		switch (imgType)
		{
		case FIT_UNKNOWN:
			flage = false;
			break;
		case FIT_BITMAP:
			depth = channel*8;
			flage = true;
			break;
		case FIT_UINT16:
			channel = 1;
			depth = channel*16;
			flage = true;
			break;
		case FIT_INT16:
			depth = channel*16;
			channel = 1;
			flage = true;
			break;
		case FIT_UINT32:
			channel = 1;
			depth = channel*32;
			flage = true;
			break;
		case FIT_INT32:
			channel = 1;
			depth = channel*32;
			flage = true;
			break;
		case FIT_FLOAT:
			channel = 1;
			depth = channel*32;			
			flage = true;
			break;
		case FIT_DOUBLE:
			channel = 1;
			depth = channel*64;
			flage = true;
			break;
		case FIT_COMPLEX:
			channel = 2;
			depth = channel*64;			
			flage = true;
			break;
		case FIT_RGB16:
			channel = 3;
			depth = channel*16;		
			flage = true;
			break;
		case FIT_RGBA16:
			channel = 4;
			depth = channel*16;
			flage = true;
			break;
		case FIT_RGBF:
			channel = 3;
			depth = channel*32;
			flage = true;
			break;
		case FIT_RGBAF:
			channel = 4;
			depth = channel*32;
			flage = true;
			break;
		default:
			flage = false;
			break;
		}

		return flage;
	}

	void fillInfo(FIBITMAP* data)
	{
		bitsSize = width*height*channel;
		redMask = FreeImage_GetRedMask(data);
		greenMask = FreeImage_GetGreenMask(data);
		blueMask = FreeImage_GetBlueMask(data);
		alphaMask=FreeImage_GetTransparencyCount(data);

		FITAG* tagMake = nullptr;
		FreeImage_GetMetadata(FIMD_EXIF_MAIN, data, "Make", &tagMake);
		if (tagMake) hInfo.make = std::string((char*)FreeImage_GetTagValue(tagMake));

		FITAG* tagModel = nullptr;
		FreeImage_GetMetadata(FIMD_EXIF_MAIN, data, "Model", &tagModel);
		if (tagModel) hInfo.model = std::string((char*)FreeImage_GetTagValue(tagModel));

		FITAG* tagFocal = nullptr;
		FreeImage_GetMetadata(FIMD_EXIF_EXIF, data, "FocalLength", &tagFocal);
		if (tagFocal)
		{
			std::string tF = FreeImage_TagToString(FIMD_EXIF_EXIF, tagFocal);

			if (!tF.empty())
			{
				size_t index = tF.find_first_of(' ');
				if (std::string::npos != index)
				{
					hInfo.focal = static_cast<float>(std::atof(tF.substr(0,index).c_str()));
				}
			}
		}
	}

public:
	FIBITMAP* data;
	unsigned long width;
	unsigned long height;
	IMAGE_TYPE imgType;
	IMAGE_FORMAT imgFormat;
	std::string imgExt;
	unsigned char channel;//颜色通道数
	unsigned char depth; //像素点的深度（也就是位数）
	unsigned int redMask;
	unsigned int greenMask;
	unsigned int blueMask;
	unsigned int alphaMask;
	unsigned long bitsSize; //以字节为单位返回图像数据区的字节数
	HeaderInfo hInfo;
};
