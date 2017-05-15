#include "opengl.h"

Image OpenGLTools::Image2GL(const Image& img)
{
	Image imgGL;
	imgGL.width = img.width;
	imgGL.height = img.height;
	FIBITMAP* t = NULL;

	if (img.depth <= 24)
	{
		t = FreeImage_ConvertTo24Bits(img.data);

		//FreeImage_Save(img.imgFormat,t,"temp/pngtest.png");

		imgGL.data = FreeImage_Allocate(img.width, img.height, 24, FI_RGBA_BLUE_MASK, FI_RGBA_GREEN_MASK, FI_RGBA_RED_MASK);

		imgGL.depth = 24;
		imgGL.channel = 3;
	}
	else
	{
		t = FreeImage_ConvertTo32Bits(img.data);

		imgGL.data = FreeImage_Allocate(img.width, img.height, 32, FI_RGBA_BLUE_MASK, FI_RGBA_GREEN_MASK, FI_RGBA_RED_MASK);
		imgGL.depth = 32;
		imgGL.channel = 4;
	}

	if (t)
	{
		for (ulong i = 0; i < img.height; i++)
		{
			for (ulong j = 0; j < img.width; j++)
			{
				BYTE *bits = (BYTE *)FreeImage_GetScanLine(t, i);
				bits += imgGL.channel*j;
				//颜色通道 FreeImage 默认是BGRA OpenGL要求是RBGA
				Pixel<BYTE> p;

				p.red = bits[FI_RGBA_BLUE];
				p.green = bits[FI_RGBA_GREEN];
				p.blue = bits[FI_RGBA_RED];

				if (imgGL.channel == 4) p.alpha = bits[FI_RGBA_ALPHA];

				imgGL.setPixel(i, j, p);
			}
		}
		FreeImage_Unload(t);
	}
	FreeImage_FlipVertical(imgGL.data);

	imgGL.imgType = FreeImage_GetImageType(imgGL.data);
	imgGL.imgFormat = img.imgFormat;
	imgGL.imgExt = img.imgExt;

	imgGL.bitsSize = imgGL.width*imgGL.height*imgGL.channel;
	imgGL.redMask = FreeImage_GetRedMask(imgGL.data);
	imgGL.greenMask = FreeImage_GetGreenMask(imgGL.data);
	imgGL.blueMask = FreeImage_GetBlueMask(imgGL.data);
	imgGL.alphaMask = FreeImage_GetTransparencyCount(imgGL.data);

	return imgGL;

}


Image OpenGLTools::ImageToSFML(const Image& src)
{
	Image img;
	img.channel = 4;
	img.imgType = src.imgType;
	img.width = src.width;
	img.height = src.height;
	img.depth = 32;
	img.imgFormat = src.imgFormat;
	img.imgExt = src.imgExt;
	img.data = FreeImage_Allocate(img.width, img.height, 32, FI_RGBA_BLUE_MASK, FI_RGBA_GREEN_MASK, FI_RGBA_RED_MASK);
	FIBITMAP* t = FreeImage_ConvertTo32Bits(src.data);

	if (t)
	{
		BYTE* bt = FreeImage_GetBits(t);
		for (ulong i = 0; i < img.height; i++)
		{
			for (ulong j = 0; j < img.width; j++)
			{

				BYTE *bits = (BYTE *)FreeImage_GetScanLine(t, i);
				bits += img.channel*j;
				//颜色通道 FreeImage 默认是BGRA,SFML要求是RBGA
				Pixel<BYTE> p;

				p.red = bits[FI_RGBA_BLUE];
				p.green = bits[FI_RGBA_GREEN];
				p.blue = bits[FI_RGBA_RED];

				if (img.channel == 4) p.alpha = bits[FI_RGBA_ALPHA];

				img.setPixel(i, j, p);
			}
		}
		FreeImage_Unload(t);
	}

	FreeImage_FlipVertical(img.data);

	img.bitsSize = img.width*img.height*img.channel;
	img.redMask = FreeImage_GetRedMask(img.data);
	img.greenMask = FreeImage_GetGreenMask(img.data);
	img.blueMask = FreeImage_GetBlueMask(img.data);
	img.alphaMask = FreeImage_GetTransparencyCount(img.data);

	return img;
}