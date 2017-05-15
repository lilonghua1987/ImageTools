#include "ImageExpandTools.h"
//
ImageExpandTools::ImageExpandTools(void)
{
}


ImageExpandTools::~ImageExpandTools(void)
{
}


//forward diff and the first column elements are zereo (left to right)
void ImageExpandTools::DiffHorizontal(const Image& img, Matrix<double>& diff)
{
	if (img.channel != diff.channel || img.height != diff.row || img.width != diff.column)
		throw("The Horizontal diff is unequal with original img !");

	for (unsigned long i = 0; i < diff.row; i++)
	{
		for (unsigned long j = 1; j < diff.column; j++)
		{
			for (uchar c = 0; c < (uchar)diff.channel; c++)
			{
				diff.at(i,j,c) = img.atVal<BYTE>(i,j,c) - img.atVal<BYTE>(i,j-1,c);
			}
		}
	}
}


//The first row elements are zereo (top to bottom)
void ImageExpandTools::DiffVertical(const Image& img, Matrix<double>& diff)
{
	if (img.channel != diff.channel || img.height != diff.row || img.width != diff.column)
		throw("The Vertical diff is unequal with original img !");

	for (unsigned long i = 1; i < diff.row; i++)
	{
		for (unsigned long j = 0; j < diff.column; j++)
		{
			for (uchar c = 0; c < (uchar)diff.channel; c++)
			{
				diff.at(i,j,c) = img.atVal<BYTE>(i,j,c) - img.atVal<BYTE>(i-1,j,c);
			}
		}
	}
}


Matrix<Pixel<uchar>> ImageExpandTools::Image2Mat(const Image &img)
{
	assert((img.channel == 1) || (img.channel == 3) || (img.channel == 4));
	Matrix<Pixel<uchar>> mat(img.height,img.width);

	for (ulong i = 0; i < img.height; i++)
	{
		for (ulong j = 0; j < img.width; j++)
		{
			mat.set(ImageTools::getColorPixel(img,i,j),i,j);
		}
	}

	return mat;
}


Matrix<float> ImageExpandTools::Image2FloatMat(const Image &img)
{
	assert((img.channel == 1));
	Matrix<float> mat(img.height,img.width);

	for (ulong i = 0; i < img.height; i++)
	{
		for (ulong j = 0; j < img.width; j++)
		{
			mat.set((float)ImageTools::getPixel(img,i,j),i,j);
		}
	}

	return mat;
}


/*
* 函数功能：将非标准图片转换为标准图片（主要是将高动态图转换为标准图片）
* 参数： 
*      src:输入图片
*      EffectiveBits：每个通道的有效表示位数（EffectiveBits>0 && EffectiveBits<=16）
* 返回值：返回变换后的图片
*/
Image ImageExpandTools::Image2Standard(const Image& src,int EffectiveBits)
{
	assert(src.depth%16 == 0);
	assert(src.depth/16 >= 1);
	assert(EffectiveBits>0 && EffectiveBits<=16);

	Image result(src.width,src.height,src.imgType,src.depth/16);

	for (ulong i = 0; i < src.height; i++)
	{
		for (ulong j = 0; j < src.width; j++)
		{
			switch (src.depth/16)
			{
			case 1:
				{
					WORD p1 = static_cast<WORD>(src.at<WORD>(i,j)*255/pow(2,EffectiveBits));
					result.at<BYTE>(i,j) = (BYTE)p1;
					break;
				}
			case 3:
				{
					FIRGB16 p = src.at<FIRGB16>(i,j);				
					Pixel<WORD> p2(p.red,p.green,p.blue);
					result.setPixel(i,j,p2*255/pow(2,EffectiveBits));
				}
			case 4:
				{
					FIRGBA16 p = src.at<FIRGBA16>(i,j);
					Pixel<WORD> p3(p.red,p.green,p.blue,p.alpha);
					result.setPixel(i,j,p3*255/pow(2,EffectiveBits));
					break;
				}
			}			
		}
	}

	//ImageTools::SaveImage(result,"temp//saveImg.bmp");
	return result;
}


Image ImageExpandTools::InPaint(Image& src ,Image& mask ,uint r)
{
	assert(src.width == mask.width);
	assert(src.height == mask.height);
	assert(mask.channel == 1);

	assert(r >= 1 && r < src.height && r < src.width);

	Image img = ImageTools::ImageClone(src);

	for (ulong i = 0; i < src.height; i++)
	{
		for (ulong j = 0; j < src.width; j++)
		{
			if ( (int)mask.at<BYTE>(i,j))
			{
				//cout<<(int)mask.at<BYTE>(i,j)<<endl;
				int sumI = 0;
				Pixel<int> sumP;
				for (int x = -(int)r; x <= (int)r; x++)
				{
					for (int y = -(int)r; y <= (int)r; y++)
					{
						if (!(int)mask.at<BYTE>(tools::bound(i+x,0,src.height-1),tools::bound(j+y,0,src.width-1)))
						{
							sumP+=src.getPixel<BYTE>(tools::bound(i+x,0,src.height-1),tools::bound(j+y,0,src.width-1));
							sumI++;
						}
					}
				}
				img.setPixel<BYTE>(i,j,sumI?(sumP/sumI):Pixel<BYTE>());
				//img.setPixel<BYTE>(i,j,Pixel<BYTE>());
			}
		}
	}

	return img;
}


void ImageExpandTools::reSizeImage(const string dir, const string fileExtension, const string saveDir, Point2D<int> top, Point2D<int> bottom)
{
	vector<string>* lf = new vector<string>();
	tools::listFileByDir(dir,fileExtension,lf);

	tools::CreatDir(saveDir.c_str());

	for (size_t i = 0; i < lf->size(); i++)
	{
		Image img(lf->at(i).c_str());
		Image saveImg = ImageTools::ImageCopy(img,top,bottom);
		string name = saveDir+"\\";
		name.append(tools::getFileName(lf->at(i)));
		std::thread t(ImageTools::SaveImage,saveImg,name.c_str());
		cout<<"current file:"<<lf->at(i)<<endl;
		t.join();	
		//_sleep(1000);
		//SaveImage(saveImg,tools::fileNameFromTime(saveDir,string(),fileExtension).c_str());		
	}

	delete lf;
}


void ImageExpandTools::inPaintAllImages(const string dir, const string fileExtension, const string saveDir)
{
	vector<string>* lf = new vector<string>();
	tools::listFileByDir(dir,fileExtension,lf);

	tools::CreatDir(saveDir.c_str());

	for (size_t i = 0; i < lf->size(); i++)
	{
		Image img(lf->at(i).c_str());
		Image erImg = ImageTools::Threshold(img,220);
		Image saveImg = ImageExpandTools::InPaint(img,erImg,25);
		string name = saveDir+"\\";
		name.append(tools::getFileName(lf->at(i)));
		std::thread t(ImageTools::SaveImage,saveImg,name.c_str());
		cout<<"current file:"<<lf->at(i)<<endl;
		t.join();	
		//_sleep(500);	
	}

	delete lf;
}


Image ImageExpandTools::ImageToJPEG(const Image& src)
{
	Image img(src.width,src.height,FIT_BITMAP,src.channel);
	img.depth = src.depth;
	img.imgFormat = FIF_JPEG;
	img.imgExt = ".jpg";

	for (ulong i = 0; i < src.height; i++)
	{
		for (ulong j = 0; j < src.width; j++)
		{
			Pixel<BYTE> p = src.getPixel<BYTE>(i,j);
			img.setPixel(i,j,p);
		}
	}

	return img;
}


Matrix<float> ImageExpandTools::loadFromFile(std::string fileName)
{
	std::ifstream fin(fileName);
	assert(!fin.fail());
	string line;
	int wideth = 0, height = 0, index = 0;

	if(!fin.eof() && fin.good())
	{
		fin >> line;
		//std::getline(fin,line);
		height = atoi(line.c_str());
	}

	if(!fin.eof() && fin.good())
	{
		fin >> line;
		wideth = atoi(line.c_str());
	}

	Matrix<float> m(height,wideth);

	while (!fin.eof() && fin.good() && (index < wideth* height)) 
	{
		fin >> line;
		m.at(index/wideth,index%wideth) = (float)atof(line.c_str());
		index++;
		//cout<<index<<endl;
	}

	fin.close();

	return m;
}


void ImageExpandTools::imShow(const string& windowName, const Image& img)
{

}


Matrix<float> ImageExpandTools::logImage(const Image& img)
{
	Matrix<float> logMat(img.height,img.width,img.channel);

	for (ulong r = 0; r < logMat.row; r++)
	{
		for (ulong c = 0; c < logMat.column; c++)
		{
			Pixel<float> pixel = img.getPixel<float>(r,c);

			if(logMat.channel > 0) logMat.at(r,c,0) = log(pixel.red);
			if(logMat.channel > 1) logMat.at(r,c,1) = log(pixel.green);
			if(logMat.channel > 2) logMat.at(r,c,2) = log(pixel.blue);
			if(logMat.channel > 3) logMat.at(r,c,3) = log(pixel.alpha);
		}
	}

	return logMat;
}


Matrix<float> ImageExpandTools::logImage(const Matrix<float>& mat)
{
	Matrix<float> logMat(mat);

	for (ulong r = 0; r < logMat.row; r++)
	{
		for (ulong c = 0; c < logMat.column; c++)
		{
			for (ulong i = 0; i < logMat.channel; i++)
			{
				logMat.at(r,c,i) = log(mat.at(r,c,i));
			}
		}
	}

	return logMat;
}