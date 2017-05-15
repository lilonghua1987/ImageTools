#include "ImageProcess.h"

ImageProcess::ImageProcess(void)
{
}


ImageProcess::~ImageProcess(void)
{
}


Image ImageProcess::imageReverse(const Image& img){

	Image reverseImg(img);


	for(unsigned y = 0; y < img.height; y++)
	{
		for(unsigned x = 0; x < img.width; x++)
		{
			reverseImg.setPixel(y,x,255 - reverseImg.getPixel<BYTE>(y,x));
		}
	}
	if (RUNTIMECHECK)
	{
		string imgName = tools::fileName(PATH,"reverseImg_",img.width,img.height,img.imgExt);
		reverseImg.save(imgName);
	}

	return reverseImg;
}


GreyHistogram ImageProcess::getGreyHistogram(const Image& img){

	GreyHistogram histogram;
	assert((img.imgType == FIT_BITMAP) && (img.channel == 1));

	for(unsigned y = 0; y < img.height; y++) 
	{
		for(unsigned x = 0; x < img.width; x++) 
		{
			bool push = true;
			for(ulong i=0;i<histogram.colorValue.size();i++){
				if(histogram.colorValue[i] == img.at<int>(y,x)){
					histogram.colorTimes[i]++;
					push = false;
					break;
				}
			}

			if(push){
				histogram.colorValue.push_back(img.at<int>(y,x));
				histogram.colorTimes.push_back(1);
			}
		}
	}

	//排序 由小到大排序
	for (ulong i = 0; i < histogram.colorValue.size(); i++)
	{
		int a = histogram.colorValue.at(i),swapI = 0,minA = a,flage = false;
		for (ulong j = i+1; j < histogram.colorValue.size(); j++)
		{
			int b = histogram.colorValue.at(j);
			if (minA > b)
			{
				minA = b;
				swapI = j;
				flage = true;
			}
		}

		if (flage)
		{
			int temp = histogram.colorValue.at(swapI);
			histogram.colorValue.at(i) = temp;
			histogram.colorValue.at(swapI) = a;

			temp = histogram.colorTimes.at(swapI);
			a = histogram.colorTimes.at(i);
			histogram.colorTimes.at(i) = temp;
			histogram.colorTimes.at(swapI) = a;
		}
	}


	histogram.pixelNU = img.width * img.height;

	for(ulong i=0;i<histogram.colorValue.size();i++){
		histogram.colorRate.push_back(histogram.colorTimes[i]*1.0/histogram.pixelNU);
	}

	return histogram;
}


void ImageProcess::printHistogram(GreyHistogram histogram){

	FILE* fp1;
	fopen_s(&fp1, "temp\\Histogram_result.txt", "w"); 
	fprintf(fp1,"greyHistogram.pixelNU=%10d\n",histogram.pixelNU);

	for(ulong i=0;i<histogram.colorValue.size();i++){
		fprintf(fp1,"\tgreyHistogram.colorValue=%10d",histogram.colorValue[i]);
		fprintf(fp1,"\tgreyHistogram.colorTimes=%10d",histogram.colorTimes[i]);
		fprintf(fp1,"\tgreyHistogram.colorRate=%10f\n",histogram.colorRate[i]);
	}
}


Image ImageProcess::HistogramEqualization(const Image& img){

	assert(img.channel == 1);

	GreyHistogram histogram = getGreyHistogram(img);

	int cSize = histogram.colorValue.size();
	int minP = histogram.colorValue.at(0);
	int maxP = histogram.colorValue.at(cSize-1);

	if (maxP<minP)
	{
		int t = maxP;
		maxP = minP;
		minP = t;
	}

	Image result(img);;

	for (ulong i = 0; i < result.height; i++)
	{
		for (ulong j = 0; j < result.width; j++)
		{
			uchar p = img.at<uchar>(i,j);
			for (int m = 0; m < cSize; m++)
			{
				if(p == histogram.colorValue.at(m))
				{
					result.at<uchar>(i,j) = static_cast<uchar>(histogram.colorRate.at(m)*(maxP-minP)+minP);
					break;
				}
			}
		}
	}

	return result;
}

Image ImageProcess::GausscianBlur(const Image& img,double sigma){

	int range = static_cast<int>(ceil(fabs(sigma)*6 + 1)),mid = 0;


	if(range%2 == 0)
		range = range-1;

	MatrixMat matrix;

	double gauss = 1.0/(2*PI*sigma*sigma),sum=0;;

	matrix.height = matrix.width = range;

	mid = (range-1)/2 +1;

	cout<<"Gauss Matrix height:"<<matrix.height<<",sigma="<<sigma<<endl;

	//cout<<"Create gauss Matrix begin"<<endl;

	for(int y=0;y<matrix.height;y++){
		vector<double> xV;
		for(int x=0;x<matrix.width;x++){

			double shaftX = 1.0*(x-mid+1)*(x-mid+1);
			double shaftY = 1.0*(y-mid+1)*(y-mid+1);
			double shaft = -(shaftX + shaftY)/(2*sigma*sigma);
			xV.push_back(gauss*exp(shaft));
			sum+=gauss*exp(shaft);
		}
		matrix.data.push_back(xV);
	}
	cout<<"Gauss Matrix sum="<<sum<<endl;
	//cout<<"Create gauss Matrix end"<<endl;

	//归一化,确保高斯权值在[0,1]之间  
	for(int i=0;i<matrix.height;i++){
		for(int j=0;j<matrix.width;j++){
			//cout<<"Gauss Matrix "<<"["<<i<<","<<j<<"]:"<<matrix.data[i][j]<<endl;
			matrix.data[i][j]/=sum;
		}
	}

	// 得到FIBITMAP的信息 
	FREE_IMAGE_TYPE img_type =  img.imgType;  //图像类型
	ulong img_width = img.width;  //图像宽度
	ulong img_height =img.height;  //图像高度

	//存储最后结果
	Image outImg(img);

		for(ulong y = 0; y < img_height; y++) {
			for(ulong x = 0; x < img_width; x++) {
				Pixel<double> tempP;
				for(int i=0;i<range;i++){
					for(int j=0;j<range;j++){
						int changeX = j-mid+1+x;
						int changeY = i-mid+1+y;

						int flageX = 0;
						int flageY = 0;

						if(changeX>=0 && changeX<(int)img_width){
							flageX = 1;
						}
						if(changeY>=0 && changeY<(int)img_height){
							flageY = 1;
						}


						if(flageX && flageY){

							tempP += img.getPixel<BYTE>(changeY,changeX)*matrix.data[i][j];
						}
					}
				}
				if((img_type == FIT_BITMAP) && (img.channel == 1))
				{
					ImageTools::setPixel(outImg,y,x,static_cast<int>(tempP.red));
				}else
				{
					ImageTools::setColorPixel(outImg,y,x,tempP);
				}
				
			}
		}

	string imgName = tools::fileName(PATH,"gaussImg_",img_width,img_height,".jpg");
	ImageTools::SaveImage(outImg, imgName.c_str());

	return outImg;
}


Image ImageProcess::GausscianSeparateBlur(const Image& img,double sigma){
	sigma = sigma > 0 ? sigma : -sigma;
	//高斯核矩阵的大小为(6*sigma+1)*(6*sigma+1)
	//ksize为奇数
	int ksize = static_cast<int>(ceil(sigma * 3) * 2 + 1);

	//计算一维高斯核
	std::vector<double> kernel(ksize);

	double scale = -0.5/(sigma*sigma);
	double cons = 1/sqrt(-scale / PI);

	double sum = 0;
	int kcenter = ksize/2;
	int i = 0, j = 0;
	for(i = 0; i < ksize; i++)
	{
		int x = i - kcenter;
		kernel.at(i) = cons * exp(x * x * scale);//一维高斯函数
		sum += kernel.at(i);
	}

	//归一化,确保高斯权值在[0,1]之间
	for(i = 0; i < ksize; i++)
	{
		kernel.at(i) /= sum;
	}

	//存储最后结果
	Image outImg(img);

	//x方向一维高斯模糊
	for(ulong y = 0; y < img.height; y++)
	{
		for(ulong x = 0; x < img.width; x++)
		{
			Pixel<float> tempP;
			sum = 0;
			for(i = -kcenter; i <= kcenter; i++)
			{
				tempP += img.get<float>(y,tools::bound(x+i,0,img.width-1))*kernel.at(kcenter+i);
				sum += kernel.at(kcenter+i);
			}

			outImg.set(y,x,tempP/sum);
		}
	}


	//y方向一维高斯模糊
	for(ulong x = 0; x < img.width; x++)
	{
		for(ulong y = 0; y < img.height; y++)
		{
			Pixel<float> tempP;
			sum = 0;
			for(i = -kcenter; i <= kcenter; i++)
			{
				tempP += img.get<float>(tools::bound(y+i,0,img.height-1),x)*kernel.at(kcenter+i);
				sum += kernel.at(kcenter+i);
			}
			outImg.set(y,x,tempP/sum);
		}
	}

	string imgName = tools::fileName(PATH,"gaussSeparateImg_",img.width,img.height,".jpg");
	outImg.save(imgName.c_str());
	return outImg;
}


Matrix<float> ImageProcess::GausscianSeparateBlur(const Matrix<float>& img,double sigma){
	assert(img.channel == 1);
	sigma = sigma > 0 ? sigma : -sigma;
	//高斯核矩阵的大小为(6*sigma+1)*(6*sigma+1)
	//ksize为奇数
	int ksize = static_cast<int>(ceil(sigma * 3) * 2 + 1);

	//计算一维高斯核
	std::vector<double> kernel(ksize);

	double scale = -0.5/(sigma*sigma);
	double cons = 1/sqrt(-scale / PI);

	double sum = 0;
	int kcenter = ksize/2;

	for(int i = 0; i < ksize; i++)
	{
		int x = i - kcenter;
		kernel.at(i) = cons * exp(x * x * scale);//一维高斯函数
		sum += kernel.at(i);

		//cout << "old " << *(kernel+i);
	}

	//归一化,确保高斯权值在[0,1]之间
	for(int i = 0; i < ksize; i++)
	{
		kernel.at(i) /= sum;
		//cout << "new " << *(kernel+i);
	}

	ulong img_width = img.column;  //图像宽度
	ulong img_height =img.row;  //图像高度

	//存储最后结果
	Matrix<float> outImg(img_height,img_width);

	//x方向一维高斯模糊
	for(ulong y = 0; y < img_height; y++)
	{
		for(ulong x = 0; x < img_width; x++)
		{
			double tempP = 0;
			sum = 0;
			for(auto i = -kcenter; i <= kcenter; i++)
			{
				tempP += img.get(y,tools::bound(x+i,0,img_width-1))*kernel.at(kcenter+i);
				sum += kernel.at(kcenter+i);
			}

			outImg.at(y,x) = static_cast<float>(tempP/sum);			
		}
	}


	//y方向一维高斯模糊
	for(ulong x = 0; x < img_width; x++)
	{
		for(ulong y = 0; y < img_height; y++)
		{
			double tempP = 0.0;
			sum = 0;
			for(auto i = -kcenter; i <= kcenter; i++)
			{
				tempP += img.get(tools::bound(y+i,0,img_height-1),x)*kernel.at(kcenter+i);
				sum += kernel.at(kcenter+i);
			}
			outImg.at(y,x) = static_cast<float>(tempP/sum);

		}
	}

	return outImg;
}


Image ImageProcess::halfSize(const Image& img){

	Image halfImg =ImageTools::ImageCopy(img,0,0,img.width/2,img.height/2);
	
	for(ulong y = 0; y < img.height/2; y++) {
		for(ulong x = 0; x < img.width/2; x++) {
			halfImg.setPixel<BYTE>(y,x,img.getPixel<BYTE>(2*y,2*x));
		}
	}

	if (RUNTIMECHECK)
	{
		string imgName = tools::fileName(PATH,"halfImg_",img.width,img.height,".jpg");
		halfImg.save(imgName);
	}

	return halfImg;
}


Image ImageProcess::doubleSize(const Image& img){

	Image doubleImg(img.width*2, img.height*2,img.imgType,img.channel);

	for(ulong y = 0; y < img.height*2; y++)
	{
		for(ulong x = 0; x < img.width*2; x++) 
		{
			doubleImg.setPixel(y,x,img.getPixel<BYTE>(y/2,x/2));
		}
	}

	if (RUNTIMECHECK)
	{
		string imgName = tools::fileName(PATH,"doubleImg_",img.width,img.height,doubleImg.imgExt);
		doubleImg.save(imgName);
	}

	return doubleImg;
}


MatrixMat ImageProcess::inverseMatrixThree(const MatrixMat& matrix){

	MatrixMat m;

	if(matrix.height == matrix.width && matrix.height == 3){
		double value;
		value = matrix.data[0][0]*matrix.data[1][1]*matrix.data[2][2];
		value += matrix.data[0][1]*matrix.data[1][2]*matrix.data[2][0];
		value += matrix.data[0][2]*matrix.data[1][0]*matrix.data[2][1];
		value -= matrix.data[0][2]*matrix.data[1][1]*matrix.data[2][0];
		value -= matrix.data[1][2]*matrix.data[2][1]*matrix.data[0][0];
		value -= matrix.data[2][2]*matrix.data[1][0]*matrix.data[0][1];

		if(value != 0){
			vector<double> row1;
			double element;
			element = matrix.data[1][1]*matrix.data[2][2] - matrix.data[1][2]*matrix.data[2][1];
			row1.push_back(element/value);
			element = matrix.data[0][1]*matrix.data[2][2] - matrix.data[2][1]*matrix.data[0][2];
			row1.push_back(-element/value);
			element = matrix.data[0][1]*matrix.data[1][2] - matrix.data[0][2]*matrix.data[2][1];
			row1.push_back(element);
			m.data.push_back(row1);

			vector<double> row2;
			element = matrix.data[1][2]*matrix.data[2][0] - matrix.data[2][2]*matrix.data[1][0];
			row2.push_back(element/value);
			element = matrix.data[0][2]*matrix.data[2][0] - matrix.data[2][2]*matrix.data[0][0];
			row2.push_back(-element/value);
			element = matrix.data[0][2]*matrix.data[1][0] - matrix.data[0][0]*matrix.data[1][2];
			row2.push_back(element);
			m.data.push_back(row2);

			vector<double> row3;
			element = matrix.data[1][0]*matrix.data[2][1] - matrix.data[2][0]*matrix.data[1][1];
			row3.push_back(element/value);
			element = matrix.data[0][0]*matrix.data[2][1] - matrix.data[2][0]*matrix.data[0][1];
			row3.push_back(-element/value);
			element = matrix.data[0][0]*matrix.data[1][1] - matrix.data[0][1]*matrix.data[1][0];
			row3.push_back(element);
			m.data.push_back(row3);
		}
	}

	return m;

}


Image ImageProcess::NCC(const Image& imgL,const Image& imgR,int W,int dispMax){

	assert(imgL.width == imgR.width && imgL.height == imgR.height);
	assert(imgL.channel == imgR.channel && imgL.channel == 1);

	Matrix<unsigned long> intImgL = ImageExpandTools::IntegralImage<unsigned long>(imgL);
	Matrix<unsigned long> intImgR = ImageExpandTools::IntegralImage<unsigned long>(imgR);

	Matrix<float> pm(imgL.height,imgL.width,dispMax+1);

	int midW = (W-1)/2;

	Matrix<float> meanImgL(intImgL.row,intImgL.column,intImgL.channel);
	Matrix<float> meanImgR(intImgR.row,intImgR.column,intImgR.channel);

	for (ulong r = 0 ; r < imgL.height; r++)
	{
		for (ulong c = 0; c < imgL.width; c++)
		{
			if (r <= (ulong)midW)
			{
				if (c <= (ulong)midW)
				{
					meanImgL.at(r,c) = (float)intImgL.at(r+midW,c+midW)/((r+midW+1)*(c+midW+1));
					meanImgR.at(r,c) = (float)intImgR.at(r+midW,c+midW)/((r+midW+1)*(c+midW+1));
				}else if ((ulong)midW < c && c < (imgL.width-1-midW))
				{
					meanImgL.at(r,c) = (float)(intImgL.at(r+midW,c+midW)-intImgL.at(r+midW,c-midW-1))/((r+midW+1)*(2*midW+1));
					meanImgR.at(r,c) = (float)(intImgR.at(r+midW,c+midW)-intImgR.at(r+midW,c-midW-1))/((r+midW+1)*(2*midW+1));
				}else
				{
					meanImgL.at(r,c) = (float)(intImgL.at(r+midW,imgL.width-1)-intImgL.at(r+midW,c-midW-1))/((r+midW+1)*(imgL.width-c+midW));
					meanImgR.at(r,c) = (float)(intImgR.at(r+midW,imgL.width-1)-intImgR.at(r+midW,c-midW-1))/((r+midW+1)*(imgL.width-c+midW));
				}
			}else if ( r > (ulong)midW && r < (imgL.height-1-midW))
			{
				if (c <= (ulong)midW)
				{
					meanImgL.at(r,c) = (float)(intImgL.at(r+midW,c+midW)-intImgL.at(r-midW-1,c+midW))/((2*midW+1)*(c+midW+1));
					meanImgR.at(r,c) = (float)(intImgR.at(r+midW,c+midW)-intImgR.at(r-midW-1,c+midW))/((2*midW+1)*(c+midW+1));
				}else if ((ulong)midW < c && c < (imgL.width-1-midW))
				{
					meanImgL.at(r,c) = (float)(intImgL.at(r+midW,c+midW)-intImgL.at(r+midW,c-midW-1)-intImgL.at(r-midW-1,c+midW)+intImgL.at(r-midW-1,c-midW-1))/((2*midW+1)*(2*midW+1));
					meanImgR.at(r,c) = (float)(intImgR.at(r+midW,c+midW)-intImgR.at(r+midW,c-midW-1)-intImgR.at(r-midW-1,c+midW)+intImgR.at(r-midW-1,c-midW-1))/((2*midW+1)*(2*midW+1));
				}else
				{
					meanImgL.at(r,c) = (float)(intImgL.at(r+midW,imgL.width-1)-intImgL.at(r+midW,c-midW-1)-intImgL.at(r-midW-1,imgL.width-1)+intImgL.at(r-midW-1,c-midW-1))/((2*midW+1)*(imgL.width-c+midW));
					meanImgR.at(r,c) = (float)(intImgR.at(r+midW,imgL.width-1)-intImgR.at(r+midW,c-midW-1)-intImgR.at(r-midW-1,imgL.width-1)+intImgR.at(r-midW-1,c-midW-1))/((2*midW+1)*(imgL.width-c+midW));
				}
			}else
			{
				if (c <= (ulong)midW)
				{
					meanImgL.at(r,c) = (float)(intImgL.at(imgL.height-1,c+midW)-intImgL.at(r-midW-1,c+midW))/((imgL.height-r+midW)*(c+midW+1));
					meanImgR.at(r,c) = (float)(intImgR.at(imgL.height-1,c+midW)-intImgR.at(r-midW-1,c+midW))/((imgL.height-r+midW)*(c+midW+1));
				}else if ((ulong)midW < c && c < (imgL.width-1-midW))
				{
					meanImgL.at(r,c) = (float)(intImgL.at(imgL.height-1,c+midW)-intImgL.at(imgL.height-1,c-midW-1)-intImgL.at(r-midW-1,c+midW)+intImgL.at(r-midW-1,c-midW-1))/((imgL.height-r+midW)*(2*midW+1));
					meanImgR.at(r,c) = (float)(intImgR.at(imgL.height-1,c+midW)-intImgR.at(imgL.height-1,c-midW-1)-intImgR.at(r-midW-1,c+midW)+intImgL.at(r-midW-1,c-midW-1))/((imgL.height-r+midW)*(2*midW+1));
				}else
				{
					meanImgL.at(r,c) = (float)(intImgL.at(imgL.height-1,imgL.width-1)-intImgL.at(imgL.height-1,c-midW-1)-intImgL.at(r-midW-1,imgL.width-1)+intImgL.at(r-midW-1,c-midW-1))/((imgL.height-r+midW)*(imgL.width-c+midW));
					meanImgR.at(r,c) = (float)(intImgR.at(imgL.height-1,imgL.width-1)-intImgR.at(imgL.height-1,c-midW-1)-intImgR.at(r-midW-1,imgL.width-1)+intImgR.at(r-midW-1,c-midW-1))/((imgL.height-r+midW)*(imgL.width-c+midW));
				}
			}
		}
	}

	for(ulong x=0;x<imgL.height;x++){
		for(ulong y=0;y<imgL.width;y++){
			float mR = meanImgR.at(x,y);
			for(int dispRange=0;dispRange<=dispMax;dispRange++){
				float r1=0.0f,r2=0.0f,r3=0.0f,r=0.0f;
				for(int h=-midW;h<=midW;h++){
					for(int w=-midW;w<=midW;w++){
						if ((dispRange+y+w)+1 > imgL.width || (dispRange+y+w) < 0)
						{
							break;
						}
						float mL = meanImgL.at(x,dispRange+y+w);
						/*if ((x+h) >= 0 && (x+h) < referenceImg_height && (y+w) >= 0 && (y+w) < referenceImg_width)
						{
						}*/
						int a,b;
						float t1 = 0.0f,t2 = 0.0f;
						a=imgL.at<BYTE>(tools::bound(x+h,0,imgL.height-1),tools::bound(dispRange+y+w,0,imgL.width-1));
						b=imgR.at<BYTE>(tools::bound(x+h,0,imgL.height-1),tools::bound(y+w,0,imgL.width-1));
						t1 = (a - mL);
						t2 = (b - mR);
						r1 += t1*t2;
						r2 += t1*t1;
						r3 += t2*t2;
					}
				}

				r2 =sqrt(r2);
				r3 = sqrt(r3);
				r =(r2*r3)? r1/(r2*r3):(r2*r3);

				pm.at(x,y,dispRange) = r;
			}

		}
	}

 	Image img(imgL.width,imgL.height,imgL.imgType,1);

	for (ulong i = 0; i < imgL.height; i++)
	{
		for (ulong j = 0; j < imgL.width; j++)
		{
			float maxF = pm.at(i,j,0);
			int bestD = 0;
			for (int d = 0; d <= dispMax; d++)
			{
				if (maxF < pm.at(i,j,d))
				{
					maxF = pm.at(i,j,d);
					bestD = d;
				}
			}

			img.at<BYTE>(i,j) = bestD;
		}
	}

	ImageTools::NormalImage(img,0,255);

	return img;
}


void ImageProcess::printMapingPointMatrix(ParallaxMatrix parallaxMatrix,FREE_IMAGE_TYPE img_type,FREE_IMAGE_FORMAT fif){

	FILE* fp6;
	fopen_s(&fp6, "temp\\parallaxMatrix.txt", "w"); 
	fprintf(fp6,"\trelation：");
	fprintf(fp6,"\tparallax：\n");

	int maxR = 0;

	for(int y=0;y<parallaxMatrix.height;y++){
		for(int x=0;x<parallaxMatrix.width;x++){
			ParallaxPoint f = parallaxMatrix.value[y][x];
			fprintf(fp6,"\t%5f",f.relation);
			fprintf(fp6,"\t%5d\n",f.parallax);

			if(maxR<f.parallax)
				maxR = f.parallax;
		}
	}
	fclose(fp6);

	Image parallaxImg;parallaxImg.createImage(parallaxMatrix.width, parallaxMatrix.height, img_type);

	for(int y=0;y<parallaxMatrix.height;y++){
		for(int x=0;x<parallaxMatrix.width;x++){
			ParallaxPoint f = parallaxMatrix.value[y][x];
			ImageTools::setPixel(parallaxImg,y,x,f.parallax * 255/maxR);
		}
	}

	string imgExt;
	switch (fif)
	{
	case FIF_BMP:
		imgExt.append(".bmp");
		break;
	case FIF_JPEG:
		imgExt.append(".jpg");
		break;
	case FIF_PNG:
		imgExt.append(".png");
		break;
	default:
		imgExt.append(".jpg");
		break;
	}
	string imgName = tools::fileNameFromTime(PATH,"parallaxMatrix",imgExt);
	ImageTools::SaveImage(parallaxImg, imgName.c_str());
}


//数组a为输入，数组y为输出，2的power次方为数组的长度 
void ImageProcess::fft(const complex<double> a[], complex<double> y[], int power) 
{ 
	if(0 == power) 
	{ 
		y[0] = a[0]; 
		return; 
	} 
	int n = 1 << power; 
	double angle = 2 * PI / n; 
	complex<double> wn(cos(angle), sin(angle)); 
	complex<double> w(1, 0); 
	complex<double> *a0 = new complex<double>[n / 2]; 
	complex<double> *a1 = new complex<double>[n / 2]; 
	complex<double> *y0 = new complex<double>[n / 2]; 
	complex<double> *y1 = new complex<double>[n / 2]; 
	for(int i = 0; i < n / 2; i ++) 
	{ 
		a0[i] = a[2 * i]; 
		a1[i] = a[2 * i + 1]; 
	} 
	//分开成两个子fft过程 
	fft(a0, y0, power - 1); 
	fft(a1, y1, power - 1); 
	complex<double> u; 
	for(int k = 0; k < n / 2; k++) //蝶形算法 
	{ 
		u = w * y1[k]; 
		y[k] = y0[k] + u; 
		y[k + n / 2] = y0[k] - u; 
		w = w * wn; 
	} 
	delete[] a0; 
	delete[] a1; 
	delete[] y0; 
	delete[] y1; 
} 


//y为输入，a为输出，2的power次方为数组的长度 
void ImageProcess::ifft(const complex<double> y[], complex<double> a[], int power) 
{ 
	int count = 1 << power; 
	complex<double> *x = new complex<double>[count]; 
	memcpy(x, y, sizeof(complex<double>) * count); 
	int i; 
	for(i = 0; i < count; i++) 
	{ 
		x[i] = complex<double>(x[i].real(), -x[i].imag()); //共轭复数 
	} 
	fft(x, a, power); //调用快速傅立叶变换算法 
	for(i = 0; i < count; i++) 
	{ 
		a[i] = complex<double>(a[i].real() / count, -a[i].imag() / count); //共轭复数 
	} 
	delete[] x; 
} 


//宽和高都截取为2的指数倍 
Image ImageProcess::adjustImageSize(Image image,ulong &wp,ulong &hp) 
{ 
	ulong w = 1; 
	ulong h = 1; 

	FREE_IMAGE_TYPE img_type =  image.imgType;  //图像类型
	ulong width = image.width;  //图像宽度
	ulong height =image.height;  //图像高度

	wp = 0, hp = 0; 
	while(w  < width) { w *= 2; (wp)++;} 
	while(h  < height) {h *= 2; (hp)++;} 

	if(w == width && h == height){
		return ImageTools::ImageClone(image);
	}else{

		Image adjustImg; adjustImg.createImage(w, h, img_type);

		if((img_type == FIT_BITMAP) && (image.channel == 1)) {
			for(ulong y = 0; y < h; y++) {
				for(ulong x = 0; x < w; x++) {
					if(y<height && x<width){
						ImageTools::setPixel(adjustImg,y,x,ImageTools::getPixel(image,x,y));
					}else{
						ImageTools::setPixel(adjustImg,y,x,0);
					}
				}
			}
		}

		string imgName = tools::fileName(PATH,"adjustImg_",width,height,".jpg");
		ImageTools::SaveImage(adjustImg, imgName.c_str());

		return adjustImg;
	}
}

complex<double> * ImageProcess::fourier(Image image)
{ 
	ulong wp = 0,hp= 0;
	Image adjustImg = adjustImageSize(image,wp,hp); //调整大小 
	ulong width = adjustImg.width;  //图像宽度
	ulong height =adjustImg.height;  //图像高度

	complex<double> *TD = new complex<double>[width * height]; //当前读取的数据为时域

	for(ulong i = 0; i < height; i++) 
	{ 
		for(ulong j = 0; j < width; j++) 
		{ 
			TD[i * width + j] = complex<double>( ImageTools:: getPixel(adjustImg,i,j), 0); 
		} 
	}

	complex<double> *FD = new complex<double>[width * height]; //申请空间保存变换后结果 

	for(ulong i = 0; i < height; i++) //在x方向上对按行进行快速傅立叶变换 
	{ 
		fft(&TD[width * i], &FD[width * i], wp); 
	} 

	//extern void *memcpy(void *dest, void *src, unsigned int count);
	//由src所指内存区域复制count个字节到dest所指内存区域
	memcpy(TD, FD, sizeof(complex<double>) * width * height);

	complex<double> *columnt = new complex<double>[height];// 存放将要进行变换的列数据
	complex<double> *columnf = new complex<double>[height]; //用于存放对应列变换后的数据

	for(ulong i = 0; i < width; i++) //调整行列数据，在y方向上按列进行快速傅立叶变换 
	{ 
		for(ulong j = 0; j < height; j++) 
		{ 
			columnt[j] = TD[j * width + i]; 
		} 
		fft(columnt, columnf, hp); 
		for(ulong j = 0; j < height; j++) 
		{ 
			FD[j * width + i] = columnf[j]; 
		} 
	} 
	delete[] columnt; 
	delete[] columnf; 

	double coef=0.02;

	for(ulong i = 0; i < height; i++) 
	{   
		for(ulong j = 0; j < width; j++) 
		{ 
			double spectral = abs(FD[i * width + j]) * coef; //灰度值调整 
			spectral = spectral > 255 ? 255 : spectral; 
			ImageTools:: setPixel(adjustImg,i,j,static_cast<int>(spectral));
		} 
	}

	string imgName = tools::fileNameFromTime(PATH,"fourierImg",".jpg");
	ImageTools::SaveImage(adjustImg, imgName.c_str());

	delete[] TD; 
	return FD;
}


void ImageProcess::inverseFourier(complex<double> *fd,int wp,int hp,FREE_IMAGE_TYPE img_type){

	ulong width = 1 << wp;  //图像宽度
	ulong height = 1 << hp;  //图像高度

	Image inverseFourierImg; inverseFourierImg.createImage(width, height, img_type);

	complex<double> *TD = fd;

	complex<double> *FD = new complex<double>[width * height]; //申请空间保存变换后结果 

	for(ulong i = 0; i < height; i++) //在x方向上对按行进行逆快速傅立叶变换 
	{ 
		ifft(&TD[width * i], &FD[width * i], wp); 
	} 

	//extern void *memcpy(void *dest, void *src, unsigned int count);
	//由src所指内存区域复制count个字节到dest所指内存区域
	memcpy(TD, FD, sizeof(complex<double>) * width * height);

	complex<double> *columnt = new complex<double>[height];// 存放将要进行变换的列数据
	complex<double> *columnf = new complex<double>[height]; //用于存放对应列变换后的数据

	for(ulong i = 0; i < width; i++) //调整行列数据，在y方向上按列进行逆快速傅立叶变换 
	{ 
		for(ulong j = 0; j < height; j++) 
		{ 
			columnt[j] = TD[j * width + i]; 
		} 
		ifft(columnt, columnf, hp); 
		for(ulong j = 0; j < height; j++) 
		{ 
			FD[j * width + i] = columnf[j]; 
		} 
	} 
	delete[] columnt; 
	delete[] columnf; 

	for(ulong i = 0; i < height; i++) 
	{   
		for(ulong j = 0; j < width; j++) 
		{ 
			double spectral = abs(FD[i * width + j]);
			spectral = spectral > 255 ? 255 : spectral; 
			ImageTools:: setPixel(inverseFourierImg,i,j,static_cast<int>(spectral));
		} 
	}

	string imgName = tools::fileNameFromTime(PATH,"inverseFourierImg",".jpg");
	ImageTools::SaveImage(inverseFourierImg, imgName.c_str());

	delete[] FD; 

}


void ImageProcess::printFourierData(complex<double> *fd,int wp,int hp){
	FILE* fp7;
	fopen_s(&fp7, "temp\\fFourierData.txt", "w"); 

	ulong width = 1 << wp; 
	ulong height = 1 << hp; 

	for(ulong i=0;i<height*width;i++){
		fprintf(fp7,"[%10f,%10f]",fd[i].real(),fd[i].imag());
		if((i+1)%wp == 0){
			fprintf(fp7,"\n");
		}
	}
	fclose(fp7);
}


void ImageProcess::dft(const complex<double> pr[],complex<double> fr[], int n) {
	int N = n;
	if(0 == N) 
	{ 
		fr[0] = pr[0]; 
		return; 
	} 
	double angle = 2 * PI / N; 

	for (int i = 0; i < N; i++) {
		double R=0,I=0;
		for (int j = 0; j < N; j++) {
			double pR=pr[j].real(),pI=pr[j].imag();
			R += pR * cos(angle * i * j) + pI * sin(angle * i * j);
			I += -pR * sin(angle * i * j) + pI * cos(angle * i * j);
		}
		fr[i] = complex<double>(R, I);
	}
}


void ImageProcess::idft(const complex<double> fr[],complex<double> pr[], int n) {
	int N = n;
	double angle = 2 * PI / N;
	for (int i = 0; i < N; i++) {
		double R=0,I=0;
		for (int j = 0; j < N; j++) {
			//x[n] = X[k] (cos(2πkn/N) + j sin(2πkn/N))
			/* Re X[k] （cos(2πkn/N) + j sin(2πkn/N)） + ---------------------(1)2Re X[ k] cos(2πkn/N)
			Im X[k] （ - sin(2πkn/N) + j cos(2πkn/N)) ---------------(2)-2 Im X[k] sin(2πkn/N)*/
			R += fr[j].real() * cos(angle * i * j) - fr[j].imag() * sin(angle * i * j);
			I += fr[j].real() * sin(angle * i * j) + fr[j].imag() * cos(angle * i * j);
		}
		pr[i] = complex<double>(R/ N, I/ N);
	}

}


complex<double> * ImageProcess::dFourier(Image image)
{ 
	ulong width = image.width;  //图像宽度
	ulong height = image.height;  //图像高度
	Image dFourierImg = ImageTools::ImageClone(image);
	//FIBITMAP *dPFourierImg = FreeImage_Copy(image,0,0,width,height);

	complex<double> *TD = new complex<double>[width * height]; //当前读取的数据为时域

	for(ulong i = 0; i < height; i++) 
	{ 
		for(ulong j = 0; j < width; j++) 
		{ 
			TD[i * width + j] = complex<double>( ImageTools::getPixel(dFourierImg,i,j), 0); 
		} 
	}

	complex<double> *FD = new complex<double>[width * height]; //申请空间保存变换后结果 

	for(ulong i = 0; i < height; i++) //在x方向上对按行进行离散傅立叶变换 
	{ 
		dft(&TD[width * i], &FD[width * i], width); 
	} 

	//extern void *memcpy(void *dest, void *src, unsigned int count);
	//由src所指内存区域复制count个字节到dest所指内存区域
	memcpy(TD, FD, sizeof(complex<double>) * width * height);

	complex<double> *columnt = new complex<double>[height];// 存放将要进行变换的列数据
	complex<double> *columnf = new complex<double>[height]; //用于存放对应列变换后的数据

	for(ulong i = 0; i < width; i++) //调整行列数据，在y方向上对按行进行离散傅立叶变换 
	{ 
		for(ulong j = 0; j < height; j++) 
		{ 
			columnt[j] = TD[j * width + i]; 
		} 
		dft(columnt, columnf, height); 
		for(ulong j = 0; j < height; j++) 
		{ 
			FD[j * width + i] = columnf[j]; 
		} 
	} 
	delete[] columnt; 
	delete[] columnf; 

	/*double pmax = 0,dmax = 0;
	for(int i = 0; i < height; i++) 
	{   
	for(int j = 0; j < width; j++) 
	{ 
	double ptemp,dtemp;
	dtemp =  sqrt(pow(FD[i * width + j].imag(),2.0)+pow(FD[i * width + j].real(),2.0))/100;
	dtemp = log(1+dtemp);
	ptemp = atan(FD[i * width + j].imag()/FD[i * width + j].real());

	if(dmax<dtemp){
	dmax = dtemp;
	}

	if(pmax<ptemp){
	pmax = ptemp;
	}
	} 
	}
	*/
	for(ulong i = 0; i < height; i++) 
	{   
		for(ulong j = 0; j < width; j++) 
		{ 
			double spectral = abs(FD[i * width + j])* 0.02; //灰度值调整 
			spectral = spectral > 255 ? 255 : spectral;

			/*ptemp = atan(FD[i * width + j].imag()/FD[i * width + j].real());*/

			ImageTools::setPixel(dFourierImg,i,j,static_cast<int>(spectral));
			/*setPixel(dPFourierImg,i,j,ptemp);*/
		} 
	}

    string imgName = tools::fileNameFromTime(PATH,"dFourierImg",".jpg");
	ImageTools::SaveImage(dFourierImg, imgName.c_str());

	delete[] TD; 
	return FD;
}


void ImageProcess::inverseDFourier(complex<double> *fd,ulong width,ulong height,FREE_IMAGE_TYPE img_type){

	Image inverseDFourierImg; inverseDFourierImg.createImage(width, height, img_type);

	complex<double> *TD = fd;

	complex<double> *FD = new complex<double>[width * height]; //申请空间保存变换后结果 

	for(ulong i = 0; i < height; i++) //在x方向上对按行进行逆快速傅立叶变换 
	{ 
		idft(&TD[width * i], &FD[width * i], width); 
	} 

	//extern void *memcpy(void *dest, void *src, unsigned int count);
	//由src所指内存区域复制count个字节到dest所指内存区域
	memcpy(TD, FD, sizeof(complex<double>) * width * height);

	complex<double> *columnt = new complex<double>[height];// 存放将要进行变换的列数据
	complex<double> *columnf = new complex<double>[height]; //用于存放对应列变换后的数据

	for(ulong i = 0; i < width; i++) //调整行列数据，在y方向上按列进行逆快速傅立叶变换 
	{ 
		for(ulong j = 0; j < height; j++) 
		{ 
			columnt[j] = TD[j * width + i]; 
		} 
		idft(columnt, columnf, height); 
		for(ulong j = 0; j < height; j++) 
		{ 
			FD[j * width + i] = columnf[j]; 
		} 
	} 
	delete[] columnt; 
	delete[] columnf; 

	for(ulong i = 0; i < height; i++) 
	{   
		for(ulong j = 0; j < width; j++) 
		{ 
			double spectral = abs(FD[i * width + j]);
			spectral = spectral > 255 ? 255 : spectral; 
			ImageTools::setPixel(inverseDFourierImg,i,j,static_cast<int>(spectral));
		} 
	}

	 string imgName = tools::fileNameFromTime(PATH,"inverseDFourierImg",".jpg");
	ImageTools::SaveImage(inverseDFourierImg, imgName.c_str());

	delete[] FD; 

}


void ImageProcess::printDFourierData(complex<double> *fd,ulong width,ulong height){
	FILE* fp8;
	fopen_s(&fp8, "temp\\dFourierData.txt", "w"); 

	for(ulong i=0;i<height*width;i++){
		fprintf(fp8,"[%10f,%10f]",fd[i].real(),fd[i].imag());
		if((i+1)%width == 0){
			fprintf(fp8,"\n");
		}
	}
	fclose(fp8);
}


void ImageProcess::parseDFourier(Image referenceImg,Image targetImg){

	complex<double> *FDL = dFourier(referenceImg);
	complex<double> *FDR = 	dFourier(targetImg);

	ulong width = referenceImg.width;  //图像宽度
	ulong height = referenceImg.height;  //图像高度

	complex<double> *TD = new complex<double>[width * height];

	Image parseDFourierImg = ImageTools::ImageCopy(referenceImg,0,0 ,width, height);

	for(ulong i = 0; i < height; i++) 
	{ 
		for(ulong j = 0; j < width; j++) 
		{ 
			/*res[i][0] = ( img2[i][0] * img1[i][0] ) - ( img2[i][1] * ( -img1[i][1] ) );
			res[i][1] = ( img2[i][0] * ( -img1[i][1] ) ) + ( img2[i][1] * img1[i][0] );
			tmp = sqrt( pow( res[i][0], 2.0 ) + pow( res[i][1], 2.0 ) );*/
			complex<double> tempR=FDR[i * width + j],tempL = FDL[i * width + j],temp;
			double model;
			temp = tempL * conj(tempR);
			model = sqrt(pow(temp.real(),2)+pow(temp.imag(),2));
			TD[i * width + j] = temp/model; 
		} 
	}

	FILE* fp9;
	fopen_s(&fp9, "temp\\parseDFourier.txt", "w"); 
	for(ulong i=0;i<height*width;i++){
		fprintf(fp9,"[%10f,%10f]",TD[i].real(),TD[i].imag());
		if((i+1)%width == 0){
			fprintf(fp9,"\n");
		}
	}
	fclose(fp9);

	complex<double> *FD = new complex<double>[width * height]; //申请空间保存变换后结果 

	for(ulong i = 0; i < height; i++) //在x方向上对按行进行逆快速傅立叶变换 
	{ 
		idft(&TD[width * i], &FD[width * i], width); 
	} 

	//extern void *memcpy(void *dest, void *src, unsigned int count);
	//由src所指内存区域复制count个字节到dest所指内存区域
	memcpy(TD, FD, sizeof(complex<double>) * width * height);

	complex<double> *columnt = new complex<double>[height];// 存放将要进行变换的列数据
	complex<double> *columnf = new complex<double>[height]; //用于存放对应列变换后的数据

	for(ulong i = 0; i < width; i++) //调整行列数据，在y方向上按列进行逆快速傅立叶变换 
	{ 
		for(ulong j = 0; j < height; j++) 
		{ 
			columnt[j] = TD[j * width + i]; 
		} 
		idft(columnt, columnf, height); 
		for(ulong j = 0; j < height; j++) 
		{ 
			FD[j * width + i] = columnf[j]; 
		} 
	} 
	delete[] columnt; 
	delete[] columnf;

	FILE* fp10;
	fopen_s(&fp10, "temp\\parseIDFourier.txt", "w"); 
	double maxR = 0,minR =255;

	for(ulong y=0;y<width * height;y++){
		if(maxR<FD[y].real())
			maxR = FD[y].real();
		if(minR>FD[y].real())
			minR = FD[y].real();
		fprintf(fp10,"[%10f,%10f]",FD[y].real(),FD[y].imag());
		//fprintf(fp10,"%10f",FD[y].real());
		if((y+1)%width == 0){
			fprintf(fp10,"\n");
		}
	}

	fprintf(fp10,"maxR=%10f,minR=%10f",maxR,minR);

	fclose(fp10);

	FILE* fp11;
	fopen_s(&fp11, "temp\\parsePixel.txt", "w"); 

	for(ulong i = 0; i < height; i++) 
	{   
		for(ulong j = int(maxR*width); j < width; j++) 
		{ 
			/*double spectral = FD[i * width + j].real(); 
			spectral = (spectral-minR)/(maxR-minR)*255;*/
			int value = ImageTools::getPixel(parseDFourierImg,i,j);
			ImageTools::setPixel(parseDFourierImg,i,j-int(maxR*width),value);

			fprintf(fp11,"%10d",value);
		} 
		fprintf(fp11,"\n");
	}
	fclose(fp11);

	string imgName = tools::fileNameFromTime(PATH,"inverseDFourierImg",referenceImg.imgExt);
	ImageTools::SaveImage(parseDFourierImg, imgName.c_str());

	delete[] FD; 

}


void ImageProcess::parseFourier(Image referenceImg,Image targetImg){

	complex<double> *FDL = fourier(referenceImg);
	complex<double> *FDR = 	fourier(targetImg);

	ulong wp = 0,hp= 0;
	Image adjustImg = adjustImageSize(referenceImg,wp,hp); //调整大小 

	ulong width = 1 << wp;  //图像宽度
	ulong height = 1 << hp;  //图像高度

	complex<double> *TD = new complex<double>[width * height];

	Image parseFourierImg;parseFourierImg.createImage(width, height, referenceImg.imgType);

	for(ulong i = 0; i < height; i++) 
	{ 
		for(ulong j = 0; j < width; j++) 
		{ 
			/*res[i][0] = ( img2[i][0] * img1[i][0] ) - ( img2[i][1] * ( -img1[i][1] ) );
			res[i][1] = ( img2[i][0] * ( -img1[i][1] ) ) + ( img2[i][1] * img1[i][0] );
			tmp = sqrt( pow( res[i][0], 2.0 ) + pow( res[i][1], 2.0 ) );*/
			complex<double> tempR=FDR[i * width + j],tempL = FDL[i * width + j],temp;
			double model;
			temp = tempR * conj(tempL);
			model = sqrt(pow(temp.real(),2)+pow(temp.imag(),2));
			TD[i * width + j] = temp/model; 
		} 
	}

	FILE* fp9;
	fopen_s(&fp9, "temp\\parseFourier.txt", "w"); 
	for(ulong i=0;i<height*width;i++){
		fprintf(fp9,"[%10f,%10f]",TD[i].real(),TD[i].imag());
		if((i+1)%width == 0){
			fprintf(fp9,"\n");
		}
	}
	fclose(fp9);

	complex<double> *FD = new complex<double>[width * height]; //申请空间保存变换后结果 

	for(ulong i = 0; i < height; i++) //在x方向上对按行进行逆快速傅立叶变换 
	{ 
		ifft(&TD[width * i], &FD[width * i], wp); 
	} 

	//extern void *memcpy(void *dest, void *src, unsigned int count);
	//由src所指内存区域复制count个字节到dest所指内存区域
	memcpy(TD, FD, sizeof(complex<double>) * width * height);

	complex<double> *columnt = new complex<double>[height];// 存放将要进行变换的列数据
	complex<double> *columnf = new complex<double>[height]; //用于存放对应列变换后的数据

	for(ulong i = 0; i < width; i++) //调整行列数据，在y方向上按列进行逆快速傅立叶变换 
	{ 
		for(ulong j = 0; j < height; j++) 
		{ 
			columnt[j] = TD[j * width + i]; 
		} 
		ifft(columnt, columnf, hp); 
		for(ulong j = 0; j < height; j++) 
		{ 
			FD[j * width + i] = columnf[j]; 
		} 
	} 
	delete[] columnt; 
	delete[] columnf;

	FILE* fp10;
	fopen_s(&fp10, "temp\\parseIFourier.txt", "w"); 
	double maxR = 0,minR =255;

	for(ulong y=0;y<width * height;y++){	
		if(maxR<FD[y].real())
			maxR = FD[y].real();
		if(minR>FD[y].real())
			minR = FD[y].real();
		fprintf(fp10,"[%10f,%10f]",FD[y].real(),FD[y].imag());
		if((y+1)%width == 0){
			fprintf(fp10,"\n");
		}
	}

	fprintf(fp10,"maxR=%10f,minR=%10f",maxR,minR);

	fclose(fp10);

	FILE* fp11;
	fopen_s(&fp11, "temp\\parseFourierPixel.txt", "w"); 

	for(ulong i = 0; i < height; i++) 
	{   
		for(ulong j = 0; j < width; j++) 
		{ 
			double spectral = FD[i * width + j].real(); 
			spectral = (spectral-minR)/(maxR-minR)*255;
			ImageTools::setPixel(parseFourierImg,i,j,static_cast<int>(spectral));

			fprintf(fp11,"%10f",spectral);
		} 
		fprintf(fp11,"\n");
	}
	fclose(fp11);

	string imgName = tools::fileNameFromTime(PATH,"parseFourierImg",referenceImg.imgExt);
	ImageTools::SaveImage(parseFourierImg, imgName.c_str());

	delete[] FD; 

}

//f(i+u,j+v) = (1-u)(1-v)f(i,j) + (1-u)vf(i,j+1) + u(1-v)f(i+1,j) + uvf(i+1,j+1)
Image ImageProcess::BilinearInterpolation(const Image& sourceImg,ulong height,ulong width){

	FREE_IMAGE_TYPE sourceImg_type =  sourceImg.imgType;
	ulong sHeight = sourceImg.height;
	ulong sWidth = sourceImg.width;

	float rationX = 1.0f*(sWidth-1)/(width-1);
	float rationY = 1.0f*(sHeight-1)/(height-1);

	Image bilinearInterImg;bilinearInterImg.createImage(width,height,sourceImg_type);

	for(ulong h=0;h<height;h++){
		float y = h*rationY;
		long y1 = static_cast<long>(floor(y)),y2,v = static_cast<long>(y-y1);
		if(y1<(long)sHeight)
			y2 = y1+1;
		else
			y2 = y1;

		float dv = 1.0f - v;

		for(ulong w=0;w<width;w++){

			float x = w*rationX;
			long x1 = static_cast<long>(floor(x)),x2,u = static_cast<long>(x - x1);
			if(x1<(long)sWidth)
				x2 = x1+1;
			else
				x2 = x1;

			float du = 1.0f - u;

			int p,p1,p2,p3,p4;
			p1 = ImageTools::getPixel(sourceImg,y1,x1);
			p2 = ImageTools::getPixel(sourceImg,y2,x1);
			p3 = ImageTools::getPixel(sourceImg,y1,x2);
			p4 = ImageTools::getPixel(sourceImg,y2,x2);

			p = static_cast<int>(du*dv*p1 + du*v*p3 + u*dv*p2 + u*v*p4);

			ImageTools::setPixel(bilinearInterImg,h,w,p);
		}
	}

	if (RUNTIMECHECK)
	{
		string imgName = tools::fileNameFromTime(PATH,"bilinearInterImg",".jpg");
		ImageTools::SaveImage(bilinearInterImg, imgName.c_str());
	}
	return bilinearInterImg;
}


Image ImageProcess::BilateralFilter(const Image& img,int w,double sigmaC,double sigmaS)
{
	sigmaC = sigmaC > 0 ? sigmaC : -sigmaC;
	sigmaS = sigmaS > 0 ? sigmaS : -sigmaS;

	assert(w > 0);

	double scaleC = -0.5/(sigmaC*sigmaC);
	double scaleS = -0.5/(sigmaS*sigmaS);

	//存储最后结果
	Image outImg = img;

	for(ulong y = 0; y < img.height; y++)
	{
		for(ulong x = 0; x < img.width; x++)
		{
			Pixel<double>  centP = img.get<double>(y,x);
			double sum = 0;
			Pixel<double> resultP;
			for(int i = -w; i <= w; i++)
			{
				for (int j = -w; j <= w; j++)
				{
					if((y+i) >= 0 && (y+i) < img.height && (x+j) >= 0 && (x+j) < img.width)
					{
						Pixel<double> tempP = img.get<double>(y+i,x+j);
						double dS = ImageExpandTools::EuclideanDistance(centP,tempP);
						double weightS = exp(dS*dS*scaleS);
						double dC = i*i + j*j;
						double weightC = exp(dC*scaleC);
						sum += weightS * weightC;
						resultP += (tempP*weightC*weightS);
					}
				}
			}
			outImg.set<int>(y,x,resultP/sum);
		}
	}

	string imgName = tools::fileName(PATH,"BilateralFilter_",img.width,img.height,img.imgExt);
	outImg.save(imgName);
	return outImg;
}



Matrix<float> ImageProcess::BilateralFilterWithLab(const Image& img,int w,double sigmaC,double sigmaS)
{
	sigmaC = sigmaC > 0 ? sigmaC : -sigmaC;
	sigmaS = sigmaS > 0 ? sigmaS : -sigmaS;

	assert(w > 0);

	double scaleC = -0.5/(sigmaC*sigmaC);
	double scaleS = -0.5/(sigmaS*sigmaS);

	Matrix<float> LabImg = ColorConversion::RGBToLab(img);
	Matrix<float> labResult(LabImg);

	for(ulong y = 0; y < img.height; y++)
	{
		for(ulong x = 0; x < img.width; x++)
		{
			Pixel<double>  centP;
			if(LabImg.channel > 0) centP.red = LabImg.at(y,x,0);
			if(LabImg.channel > 1) centP.green = LabImg.at(y,x,1);
			if(LabImg.channel > 2) centP.blue = LabImg.at(y,x,2);
			if(LabImg.channel > 3) centP.alpha = LabImg.at(y,x,3);

			double sum = 0;
			Pixel<double> resultP;
			for(int i = -w; i <= w; i++)
			{
				for (int j = -w; j <= w; j++)
				{
					if((y+i) >= 0 && (y+i) < img.height && (x+j) >= 0 && (x+j) < img.width)
					{
						Pixel<double> tempP;
						if(LabImg.channel > 0) tempP.red = LabImg.at(y+i,x+j,0);
						if(LabImg.channel > 1) tempP.green = LabImg.at(y+i,x+j,1);
						if(LabImg.channel > 2) tempP.blue = LabImg.at(y+i,x+j,2);
						if(LabImg.channel > 3) tempP.alpha = LabImg.at(y+i,x+j,3);

						double d = ImageExpandTools::EuclideanDistance(centP,tempP);
						double weightS = exp(d*d*scaleS);
						double dC = i*i + j*j;
						double weightC = exp(dC*scaleC);
						sum += weightS * weightC;
						resultP += (tempP*weightC*weightS);
					}
				}
			}

			resultP = resultP/sum;

			if(labResult.channel > 0) labResult.at(y,x,0) = (float)resultP.red;
			if(labResult.channel > 1) labResult.at(y,x,1) = (float)resultP.green;
			if(labResult.channel > 2) labResult.at(y,x,2) = (float)resultP.blue;
			if(labResult.channel > 3) labResult.at(y,x,3) = (float)resultP.alpha;
		}
	}

	if (!RUNTIMECHECK)
	{
		Image resultImg(img.width,img.height,img.imgType,img.channel);
		ColorConversion::LabtoRGB(labResult,resultImg);
		string imgName = tools::fileName(PATH,"BilateralFilter_",img.width,img.height,img.imgExt);
		resultImg.save(imgName);
	}

	return labResult;
}


void ImageProcess::RecursiveBilateralFilter(const Image& img,Matrix<double>& imgFilter,double sigmaS,double sigmaR)
{
	double range_table[255+1];//compute a lookup table
	double inv_sigma_range;
	inv_sigma_range=1.0/(sigmaR*255);
	for(int i=0;i<=255;i++) range_table[i]=exp(-i*inv_sigma_range);
	double alpha=exp(-sqrt(2.0)/(sigmaS*img.width));//filter kernel size

	Matrix<double> temp(img.height,img.width,img.channel);

	Pixel<BYTE> pixel1,pixel2;
	Pixel<double> pixelT,pixelT2;
	for(ulong y = 0; y < img.height; y++)/*horizontal filtering*/
	{
		pixelT = pixel1 = img.get<BYTE>(y,0);

        temp.at(y,0,0) = pixel1.red;
		temp.at(y,0,1) = pixel1.green;
		temp.at(y,0,2) = pixel1.blue;

		for(ulong x = 1; x < img.width; x++) //from left to right
		{
			pixel2 = img.get<BYTE>(y,x);
			BYTE dr=abs(pixel2.red-pixel1.red);
			BYTE dg=abs(pixel2.green-pixel1.green);
			BYTE db=abs(pixel2.blue-pixel1.blue);
			int range_dist=(((dr<<1)+dg+db)>>2);
			double weight=range_table[range_dist];
			double alpha_=weight*alpha;
			double inv_alpha_=1-alpha_;

			temp.at(y,x,0) = pixelT.red = inv_alpha_ * (double)pixel2.red + alpha_* pixelT.red;
			temp.at(y,x,1) = pixelT.green = inv_alpha_ * (double)pixel2.green + alpha_* pixelT.green;
			temp.at(y,x,2) = pixelT.blue = inv_alpha_ * (double)pixel2.blue + alpha_* pixelT.blue;

			pixel1 = pixel2;
		}
		ulong w1 = img.width-1;

		pixelT = pixel1 = img.get<BYTE>(y,w1);
		temp.at(y,w1,0) = 0.5 * (temp.at(y,w1,0) + pixelT.red);
		temp.at(y,w1,1) = 0.5 * (temp.at(y,w1,1) + pixelT.green);
		temp.at(y,w1,2) = 0.5 * (temp.at(y,w1,2) + pixelT.blue);

		pixelT.green = pixelT.blue = pixelT.red;

		for(int x= (int)img.width - 2; x >= 0; x--) //from right to left
		{
			pixel2 = img.get<BYTE>(y,x);
			BYTE dr=abs(pixel2.red-pixel1.red);
			BYTE dg=abs(pixel2.green-pixel1.green);
			BYTE db=abs(pixel2.blue-pixel1.blue);
			int range_dist=(((dr<<1)+dg+db)>>2);
			double weight=range_table[range_dist];
			double alpha_=weight*alpha;
			double inv_alpha_=1-alpha_;

			pixelT.red = inv_alpha_ * (double)pixel2.red + alpha_* pixelT.red;
			pixelT.green = inv_alpha_ * (double)pixel2.green + alpha_* pixelT.green;
			pixelT.blue = inv_alpha_ * (double)pixel2.blue + alpha_* pixelT.blue;

			temp.at(y,x,0) = 0.5 * (temp.at(y,x,0) + pixelT.red);
			temp.at(y,x,1) = 0.5 * (temp.at(y,x,1) + pixelT.green);
			temp.at(y,x,2) = 0.5 * (temp.at(y,x,2) + pixelT.blue);

			pixel1 = pixel2;
		}
	}
	alpha=exp(-sqrt(2.0)/(sigmaS*img.height));//filter kernel size
	imgFilter = temp;
	for(ulong y = 1; y < img.height; y++)
	{
		for(ulong x = 0; x < img.width; x++)
		{
			pixelT = pixel1 = img.get<BYTE>(y-1,x);
			pixel2 = img.get<BYTE>(y,0);
			BYTE dr=abs(pixel2.red-pixel1.red);
			BYTE dg=abs(pixel2.green-pixel1.green);
			BYTE db=abs(pixel2.blue-pixel1.blue);
			int range_dist=(((dr<<1)+dg+db)>>2);
			double weight=range_table[range_dist];
			double alpha_=weight*alpha;
			double inv_alpha_=1-alpha_;

			for (auto c = 0; c < img.channel; c++)
			{
				temp.at(y,x,c) = inv_alpha_ * temp.at(y,x,c) + alpha_ * temp.at(y-1,x,c);
			}
		}
	}

	for(ulong x = 0; x < img.width; x++)
	{
		for(auto c = 0; c < 3; c++)
		{
			temp.at(temp.row-1,x,c) = 0.5 * (temp.at(temp.row-1,x,c) + imgFilter.at(temp.row-1,x,c));
		}
	}

	for(int y = (int)img.height-2; y >= 0; y--)
	{
		for(int x = 0; x < (int)img.width; x++)
		{
			pixel1 = img.get<BYTE>(y+1,x);
			pixel2 = img.get<BYTE>(y,0);
			BYTE dr=abs(pixel2.red-pixel1.red);
			BYTE dg=abs(pixel2.green-pixel1.green);
			BYTE db=abs(pixel2.blue-pixel1.blue);
			int range_dist=(((dr<<1)+dg+db)>>2);
			double weight=range_table[range_dist];
			double alpha_=weight*alpha;
			double inv_alpha_=1-alpha_;

			for(auto c=0;c<3;c++)
			{
				imgFilter.at(y,x,c) = inv_alpha_ * (double)img.atVal<BYTE>(y,x,c) + alpha_ * imgFilter.at(y+1,x,c);
				temp.at(y,x,c) = 0.5 * (temp.at(y,x,c) + imgFilter.at(y,x,c));
			}
		}
	}

	imgFilter = temp;
}


Image ImageProcess::MeanShift(const Image& img, const float Spatial, const float Color, const unsigned int itre)
{
	const float radiusD2 = Spatial;
	const float radiusR2 = Color;

	const int     N = 20;
	const double  H = 0.499999f;
	const float distR0 = 0.0001f;

	float newL = 0.0f ,  newU = 0.0f , newV = 0.0f;
	float mL   = 0.0f ,  mU   = 0.0f , mV   = 0.0f;

	float RadioR=0.0f;
	int   wSum   = 0    , xSum   = 0   , ySum = 0 ;
	float dif0, dif1, dif2, dif3, dif4, difS;
	float avgL, avgU, avgV;

	Matrix<float> luv(img.height,img.width,img.channel);
	ColorConversion::RGBtoLUV(img,luv);
	Matrix<float> out_luv(img.height,img.width,img.channel);

	for(int winY = 0; winY < (int)luv.row; winY++)
	{
		for(int winX = 0; winX < (int)luv.column; winX++)
		{

			mL = luv.at(winY,winX,0); mU = luv.at(winY,winX,1); mV = luv.at(winY,winX,2);

			int newX = winX;
			int newY = winY;
			unsigned int ITER = 0;
			float distR2 = 10;
			while(distR2 > distR0 && ITER < itre)
			{
				wSum = 0; newL = 0.0f; newU = 0.0f; newV = 0.0f;
				xSum = 0; ySum = 0;
				for(int i = -N; i <= N; i++)
				{
					int r = newY + i;
					if(r < 0 ||r >= (int)luv.row) continue;
					for(int j = -N; j <= N; j++)
					{
						int c = newX + j;
						if(c <0 || c >= (int)luv.column) continue;

						difS=(i*i + j*j)/(radiusD2*radiusD2);
						if(difS < 1.0f)
						{
							dif0 = luv.at(r,c) - mL;
							if(dif0 < 0.0f) dif0 = -1*dif0;

							dif1 = luv.at(r,c,1) - mU; 
							dif2 = luv.at(r,c,2) - mV;

							RadioR = (dif0*dif0 + dif1*dif1 + dif2*dif2)/(radiusR2*radiusR2);

							if(RadioR < 1.0f)
							{
								wSum++; xSum+=j; ySum+=i;
								newL+=luv.at(r,c); newU+=luv.at(r,c,1); newV+=luv.at(r,c,2);
							}
						}//if-else
					} //for j --loop
				} // for i --loop

				avgL = (float)newL/wSum;         avgU = (float)newU/wSum;		 avgV = (float)newV/wSum;

				if(xSum>=0) xSum = (int)((float)xSum/wSum + H);
				else xSum = (int)((float)xSum/wSum - H);
				if(ySum>=0) ySum = (int)((float)ySum/wSum + H);
				else ySum = (int)((float)ySum/wSum - H);

				dif0 = avgL - mL;              dif1 = avgU - mU;
				dif2 = avgV - mV;			  dif3 = (float)(xSum - newX);
				dif4 = (float)(ySum - newY);

				distR2 = (dif0*dif0 + dif1*dif1 + dif2*dif2)+ (dif3*dif3 + dif4*dif4);
				distR2=sqrt(distR2);

				mL=avgL; mU= avgU; mV = avgV;
				newX+= xSum; newY += ySum;
				ITER++;
			}// end while
			out_luv.at(winY,winX)   = mL; 
			out_luv.at(winY,winX,1) = mU;
			out_luv.at(winY,winX,2) = mV;
		}//for--loop
	}//for--loop

	Image destImg(img);
	ColorConversion::LUVtoRGB(out_luv,destImg);
	return destImg;
}