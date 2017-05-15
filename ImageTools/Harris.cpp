#include "Harris.h"
#include "tools.h"

Harris::Harris(void)
{
}


Harris::~Harris(void)
{
}


Harris::Harris(Image img)
{
	_IMG_WIDTH_ = img.width;
	_IMG_HEIGHT_ = img.height;


}

double** Harris::difference(Image img){

	const int img_width = img.width;  //图像宽度
	const int img_height = img.height;  //图像高度

	double **relationMatrix =(double **)malloc(sizeof(double *)*3);

	double *mX = (double *)malloc(sizeof(double)*img_height*img_width);

	double *mY = (double *)malloc(sizeof(double)*img_height*img_width);

	double *mZ = (double *)malloc(sizeof(double)*img_height*img_width);

	Image differenceXImg = ImageTools::ImageCopy(img,0,0,img_width,img_height);
    Image differenceYImg = ImageTools::ImageCopy(img,0,0,img_width,img_height);

	for(long h=0;h<img_height;h++){
		for(long w=0;w<img_width;w++){
			int valueX = 0, valueY = 0, valueZ = 0;

			if(w<img_width-1){
				valueX = ImageTools::getPixel(img,h,w) - ImageTools::getPixel(img,h,w+1);
			}else{
				valueX = ImageTools::getPixel(img,h,w) - ImageTools::getPixel(img,h,w-1);
			}
			mX[h*img_width+w] = valueX * valueX;

			ImageTools::setPixel(differenceXImg,h,w,valueX);

			if(h<img_height-1){
				valueY = ImageTools::getPixel(img,h,w) - ImageTools::getPixel(img,h+1,w);
			}else{
				valueY = ImageTools::getPixel(img,h,w) - ImageTools::getPixel(img,h-1,w);
			}
			mY[h*img_width+w] = valueY * valueY;

			ImageTools::setPixel(differenceYImg,h,w,valueY);

			valueZ = valueX * valueY;

			mZ[h*img_width+w] = valueZ;

		}
	}

	relationMatrix[0] = mX;
	relationMatrix[1] = mY;
	relationMatrix[2] = mZ;
	
	string imgXName = tools::fileNameFromTime(PATH,"differenceXImg",".jpg");
	string imgYName = tools::fileNameFromTime(PATH,"differenceYImg",".jpg");
	ImageTools::SaveImage(differenceXImg, imgXName.c_str());
	ImageTools::SaveImage(differenceYImg, imgYName.c_str());

	return relationMatrix;
}


double* Harris::GausscianBlur2D(double *matrixData,int height,int width,int window,double sigma){

	int range = (int)floor(fabs(sigma)*3)*2 + 1;

	if(window>range)
		range = window;

	int mid = range/2;

	MatrixMat matrix;

	double gauss = 1.0/(2*PI*sigma*sigma),sum=0;;

	matrix.height = matrix.width = range;

	cout<<"Gauss Matrix height:"<<matrix.height<<",sigma="<<sigma<<endl;

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

	//归一化,确保高斯权值在[0,1]之间  
	for(int i=0;i<matrix.height;i++){
		for(int j=0;j<matrix.width;j++){
			matrix.data[i][j]/=sum;
		}
	}

	double *tempM =new double[height*width];

	memcpy(tempM, matrixData, sizeof(double) * width * height);

	for(int y = 0; y < height; y++) {
		for(int x = 0; x < width; x++) {
			double tempPix = 0,m=0;
			for(int i=0;i<range;i++){
				for(int j=0;j<range;j++){
					int changeX = j-mid+1+x;
					int changeY = i-mid+1+y;

					int flageX = 0;
					int flageY = 0;

					if(changeX>=0 && changeX<width){
						flageX = 1;
					}
					if(changeY>=0 && changeY<height){
						flageY = 1;
					}


					if(flageX && flageY){
						tempPix+=tempM[changeY*width+changeX]*matrix.data[i][j];
					}
				}
			}
			matrixData[y*width+x] = tempPix;
		}
	}

	return matrixData;
}


double* Harris::GausscianSeparateBlur(double *matrixData,int height,int width,int window,double sigma){

	int ksize = (int)floor(fabs(sigma) * 3) * 2 + 1;

	ksize = window > ksize ? window : ksize; 

	//计算一维高斯核
	double *kernel = new double[ksize];

	double scale = -0.5/(sigma*sigma);
	double cons = 1/sqrt(-scale / PI);

	double sum = 0;
	int kcenter = ksize/2;
	int i = 0, j = 0;
	for(i = 0; i < ksize; i++)
	{
		int x = i - kcenter;
		*(kernel+i) = cons * exp(x * x * scale);//一维高斯函数
		sum += *(kernel+i);
	}

	//归一化,确保高斯权值在[0,1]之间
	for(i = 0; i < ksize; i++)
	{
		*(kernel+i) /= sum;
	}

	double *tempM =new double[height*width];

	memcpy(tempM, matrixData, sizeof(double) * width * height);

	//x方向一维高斯模糊
	for(int y = 0; y < height; y++){
		for(int x = 0; x < width; x++){
			double mul = 0;
			sum = 0;
			double bmul = 0, gmul = 0, rmul = 0;
			for(i = -kcenter; i <= kcenter; i++)
			{
				if((x+i) >= 0 && (x+i) < width)
				{	
					mul += tempM[y*width+x+i]*(*(kernel+kcenter+i));
					sum += (*(kernel+kcenter+i));
				}
			}
			
		    matrixData[y*width+x] = mul/sum;
			
		}
	}

	memcpy(tempM, matrixData, sizeof(double) * width * height);
	
	//y方向一维高斯模糊
	for(int x = 0; x < width; x++)
	{
		for(int y = 0; y < height; y++)
		{
			double mul = 0;
			sum = 0;
			double bmul = 0, gmul = 0, rmul = 0;
			for(i = -kcenter; i <= kcenter; i++)
			{
				if((y+i) >= 0 && (y+i) < height)
				{
					mul += tempM[(y+i)*width+x]*(*(kernel+kcenter+i));
					sum += (*(kernel+kcenter+i));
				}
			}
			
			matrixData[y*width+x] = mul/sum;
		
		}
	}

	return matrixData;
}


AngularPoint** Harris::checkAngular(double** matrixData , int size , double thresh){

	AngularPoint **angular;

	angular=(AngularPoint **)malloc(sizeof(AngularPoint *)*_IMG_HEIGHT_);

	for(int i=0;i<3;i++){
		GausscianSeparateBlur(matrixData[i],_IMG_HEIGHT_,_IMG_WIDTH_,0,0.8);
	}

	for(int h=0;h<_IMG_HEIGHT_;h++){
		angular[h]=(AngularPoint *)malloc(sizeof(AngularPoint)*_IMG_WIDTH_);
		for(int w=0;w<_IMG_WIDTH_;w++){
			double A = matrixData[0][h*_IMG_WIDTH_+w];
			double B =  matrixData[1][h*_IMG_WIDTH_+w];
			double C = matrixData[2][h*_IMG_WIDTH_+w] * matrixData[2][h*_IMG_WIDTH_+w];
			
			/*double r =(r1-r2)*(r1-r2) - k*r3*r3;*/
			double R = 0;

			if(A>0 || B>0)
				R = fabs((A*B - C*C)/(A+B));

			angular[h][w].r = R;
			angular[h][w].y = h;
			angular[h][w].x = w;
			angular[h][w].isMax = false;

			if(R>thresh){
				angular[h][w].isAngular = true;
			}else{
				angular[h][w].isAngular = false;
			}
		}
	}

	if(size%2 == 0)
		size+=1;

	for(int h=size/2;h<_IMG_HEIGHT_-size/2;h++){
		for(int w=size/2;w<_IMG_WIDTH_-size/2;w++){
			double maxR = angular[h-size/2][w-size/2].r;
			int indexX = w-size/2,indexY = h-size/2;
			for(int s=-size/2;s<=size/2;s++){
				if(maxR<angular[h+s][w+s].r){
					maxR = angular[h+s][w+s].r ;
					indexX = w+s;
					indexY = h+s;
				}
			}

			angular[indexY][indexX].isMax = true;
		}
	}

	return angular;
}


void Harris::printAngular(AngularPoint** angular){

	FILE* fp1;
	fopen_s(&fp1, "temp\\angular.txt", "w"); 

	for(int h=0;h<_IMG_HEIGHT_;h++){
		for(int w=0;w<_IMG_WIDTH_;w++){
			if(angular[h][w].isMax && angular[h][w].isAngular){
				AngularPoint point = angular[h][w];
				fprintf(fp1,"[%10f,%d,%d]",point.r,point.y,point.x);
			}
		}
			
		fprintf(fp1,"\n");
	}

	fclose(fp1);
}


vector<vector<AngularPoint>> Harris::getAngular(AngularPoint **angular){
	
	vector<vector<AngularPoint>> angularPoint;

	for(int h=0;h<_IMG_HEIGHT_;h++){
		vector<AngularPoint> angularPointLayer;
		for(int w=0;w<_IMG_WIDTH_;w++){
			if(angular[h][w].isMax && angular[h][w].isAngular){
				angularPointLayer.push_back(angular[h][w]);
			}
		}
		angularPoint.push_back(angularPointLayer);
	}

	return angularPoint;
}


void Harris::printAngular(vector<vector<AngularPoint>> angularPoint,Image img){

	Image angularImg = ImageTools::ImageCopy(img,0,0,_IMG_WIDTH_,_IMG_HEIGHT_);
	if(!angularImg.data) return;

	/*FILE* fp2;
	fopen_s(&fp2, "temp\\AngularPoint.txt", "w");*/ 

	for(size_t i=0;i<angularPoint.size();i++){
		for(size_t j=0;j<angularPoint.at(i).size();j++){
			AngularPoint point = angularPoint.at(i).at(j);
			ImageTools::setColorPixel(angularImg,point.y,point.x,Pixel<uchar>(255,255,255));
		    /*fprintf(fp2,"[%10f,%d,%d]",point.r,point.y,point.x);*/
		}

		/*fprintf(fp2,"\n");*/
	}

	/*fclose(fp2);*/

	string imgName;
	imgName.append(PATH);
	imgName.append("angularImg");
	imgName.append(".jpg");
	ImageTools::SaveImage(angularImg, imgName.c_str());
}


ParallaxMatrix Harris::NCCByAngular(Image referenceImg,Image targetImg,int H,int W,int dispMax){

	double **matrix = difference(referenceImg);

	AngularPoint **angular = checkAngular(matrix,3,500);

	printAngular(angular);

	vector<vector<AngularPoint>> angularPoint = getAngular(angular);

	printAngular(angularPoint,referenceImg);

	FREE_IMAGE_TYPE searchImg_type =  targetImg.imgType;  //图像类型
	int searchImg_width = targetImg.width;  //图像宽度
	int searchImg_height = targetImg.height;  //图像高度

	ParallaxMatrix pm;

	int midH = (H-1)/2,midW = (W-1)/2;

	for(int y1=midH;y1<searchImg_height-midH;y1++){
		vector<ParallaxPoint> pl;
		int sizeX = angularPoint.at(y1).size(),maxR = 0;
		for(int x1=midW;x1<searchImg_width-midW;x1++){
			ParallaxPoint point;
			point.relation = 0;

			bool pass = false;
			if(sizeX > x1){
				pass = angularPoint.at(y1).at(x1).x<(searchImg_width-midW)?true:false;
			}

			if(pass){
				int index = angularPoint.at(y1).at(x1).x;
				int dispRange = index-dispMax;
				if(dispRange <midW)
					dispRange = midW;
				for(;dispRange<=index;dispRange++){
					double r=0,r1=0,r2=0,r3=0;
					for(int h=-midH;h<=midH;h++){
						for(int w=-midW;w<=midW;w++){
							int a,b;
							a=ImageTools::getPixel(referenceImg,y1+h,index+w);
							b=ImageTools::getPixel(targetImg,y1+h,dispRange+w);
							r1 += a*a;
							r2 += b*b;
							r3 += a*b;
						}
					}

					r1 =sqrt(r1);
					r2 = sqrt(r2);
					r = r3/(r2*r1);

					if(point.relation<r){
						point.relation = r;
						point.parallax = index-dispRange;

					}

					if(maxR<point.parallax)
						maxR = point.parallax ;
				}
			}else
			{
				point.parallax = maxR/2;
			}

			pl.push_back(point);
		}
		pm.value.push_back(pl);
	}

	pm.height = searchImg_height - H;
	pm.width = searchImg_width - W;

	return pm;
}