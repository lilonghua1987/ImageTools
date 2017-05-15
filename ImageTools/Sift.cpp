#include "Sift.h"
#include "tools.h"

using namespace std;

Sift::Sift(void)
{
	// SIFT 
	octave = 5;
	sub_level = 6;
	sigma = 1.6;
	k = 2;
	S = 3;
}


Sift::Sift(Image img):img(img)
{
}


Sift::Sift(Image img, double sigma, int k, int cotave, int sub_level, int S)
	:img(img),
	sigma(sigma),
	k(k),
	octave(cotave),
	sub_level(sub_level),
	S(S)
{
}


Sift::~Sift(void)
{
}


void Sift::init(double sigma, int k, int octave, int sub_level, int S)
{
	this->sigma = sigma;
	this->k = k;
	this->octave = octave;
	this->sub_level = sub_level;
	this->S = S;
}


GaussianPyramid Sift::createGaussianPyramid(Image img){

	GaussianPyramid pyramid;

	int img_width = img.width;  //图像宽度
	int img_height = img.height;  //图像高度
	Image tempImg = ImageTools::ImageClone(img);

	for(int i=0;i<sub_level;i++){
		pyramid.factor.push_back(pow(k,1.0*i/S)*sigma);
	}

	for(int o=0;o<octave;o++){
		vector<Image*> po;
		if(o > 0){
			int img_width = pyramid.img[o-1][sub_level-3]->width;  //图像宽度
			int img_height = pyramid.img[o-1][sub_level-3]->height;  //图像高度
			tempImg = ImageTools::ImageClone(*pyramid.img[o-1][sub_level-3]);
			tempImg = ImageProcess::halfSize(tempImg);
		}else{
			//tempImg = GausscianBlur(tempImg,0.5);
			tempImg = ImageProcess::GausscianSeparateBlur(tempImg,0.5);
			tempImg = ImageProcess::doubleSize(tempImg);
		}
		for(int s=0;s<sub_level;s++){
			cout<<"sigma="<<pyramid.factor[s]<<"["<<o+1<<","<<s+1<<"]"<<endl;

			if(s > 0){
				//tempImg = GausscianBlur(tempImg,pyramid.factor[s]);
				tempImg = ImageProcess::GausscianSeparateBlur(tempImg,pyramid.factor[s]);
			}
			Image* copyImg = new Image(ImageTools::ImageClone(tempImg));
			po.push_back(copyImg);

			string imgName = tools::fileName(PATH,"ProduceGaussianPyramid_",img_width,s+1,img.imgExt);

			ImageTools::SaveImage(*copyImg, imgName.c_str());
		}
		pyramid.img.push_back(po);
	}

	return pyramid;
}


void Sift::printGaussianPyramid(GaussianPyramid pyramid){

	for(ulong o=0;o<pyramid.img.size();o++){
		for(ulong s=0;s<pyramid.img[o].size();s++){
			string imgName = tools::fileName(PATH,"ShowGaussianPyramid_",o+1,s+1,img.imgExt);

			cout<<"sigma="<< pyramid.factor[s]<<endl;

			ImageTools::SaveImage(*pyramid.img[o][s], imgName.c_str());

		}
	}
}


DOGPyramid Sift::createDOGPyramid(GaussianPyramid pyramid){

	DOGPyramid dogPyramid;

	FILE* fp4;
	fopen_s(&fp4, "temp\\dog_pixel.txt", "w"); 
	fprintf(fp4,"\tabovePixel：");
	fprintf(fp4,"\tunderPixel：");
	fprintf(fp4,"\tabsPixel：\n");

	for(ulong o=0;o<pyramid.img.size();o++){
		vector<Image*> dogO;
		FREE_IMAGE_TYPE img_type =  pyramid.img[o][0]->imgType;
		int bpp = pyramid.img[o][0]->channel;
		int img_width = pyramid.img[o][0]->width;  //图像宽度
		int img_height = pyramid.img[o][0]->height;  //图像高度
		for(ulong s=1;s<pyramid.img[o].size();s++){
			Image *dogImg = pyramid.img[o][s];
			Image *underImg = pyramid.img[o][s-1];
			if((img_type == FIT_BITMAP) && (bpp == 1)) {
				for(int y = 0; y < img_height; y++) {
					for(int x = 0; x < img_width; x++) {
						int bitsAbove = ImageTools::getPixel(*dogImg,y,x);
						int bitsUnder = ImageTools::getPixel(*underImg,y,x);
						if(o<1 && s<2){
							fprintf(fp4,"\t%10d",bitsAbove);
							fprintf(fp4,"\t%10d",bitsUnder);
						}

						/*if((bitsAbove[x] - bitsUnder[x]) < 0)
						bitsAbove[x] = 0;
						else*/
						bitsAbove = bitsAbove - bitsUnder;
						ImageTools::setPixel(*dogImg,y,x,bitsAbove);
						if(o<1 && s<2){
							fprintf(fp4,"\t%10d\n",bitsAbove);
						}
					}
				}
			}
			dogO.push_back(dogImg);

			string imgName = tools::fileName(PATH,"ProduceDOGPyramid_",img_width,s,img.imgExt);
			ImageTools::SaveImage(*dogImg, imgName.c_str());
		}
		dogPyramid.img.push_back(dogO);
	}

	fclose(fp4);

	return dogPyramid;
}


void Sift::printDOGPyramid(DOGPyramid dog){

	for(ulong o=0;o<dog.img.size();o++){
		for(ulong s=0;s<dog.img[o].size();s++){
			string imgName = tools::fileName(PATH,"ShowDOGPyramid_",o+1,s+1,img.imgExt);

			ImageTools::SaveImage(*dog.img[o][s], imgName.c_str());
		}
	}
}


ExtremeMatrix Sift::createExtremeMatrix(DOGPyramid dog){

	ExtremeMatrix matrix;
	long int sum=0;
	for(ulong o=0;o<dog.img.size();o++){
		FREE_IMAGE_TYPE img_type =  dog.img[o][0]->imgType;
		int bpp = dog.img[o][0]->channel;
		ulong img_width = dog.img[o][0]->width;  //图像宽度
		ulong img_height = dog.img[o][0]->height;  //图像高度
		ExtremePointMatrix pointMatrix;
		for(ulong s=1;s<dog.img[o].size()-1;s++){
			cout<<"获取第:["<<o+1<<","<<s<<"]"<<endl;
			vector<FeaturePoint> point;
			Image *imgAbove = dog.img[o][s+1];
			Image *imgMid = dog.img[o][s];
			Image *imgUnder = dog.img[o][s-1];
			if((img_type == FIT_BITMAP) && (bpp == 8)) {
				for(ulong y = 1; y < img_height-1; y++) {
					for(ulong x = 1; x < img_width-1; x++) {
						FeaturePoint feature;
						int bitsAbove=0,bitsMid =0,bitsUnder=0,bits=ImageTools::getPixel(*imgMid, y,x);
						vector<int> dot;
						for(ulong i=y-1;i<=y+1;i++){
							for(ulong j=x-1;j<=x+1;j++){
								bitsAbove = ImageTools::getPixel(*imgAbove, i,j);
								bitsMid = ImageTools::getPixel(*imgMid, i,j);
								bitsUnder = ImageTools::getPixel(*imgUnder, i,j);
								dot.push_back(bitsAbove);
								dot.push_back(bitsMid);
								dot.push_back(bitsUnder);
							}
						}
						sort(dot.begin(), dot.end() );

						if(bits == dot[0]){
							feature.indexO = o;
							feature.indexS = s;
							feature.indexX = x;
							feature.indexY = y;
							feature.pixel = bits;
							feature.min = true;
							point.push_back(feature);
							sum++;
							//cout<<"bits="<<bits<<endl;
						}else if(bits == dot[dot.size()-1]){
							feature.indexO = o;
							feature.indexS = s;
							feature.indexX = x;
							feature.indexY = y;
							feature.pixel = bits;
							feature.min = false;
							point.push_back(feature);
							sum++;
							//cout<<"bits="<<bits<<endl;
						}

					}
				}
			}
			pointMatrix.layer.push_back(point);
		}
		matrix.group.push_back(pointMatrix);
	}
	cout<<"sum="<<sum<<endl;

	return matrix;
}



void Sift::printExtremeMatrix(ExtremeMatrix extreamMatrix){

	FILE* fp5;
	fopen_s(&fp5, "temp\\extreamMatrix.txt", "w"); 
	fprintf(fp5,"\tIndexX：");
	fprintf(fp5,"\tIndexY：");
	fprintf(fp5,"\tPixel：");
	fprintf(fp5,"\tmin：\n");

	for(ulong o=0;o<extreamMatrix.group.size();o++){
		ExtremePointMatrix  pointMatrix  = extreamMatrix.group[o];
		for(ulong s=0;s<pointMatrix.layer.size();s++){
			for(ulong i=0;i<pointMatrix.layer[s].size();i++){
				FeaturePoint f = pointMatrix.layer[s][i];
				fprintf(fp5,"\t%10d",f.indexX);
				fprintf(fp5,"\t%10d",f.indexY);
				fprintf(fp5,"\t%5f",f.pixel);
				fprintf(fp5,"\t%5d\n",f.min);
			}
		}
	}
	fclose(fp5);
}