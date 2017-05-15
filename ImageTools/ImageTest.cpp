#include "ImageProcess.h"
#include "Fouier.h"
#include "Harris.h"
#include "sharpen.h"
#include "AdaptiveWeight.h"
#include "Sift.h"
#include "tools.h"
#include "ImageSegment.h"
#include "BPStereo.h"
#include "SLIC.h"
#include "Reconstruction.h"
#include "LFrame.h"
//#include <Eigen/Dense>
#include "StereoMatch.h"
#include"MrfStereo.h"
#include "PSO/PSOStereoMatch.h"
#include "Elas.h"
#include "Algorithm.h"
#include "CSVFile.h"


//memory leak check
//#ifdef  _DEBUG
//# include <vld.h>
//#endif



//#include <thread>

#define SOURCE "img\\"

GLUquadricObj * drwCylindrical()
{
	GLUquadricObj *quadratic;					// 二次几何体
    GLuint  object=0;						// 二次几何体标示符

	quadratic=gluNewQuadric();				// 创建二次几何体
	gluQuadricNormals(quadratic, GLU_SMOOTH);		// 使用平滑法线
	gluQuadricTexture(quadratic, GL_TRUE);		// 使用纹理

	return quadratic;
}

int main(int argc, char** argv)
{	
	string imgLeft,imgRight;
	string dir("Cones/");
	imgLeft.append(SOURCE);
	imgLeft.append(dir);
	imgLeft.append("imL.png");
	imgRight.append(SOURCE);
	imgRight.append(dir);
	imgRight.append("imR.png");

	Image exifImg(imgLeft);

	cout << "Make :" << exifImg.hInfo.getMake() << " Model :" << exifImg.hInfo.getModel() << " Focal : " << exifImg.hInfo.getFocal() << endl;

	FITAG* tagFocal = nullptr;
	string tF;
	FreeImage_GetMetadata(FIMD_EXIF_EXIF, exifImg.data, "FocalLength", &tagFocal);
	if (tagFocal)
	{
		tF = FreeImage_TagToString(FIMD_EXIF_EXIF, tagFocal);
		cout << tF << endl;
	}
	size_t index = tF.find_first_of(' ');
	cout << tF.substr(0,index) << endl;

	FITAG *tag = NULL;
	FIMETADATA *mdhandle = NULL;

	mdhandle = FreeImage_FindFirstMetadata(FIMD_EXIF_EXIF, exifImg.data, &tag);

	if (mdhandle) {
		do {
			// get the tag key
			const char *key = FreeImage_GetTagKey(tag);
			// convert the tag value to a string
			const char *value = FreeImage_TagToString(FIMD_EXIF_EXIF, tag);

			// print the tag 
			// note that most tags do not have a description, 
			// especially when the metadata specifications are not available
			cout << FreeImage_GetTagKey(tag) << "=" << value << "\n";

		} while (FreeImage_FindNextMetadata(mdhandle, &tag));
	}

	FreeImage_FindCloseMetadata(mdhandle);

		


	/*Image img(imgLeft);
	DWORD hist[256] = {0};
	ImageTools::getHistoGram(img,hist,FICC_RGB);*/

	/*RunTimer rt;
	rt.start();
	PSOStereoMatch pMatch(imgLeft,imgRight,-19);
	pMatch.RunMatch();
	rt.timeDisplay("RunMatch");*/

	/*int dispNu = 15;
	int scale = 16;
	double eps = 0.01 * 0.01;

	Image gImg("E:/vs2012workspace/KZ2/result/tsukuba/imL.png");
	Image dImg("E:/vs2012workspace/KZ2/result/tsukuba/disp.png");

	Matrix<double> img(gImg.height, gImg.width, gImg.channel);
	ImageExpandTools::Image2Mat(gImg, img);
	img /= 255.0;

	Matrix<double> disp(dImg.height, dImg.width, dImg.channel);
	ImageExpandTools::Image2Mat(dImg, disp);
	disp /= (double)scale;
	uint r = (uint)std::ceil(std::max(img.row/40.0,img.column/40.0));

	disp = Filter::weightedMedianFilter(disp,img,r,dispNu,eps);
	disp = Filter::medianFilter(disp);

	Image dOut = ImageExpandTools::Mat2Image(disp,scale);

	dOut.save("E:/vs2012workspace/KZ2/result/tsukuba/wmDisp.png");*/

	//RunTimer rt;
	//rt.start();
	////runStereoMatch(imgLeft, imgRight, "temp/teddy/");
	////MrfStereoMatch(imgLeft,imgRight);
	//ElasStereoMatch(imgLeft, imgRight);
	//rt.timeDisplay("StereoMatch");

	/*Matrix<float> luv = ColorConversion::RGBToLUV(imgLeft);
	Image luvImg(luv.column,luv.row,FIT_BITMAP,luv.channel);
	ColorConversion::LUVtoRGB(luv,luvImg);
	luvImg.save("temp/luvimg");*/

	/*Matrix<float> lab = ColorConversion::RGBToLab(imgLeft);
	Image labImg(lab.column,lab.row,FIT_BITMAP,lab.channel);
	ColorConversion::LabtoRGB(lab,labImg);
	labImg.save("temp/labImg");*/

	/*Image meanshift = ImageProcess::MeanShift(imgLeft, 20, 30);
	meanshift.save("temp/meanshift");*/

	//Image imL(imgLeft);
	//ImageTools::ImageAdjustBrightness(imL,-20);
	//imL.save("temp/view1");
	//Image imR(imgRight);
	//ImageTools::ImageAdjustBrightness(imR,20);
	//imR.save("temp/view5");

	/*Image ncc = ImageProcess::NCC(ImageTools::ImageToGrey(imgLeft),ImageTools::ImageToGrey(imgRight),25,20);
	ImageTools::SaveImage(ncc,"temp/ncc.jpg");*/

	//ImageExpandTools::build3d(Image("img/tree_left_d.jpg"),"img/tree_left.jpg",120.f,"temp/stereo.ply");

	//darkChannel(imgLeft);

	//StereoMatch str(imgLeft,imgRight);

	/*Image img(imgLeft);
	ImageExpandTools::ImageAdjustContrast(img,150.0);
	ImageTools::SaveImage(img,"temp/plane_s.bmp");*/

	//Reconstruction::shading(ImageExpandTools::ImageToGrey(imgRight),20,39.7,92.3);

	/*Matrix<float> m = ImageExpandTools::loadFromFile("img/plane_psp_00_image.txt");
	ImageTools::SaveImage(ImageExpandTools::Mat2Image(m,true),"temp//plane.jpg");*/

	//LFrame fram(1024,768,"李龙华");
	//fram.draw();
	////fram.~LFrame();

	//LFrame fram2(600,480,"test");
	//fram2.draw();


	/*Image img(imgLeft.c_str());
	Image erImg = ImageExpandTools::Threshold(img,220);
	Image saveImage = ImageExpandTools::InPaint(img,erImg,25);
	ImageTools::SaveImage(saveImage,"temp//InPaint.jpg");*/
	//ImageExpandTools::inPaintAllImages("E:\\百度云\\计算机视觉代码\\ImageTools\\ImageTools\\images",".jpg","inpaint");

	//readBmp("img//00003200.bmp");

	/*Image img("img//00003200.bmp");
	Image saveImg = ImageExpandTools::Image2Standard(img,10);
	ImageTools::SaveImage(saveImg,"temp//saveImg1.bmp");*/

	//vector<string>* lf = new vector<string>();
	//string fileExtension(".bmp"),saveDir("F:\\images");
	//tools::listFileByDir("F:\\新建文件夹",fileExtension,lf);

	//tools::CreatDir(saveDir.c_str());

	//for (int i = 1; i < lf->size(); i++)
	//{	
	//	/*std::thread t(threadReSizeImage,lf->at(i));*/
	//	Image img(lf->at(0).c_str());
	//    Image saveImg = ImageExpandTools::Image2Standard(img,10);
	//	ImageTools::SaveImage(saveImg,tools::fileNameFromTime(saveDir,string(),fileExtension).c_str());	
	//	/*t.join();*/					
	//}

	//delete lf;

	//ImageExpandTools::reSizeImage("E:\\vs2012workspace\\MotionDetection\\Release\\temp\\images",".jpg","images",Point2D<int>(0,100),Point2D<int>(1280,1024));

	/*vector<string>* lf = new vector<string>();
	string fileExtension(".jpg"),saveDir("images");
	tools::listFileByDir("E:\\vs2012workspace\\MotionDetection\\Release\\temp\\images",fileExtension,lf);

	tools::CreatDir(saveDir.c_str());
	Point2D<int> top(0,100), bottom(1280,1024);

	for (int i = 0; i < lf->size(); i++)
	{	
	std::thread t(threadReSizeImage,lf->at(i));
	string name = "images\\";
	name.append(tools::getFileName(lf->at(i)));
	Image img(lf->at(i).c_str());
	Image saveImg = ImageTools::ImageCopy(img,top,bottom);
	ImageTools::SaveImage(saveImg,name.c_str());	
	t.join();		
	_sleep(1000);
	}

	delete lf;*/

	//int numlabels(0);
	//int m_spcount = 500,m_compactness = 10;
	//Image img(imgLeft);
	//std::vector<int> labels;
	//SLIC slic;
	//slic.PerformSLICO_ForGivenK(img, labels, numlabels, m_spcount, m_compactness);//for a given number K of superpixels
	//////slic.PerformSLICO_ForGivenStepSize(img, labels, numlabels, m_stepsize, m_compactness);//for a given grid step size
	//////slic.DrawContoursAroundSegments(img, labels, width, height, 0);//for black contours around superpixels
	//slic.DrawContoursAroundSegmentsTwoColors(img, labels);//for black-and-white contours around superpixels
	//////slic.SaveSuperpixelLabels(labels,width,height,picvec[k],saveLocation);
	//ImageTools::SaveImage(img,"temp//slic_new.jpg");

	/*BPStereo bp(50);
	Image result = bp.bpMatch(ImageTools::ImageToGrey(imgLeft),ImageTools::ImageToGrey(imgRight));
	ImageTools::SaveImage(result,"temp/bp.jpg");*/
	/*Image result = bp.restore(ImageTools::ImageToGrey(imgLeft.c_str()));
	ImageTools::SaveImage(result,"temp//restore.jpg");*/

	/*Image he = ImageProcess::HistogramEqualization(ImageExpandTools::ImageToGrey(imgLeft.c_str()));
	ImageTools::SaveImage(he,"temp//HistogramEqualization.jpg");*/

	/*sharpen sharp;
	Image img = sharp.laplacian(imgRight.c_str(),8);
	ImageTools::SaveImage(img,"temp//laplacian.jpg");*/

	/*Image img = ImageProcess::NCC(ImageExpandTools::ImageToGrey(imgRight.c_str()),ImageExpandTools::ImageToGrey(imgLeft.c_str()),21,21,19);
	ImageTools::SaveImage(img,"temp//ncc.jpg");*/

	//Image ImageRotateEx = ImageExpandTools::ImageRotateEx(imgL,30,Point2D<double>(30,60),Point2D<double>(imgL.width/2,imgL.height/2),false);

	//string imgName = tools::fileName(PATH,"ImageRotateEx",ImageRotateEx.width,ImageRotateEx.height,ImageRotateEx.imgExt);
	//ImageTools::SaveImage(ImageRotateEx, imgName.c_str());

	/*Image imgR(imgRight.c_str());
	Image greyImgR = ImageExpandTools::ImageToGrey(imgR);*/

	//AdaptiveWeight a = AdaptiveWeight(Image(imgLeft.c_str()),Image(imgRight.c_str()));
	/*Image outImg = ImageProcess::GausscianSeparateBlur(imgL,1.6);*/
	/*Image bilinearInterImg = ImageProcess::BilinearInterpolation(greyImg,imgL.height*2.5,imgL.width*2.5);*/

	//sharpen sharp =sharpen();
	/*sharp.Canny(ImageExpandTools::ImageToGrey(imgLeft.c_str()),0.1,0.9,0.76);*/

	/*double *grad = sharp.laplacianGrad(greyImg);
	sharp.printLaplacianGrad(grad,greyImg);
	sharp.laplacian(ImageProcess::GausscianSeparateBlur(greyImg,1.6),4);*/

	//ImageProcess::halfSize(greyImg);
	/*ImageProcess::doubleSize(greyImg);*/

	/*Harris harris = Harris(ImageExpandTools::ImageToGrey(imgLeft.c_str()));
	double **matrix = harris.difference(ImageExpandTools::ImageToGrey(imgLeft.c_str()));
	AngularPoint **angular = harris.checkAngular(matrix,3,5000);
	harris.printAngular(angular);
	vector<vector<AngularPoint>> point = harris.getAngular(angular);
	harris.printAngular(point,imgLeft.c_str());*/

	/*ParallaxMatrix mm = harris.NCCByAngular(greyImg,greyImgR,21,21,70);
	ImageProcess::printMapingPointMatrix(mm,greyImg.imgType,greyImg.fif);*/

	//ImageProcess::imageReverse(greyImg);

	/*GreyHistogram greyHistogram = ImageProcess::getGreyHistogram(ImageExpandTools::ImageToGrey(imgLeft.c_str()));
	ImageProcess::printHistogram(greyHistogram);*/

	/*ImageProcess::GausscianBlur(greyImg,1.6);
	ImageProcess::GausscianSeparateBlur(imgL,1.6);*/

	/*Sift sift;
	GaussianPyramid pyramid = sift.createGaussianPyramid(greyImg);
	sift.printGaussianPyramid(pyramid);
	DOGPyramid dog = sift.createDOGPyramid(pyramid);
	sift.printDOGPyramid(dog);
	ExtremeMatrix matrix = sift.createExtremeMatrix(dog);
	sift.printExtremeMatrix(matrix);*/

	/*ImageProcess::fourier(greyImg);*/

	/*ImageProcess::parseDFourier(greyImg,greyImgR);*/

	/*ImageProcess::parseFourier(greyImg,greyImgR);*/

	//ParallaxMatrix mm = ImageProcess::NCC(greyImg,greyImgR,21,21,100);
	//cout<<"The NCC is end,the print is begine"<<endl;
	//ImageProcess::printMapingPointMatrix(mm,greyImg.imgType,greyImg.fif);

	/*RunTimer t;t.start();
	ImageProcess::BilateralFilterWithLab(imgLeft,5,2,30);
	t.timeDisplay("BilateralFilterLab");*/

	/*RunTimer t;t.start();
	printf("processing\n");
	int num_ccs;
	ImageSegment(imgLeft.c_str(), 1.5, 500, 25, num_ccs); 
	printf("got %d components\n", num_ccs);
	t.timeDisplay("ImageSegment");*/

    /*int a[22] = {1,6,2,8,4,3,11,7,44,77,22,11,33,10,99,100,77,66,55,44,22,0};
	alg::quickSort(a,22,false);
	alg::shellSort(a,22,false);*/

	std::system("pause");

	return 0;
}