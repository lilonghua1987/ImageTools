#pragma once

#define _CRT_SECURE_NO_WARNINGS

#include <string>
#include <time.h>
#include <fstream>
#include <sstream>
#include <io.h> 
#include <vector>
#include <iostream>
#include "CoreStruct.h"

#ifndef PI
#define PI 3.14159265358979323846264338327950288419716939937510582
#endif

#ifdef _WIN32
#include <direct.h>
#include <io.h>
#elif _LINUX
#include <stdarg.h>
#include <sys/stat.h>
#endif

#ifdef _WIN32
#define ACCESS _access
#define MKDIR(a) _mkdir((a))
#elif _LINUX
#define ACCESS access
#define MKDIR(a) mkdir((a),0755)
#endif

using namespace std;

/*
计算运行时的时间
*/
class   RunTimer		
{
public: 
	void start();
	float stop();
	void timeDisplay(std::string disp);
	void fpsDisplay(std::string disp);
private:
	clock_t m_begin; 
	clock_t m_end;
};


/************************************************************************/
/*   一些工具方法                                                                   */
/************************************************************************/
class tools
{
public:
	tools(void);
	~tools(void);

	/************************************************************************/
	/* @prefix   文件名字                                                                   */
	/* @suffix 文件的存储格式        example:.jpg                                           */
	/************************************************************************/
	 static string fileNameFromTime(string prefix,string suffix);
	 static string fileNameFromTime(const char* prefix,const char* suffix);
	 static string fileNameFromTime(string path,string name,string suffix);
	 static string fileNameFromTime(const char* path,const char* name,const char* suffix);
	 static string fileName(string prefix,long w,long h,string suffix);
	 static string fileName(const char* prefix,long w,long h,const char* suffix);
	 static string fileName(string path,string name,long w,long h,string suffix);
	 static string fileName(const char* path,const char* name,long w,long h,const char* suffix);
	 static string getFileName(const string name);

	 static void getWords(std::vector<std::string>& words, const std::string& str, const char splite);

     static long bound(const long x, const long min, const long max)
	 {
		  return (x < min ? min : (x > max ? max : x));
     };

	 static long lowBound(const long x, const long min)
	 {
		 return (x < min) ? min:x;
	 }

	 static long upBound(const long x, const long max)
	 {
		 return (x > max) ? max:x;
	 }

	 static Point2D<long> bound(const Point2D<long> x, const Point2D<long> min, const Point2D<long> max)
	 {
		 return upBound(lowBound(x,min),max);
	 }

	 static Point2D<long> lowBound(const Point2D<long> point, const Point2D<long> pointMin)
	 {
		 return Point2D<long>(lowBound(point.x,pointMin.x),lowBound(point.y,pointMin.y));
	 }

	 static Point2D<long> upBound(const Point2D<long> point, const Point2D<long> pointMax)
	 {
		 return Point2D<long>(upBound(point.x,pointMax.x),upBound(point.y,pointMax.y));
	 }

	 template<typename T>
	 static T Max(const T &x, const T &y)
	 {
		  return (x < y ? y : x);
     };

	 template<typename T>
	 static T Min(const T &x, const T &y)
	 {
		  return (x > y ? y : x);
     };

	 template<typename T>
	 static int Round(const T& x)
	 {
		 return ((int)( (x) >= 0 ? (x) + .5 : (x) - .5));
	 };

	 template<typename T1,typename T2>
	 static double EuclideanDistance(const Pixel<T1> &p1, const Pixel<T2> &p2)
	 {
		 double d1 = p1.red - p2.red;
		 double d2 = p1.green - p2.green;
		 double d3 = p1.blue - p2.blue;
		 double d4 = p1.alpha - p2.alpha;
		 return sqrt(d1*d1 + d2*d2 + d3*d3 + d4*d4);
	 };

	 template <typename T1,typename T2>
	 static double EuclideanDistance(const Point2D<T1> &p1, const Point2D<T2> &p2)
	 {
		 double d1 = p1.x - p2.x;
		 double d2 = p1.y - p2.y;
		 return sqrt(d1*d1 + d2*d2);
	 };

	 static int CreatDir(const char * const pDir);

	 static void  listFileByDir(string dir, string fileExtension, vector<string>* fileList);


	 //math
	 // Converts the given polar coordinates of a point to cartesian ones
	 template<class T1, class T2>
	 static void polar2Cartesian(T1 r, T1 t, T2 &y, T2 &x)
	 {
		 x = (T2)(r * cos(t));
		 y = (T2)(r * sin(t));
	 }

	 template<class T>
	 static void layeredGradient(T* data, int h, int w, int layer_no, T* layers, T* workspace = 0, int lwork = 0)
	 {
			 int data_size = h * w;
			 assert(layers != NULL);
			 memset(layers, 0, sizeof(T)*data_size*layer_no);

			 bool empty = false;
			 T* work = NULL;
			 if (lwork < 3 * data_size) {
				 work = new T[3 * data_size];
				 empty = true;
			 }

			 // // smooth the data matrix
			 // T* bdata = blur_gaussian_2d<T,T>( data, h, w, 0.5, 5, false);
			 float kernel[5]; gaussian_1d(kernel, 5, 0.5, 0);
			 memcpy(work, data, sizeof(T)*data_size);
			 convolve_sym(work, h, w, kernel, 5);

			 T *dx = work + data_size;
			 T *dy = work + 2 * data_size;
			 gradient(work, h, w, dy, dx);

#if defined(WITH_OPENMP)
#pragma omp parallel for
#endif
			 for (int l = 0; l<layer_no; l++)
			 {
				 float angle = 2 * l*pi() / layer_no;
				 float kos = cos(angle);
				 float zin = sin(angle);

				 T* layer_l = layers + l*data_size;

				 for (int index = 0; index<data_size; index++)
				 {
					 float value = kos * dx[index] + zin * dy[index];
					 if (value > 0) layer_l[index] = value;
					 else            layer_l[index] = 0;
				 }
			 }
			 if (empty) delete[]work;
	 }
};

