#include "Fouier.h"
#include<iostream>
using namespace std;
#include<cmath>
#include<ctime>
#define DATALEN 32
#define KEYVALUE 10000 //生成随机浮点数的值，保证分子和分母在这个值之内
#include "tools.h"

Fouier::Fouier(void)
{
   /* data=new complex<double>[DATALEN];
    srand(unsigned int(time(0)));
    cout<<"源数据："<<endl;
    for(int i=0;i<DATALEN;i++)
    {
        data[i]=(rand()%(KEYVALUE))/(double)(rand()%(KEYVALUE)+1);
        if(i%5==0&&i!=0)
            cout<<endl;
        cout<<data[i]<<" ";
    }
    cout<<endl;*/
}


Fouier::~Fouier(void)
{
    delete [] data;
}
complex<double> Fouier:: W(int k,int n)//欧拉公式
{
    double alpha=-2*PI*k/n;
    return complex<double>(cos(alpha),sin(alpha));
}
void Fouier::fft(int start,int step,int len)
{
    if(len==1)//一个元素
    {
        //一个元素不需要变换
        return ;
    }
    fft(start,step*2,len/2);//X1(k)
    fft(start+step,step*2,len/2);//X2(k)
    complex<double> X1,X2;
    for(int i=0;i<len/2;i++)
    {
        X1=data[start+step*i*2];
        X2=data[start+step*(i*2+1)];
        //计算X(k):k=0~N/2-1
        data[start+step*i]=X1+W(i,len)*X2;
        //计算X(k):k=N/2~N-1
        data[start+step*(i+len/2)]=X1-W(i,len)*X2;
    }
}
void Fouier::fft()
{
    fft(0,1,DATALEN);
    cout<<"变换后数据："<<endl;
    for(int i=0;i<DATALEN;i++)
    {
        if(i%5==0&&i!=0)
            cout<<endl;
        cout<<data[i]<<" ";
    }
}


/*************************************************************************
 *
 * 函数名称：
 *   FFT()
 *
 * 参数:
 *   complex<double> * TD	- 指向时域数组的指针
 *   complex<double> * FD	- 指向频域数组的指针
 *   r						－2的幂数，即迭代次数
 *
 * 返回值:
 *   无。
 *
 * 说明:
 *   该函数用来实现快速付立叶变换。
 *
 ************************************************************************/
void Fouier::fft(complex<double> * TD, complex<double> * FD, int r)
{
	// 付立叶变换点数
	long	count;
	
	// 循环变量
	int		i,j,k;
	
	// 中间变量
	int		bfsize,p;
	
	// 角度
	double	angle;
	
	complex<double> *W,*X1,*X2,*X;
	
	// 计算付立叶变换点数
	count = 1 << r;
	
	// 分配运算所需存储器
	W  = new complex<double>[count / 2];
	X1 = new complex<double>[count];
	X2 = new complex<double>[count];
	
	// 计算加权系数
	for(i = 0; i < count / 2; i++)
	{
		angle = -i * PI * 2 / count;
		W[i] = complex<double> (cos(angle), sin(angle));
	}
	
	// 将时域点写入X1
	memcpy(X1, TD, sizeof(complex<double>) * count);
	
	// 采用蝶形算法进行快速付立叶变换
	for(k = 0; k < r; k++)
	{
		for(j = 0; j < 1 << k; j++)
		{
			bfsize = 1 << (r-k);
			for(i = 0; i < bfsize / 2; i++)
			{
				p = j * bfsize;
				X2[i + p] = X1[i + p] + X1[i + p + bfsize / 2];
				X2[i + p + bfsize / 2] = (X1[i + p] - X1[i + p + bfsize / 2]) * W[i * (1<<k)];
			}
		}
		X  = X1;
		X1 = X2;
		X2 = X;
	}
	
	// 重新排序
	for(j = 0; j < count; j++)
	{
		p = 0;
		for(i = 0; i < r; i++)
		{
			if (j&(1<<i))
			{
				p+=1<<(r-i-1);
			}
		}
		FD[j]=X1[p];
	}
	
	// 释放内存
	delete W;
	delete X1;
	delete X2;
}


/*************************************************************************
 *
 * 函数名称：
 *   fourier()
 *
 * 参数:
 *   FIBITMAP *image    - 指向源DIB图像指针
 *
 * 返回值:
 *   complex<double> *            
 *
 * 说明:
 *   该函数用来对图像进行付立叶变换。
 *
 ************************************************************************/
complex<double> * Fouier::fourier(const Image& image)
{
	int width = image.width;  //图像宽度
	int height = image.height;  //图像高度
	
	// 中间变量
	double	dTemp;
	
	// 循环变量
	long	i;
	long	j;
	
	// 进行付立叶变换的宽度和高度（2的整数次方）
	long	w;
	long	h;
	
	int		wp;
	int		hp;
	
	// 赋初值
	w = 1;
	h = 1;
	wp = 0;
	hp = 0;
	
	// 计算进行付立叶变换的宽度和高度（2的整数次方）
	while(w * 2 <= width)
	{
		w *= 2;
		wp++;
	}
	
	while(h * 2 <= width)
	{
		h *= 2;
		hp++;
	}
	
	Image fourierImg(w, h, image.imgType);
	// 分配内存
	complex<double> *TD = new complex<double>[w * h];
	complex<double> *FD = new complex<double>[w * h];
	
	// 行
	for(i = 0; i < h; i++)
	{
		// 列
		for(j = 0; j < w; j++)
		{	
			// 给时域赋值
			if(i<height && j<width){
				TD[j + w * i] = complex<double>(image.at<BYTE>(i,j), 0);
			}else{
				TD[j + w * i] = complex<double>(0, 0);
			}
		}
	}
	
	for(i = 0; i < h; i++)
	{
		// 对y方向进行快速付立叶变换
		fft(&TD[w * i], &FD[w * i], wp);
	}
	
	// 保存变换结果
	for(i = 0; i < h; i++)
	{
		for(j = 0; j < w; j++)
		{
			TD[i + h * j] = FD[j + w * i];
		}
	}
	
	for(i = 0; i < w; i++)
	{
		// 对x方向进行快速付立叶变换
		fft(&TD[i * h], &FD[i * h], hp);
	}
	
	// 行
	for(i = 0; i < h; i++)
	{
		// 列
		for(j = 0; j < w; j++)
		{
			// 计算频谱
			dTemp = sqrt(FD[j * h + i].real() * FD[j * h + i].real() + 
				         FD[j * h + i].imag() * FD[j * h + i].imag()) / 100;
			
			// 判断是否超过255
			if (dTemp > 255) dTemp = 255;
			
			// 指向DIB第(i<h/2 ? i+h/2 : i-h/2)行，第(j<w/2 ? j+w/2 : j-w/2)个象素的指针
			// 此处不直接取i和j，是为了将变换后的原点移到中心
			//lpSrc = (unsigned char*)lpDIBBits + lLineBytes * (lHeight - 1 - i) + j;
			
			// 更新源图像
			fourierImg.at<BYTE>((h - 1 - (i<h/2 ? i+h/2 : i-h/2)),(j<w/2 ? j+w/2 : j-w/2)) = (BYTE)(dTemp);
		}
	}
	string filename = tools::fileNameFromTime(PATH,"testfourierImg",".jpg");
	fourierImg.save(filename.c_str());
	
	delete TD;

	return FD;
}