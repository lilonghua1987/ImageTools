#include "Fouier.h"
#include<iostream>
using namespace std;
#include<cmath>
#include<ctime>
#define DATALEN 32
#define KEYVALUE 10000 //���������������ֵ����֤���Ӻͷ�ĸ�����ֵ֮��
#include "tools.h"

Fouier::Fouier(void)
{
   /* data=new complex<double>[DATALEN];
    srand(unsigned int(time(0)));
    cout<<"Դ���ݣ�"<<endl;
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
complex<double> Fouier:: W(int k,int n)//ŷ����ʽ
{
    double alpha=-2*PI*k/n;
    return complex<double>(cos(alpha),sin(alpha));
}
void Fouier::fft(int start,int step,int len)
{
    if(len==1)//һ��Ԫ��
    {
        //һ��Ԫ�ز���Ҫ�任
        return ;
    }
    fft(start,step*2,len/2);//X1(k)
    fft(start+step,step*2,len/2);//X2(k)
    complex<double> X1,X2;
    for(int i=0;i<len/2;i++)
    {
        X1=data[start+step*i*2];
        X2=data[start+step*(i*2+1)];
        //����X(k):k=0~N/2-1
        data[start+step*i]=X1+W(i,len)*X2;
        //����X(k):k=N/2~N-1
        data[start+step*(i+len/2)]=X1-W(i,len)*X2;
    }
}
void Fouier::fft()
{
    fft(0,1,DATALEN);
    cout<<"�任�����ݣ�"<<endl;
    for(int i=0;i<DATALEN;i++)
    {
        if(i%5==0&&i!=0)
            cout<<endl;
        cout<<data[i]<<" ";
    }
}


/*************************************************************************
 *
 * �������ƣ�
 *   FFT()
 *
 * ����:
 *   complex<double> * TD	- ָ��ʱ�������ָ��
 *   complex<double> * FD	- ָ��Ƶ�������ָ��
 *   r						��2������������������
 *
 * ����ֵ:
 *   �ޡ�
 *
 * ˵��:
 *   �ú�������ʵ�ֿ��ٸ���Ҷ�任��
 *
 ************************************************************************/
void Fouier::fft(complex<double> * TD, complex<double> * FD, int r)
{
	// ����Ҷ�任����
	long	count;
	
	// ѭ������
	int		i,j,k;
	
	// �м����
	int		bfsize,p;
	
	// �Ƕ�
	double	angle;
	
	complex<double> *W,*X1,*X2,*X;
	
	// ���㸶��Ҷ�任����
	count = 1 << r;
	
	// ������������洢��
	W  = new complex<double>[count / 2];
	X1 = new complex<double>[count];
	X2 = new complex<double>[count];
	
	// �����Ȩϵ��
	for(i = 0; i < count / 2; i++)
	{
		angle = -i * PI * 2 / count;
		W[i] = complex<double> (cos(angle), sin(angle));
	}
	
	// ��ʱ���д��X1
	memcpy(X1, TD, sizeof(complex<double>) * count);
	
	// ���õ����㷨���п��ٸ���Ҷ�任
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
	
	// ��������
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
	
	// �ͷ��ڴ�
	delete W;
	delete X1;
	delete X2;
}


/*************************************************************************
 *
 * �������ƣ�
 *   fourier()
 *
 * ����:
 *   FIBITMAP *image    - ָ��ԴDIBͼ��ָ��
 *
 * ����ֵ:
 *   complex<double> *            
 *
 * ˵��:
 *   �ú���������ͼ����и���Ҷ�任��
 *
 ************************************************************************/
complex<double> * Fouier::fourier(const Image& image)
{
	int width = image.width;  //ͼ����
	int height = image.height;  //ͼ��߶�
	
	// �м����
	double	dTemp;
	
	// ѭ������
	long	i;
	long	j;
	
	// ���и���Ҷ�任�Ŀ�Ⱥ͸߶ȣ�2�������η���
	long	w;
	long	h;
	
	int		wp;
	int		hp;
	
	// ����ֵ
	w = 1;
	h = 1;
	wp = 0;
	hp = 0;
	
	// ������и���Ҷ�任�Ŀ�Ⱥ͸߶ȣ�2�������η���
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
	// �����ڴ�
	complex<double> *TD = new complex<double>[w * h];
	complex<double> *FD = new complex<double>[w * h];
	
	// ��
	for(i = 0; i < h; i++)
	{
		// ��
		for(j = 0; j < w; j++)
		{	
			// ��ʱ��ֵ
			if(i<height && j<width){
				TD[j + w * i] = complex<double>(image.at<BYTE>(i,j), 0);
			}else{
				TD[j + w * i] = complex<double>(0, 0);
			}
		}
	}
	
	for(i = 0; i < h; i++)
	{
		// ��y������п��ٸ���Ҷ�任
		fft(&TD[w * i], &FD[w * i], wp);
	}
	
	// ����任���
	for(i = 0; i < h; i++)
	{
		for(j = 0; j < w; j++)
		{
			TD[i + h * j] = FD[j + w * i];
		}
	}
	
	for(i = 0; i < w; i++)
	{
		// ��x������п��ٸ���Ҷ�任
		fft(&TD[i * h], &FD[i * h], hp);
	}
	
	// ��
	for(i = 0; i < h; i++)
	{
		// ��
		for(j = 0; j < w; j++)
		{
			// ����Ƶ��
			dTemp = sqrt(FD[j * h + i].real() * FD[j * h + i].real() + 
				         FD[j * h + i].imag() * FD[j * h + i].imag()) / 100;
			
			// �ж��Ƿ񳬹�255
			if (dTemp > 255) dTemp = 255;
			
			// ָ��DIB��(i<h/2 ? i+h/2 : i-h/2)�У���(j<w/2 ? j+w/2 : j-w/2)�����ص�ָ��
			// �˴���ֱ��ȡi��j����Ϊ�˽��任���ԭ���Ƶ�����
			//lpSrc = (unsigned char*)lpDIBBits + lLineBytes * (lHeight - 1 - i) + j;
			
			// ����Դͼ��
			fourierImg.at<BYTE>((h - 1 - (i<h/2 ? i+h/2 : i-h/2)),(j<w/2 ? j+w/2 : j-w/2)) = (BYTE)(dTemp);
		}
	}
	string filename = tools::fileNameFromTime(PATH,"testfourierImg",".jpg");
	fourierImg.save(filename.c_str());
	
	delete TD;

	return FD;
}