#pragma once
#include"complex"
#include "ImageProcess.h"

class Fouier
{
	complex<double> * data;
    void fft(int start,int step,int len);
    complex<double> W(int k,int n);//e^(-i*2*pi*k/n)
public:
    Fouier(void);
    ~Fouier(void);
    
    void fft();

	void fft(complex<double> * TD, complex<double> * FD, int r);

	complex<double> * fourier(const Image& image);
};