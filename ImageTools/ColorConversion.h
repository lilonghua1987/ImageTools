#pragma once

#include "Image.h"
#include "Matrix.h"

using namespace std;

class ColorConversion
{
public:
	ColorConversion(void);
	~ColorConversion(void);

private:
	static const double epsilon;	//actual CIE standard
	static const double kappa;		//actual CIE standard
	static double F(double input);
	static int RoundC(const double in_x)
	{
		if (in_x < 0)
			return (int)(in_x - 0.5);
		else
			return (int)(in_x + 0.5);
	}

	static void RGBtoXYZ(double R, double G, double B, double &X, double &Y, double &Z);
	static void XYZtoLab(double X, double Y, double Z, double &L, double &a, double &b);
	static void RGBtoLab(double R, double G, double B, double &L, double &a, double &b);

	static void XYZtoLUV(double X, double Y, double Z, double &L, double &U, double &V);
	static void RGBtoLUV(double R, double G, double B, double &L, double &U, double &V);

	static double inverseF(double input, double factor);
	static void LabtoXYZ(double L, double a, double b, double &X, double &Y, double &Z);
	static void XYZtoRGB(double X, double Y, double Z, double &R, double &G, double &B);
	static void LabtoRGB(double L, double a, double b, double &R, double &G, double &B);

	static void LUVtoXYZ(double L, double U, double V, double &X, double &Y, double &Z);
	static void LUVtoRGB(double L, double U, double V, double &R, double &G, double &B);

	static double EuclideanDistance(double a1, double a2, double a3, double b1, double b2, double b3);

public:
	static Matrix<Pixel<float> > ColorSpaceConversionFromRGB2Lab(const Image &pRGB);
	static Matrix<float> RGBToLab(const Image& img);
	static void RGBToLab(const Image& img, Matrix<float>& dest);
	static void LabtoRGB(const Matrix<float> &imgLab,Image &dest);

	static Matrix<float> RGBToLUV(const Image& img);
	static void RGBtoLUV(const Image& img, Matrix<float>& dest);
	static void LUVtoRGB(const Matrix<float> &luvMat,Image &dest);
};


inline double ColorConversion::F(double input)
{
	if(input > epsilon)
		return (pow(input, 1.0/3.0));
	else
		return (kappa*input + 16.0)/116.0;
}

inline double ColorConversion::EuclideanDistance(double a1, double a2, double a3, double b1, double b2, double b3)
{
	double d1 = a1 - b1;
	double d2 = a2 - b2;
	double d3 = a3 - b3;
	return sqrt(d1*d1 + d2*d2 + d3*d3);
}


// RGB -> XYZ
inline void ColorConversion::RGBtoXYZ(double R, double G, double B, double &X, double &Y, double &Z)
{
	X = R*0.4124564 + G*0.3575761 + B*0.1804375;
	Y = R*0.2126729 + G*0.7151522 + B*0.0721750;
	Z = R*0.0193339 + G*0.1191920 + B*0.9503041;
}


// XYZ -> CIELab
inline void ColorConversion::XYZtoLab(double X, double Y, double Z, double &L, double &a, double &b)
{
	const double Xo=244.66128;
	const double Yo=255.0;
	const double Zo=277.63227;
	L=116*F(Y/Yo)-16;
	a=500*(F(X/Xo)-F(Y/Yo));
	b=200*(F(Y/Yo)-F(Z/Zo));
}



// XYZ -> CIELUV
inline void ColorConversion::XYZtoLUV(double X, double Y, double Z, double &L, double &U, double &V)
{
	const double Yn = 1.00000;
	const double Un_prime = 0.19784977571475;
	const double Vn_prime = 0.46834507665248;
	const double Lt	= 0.008856*255.0;
	
	double L0 = Y / Yn;

	if(L0 > Lt)
		L = (116.0 * (pow(L0, 1.0/3.0)) - 16.0);
	else
		L = (903.3 * L0);

	double constant	= X + 15 * Y + 3 * Z;
	double	u_prime, v_prime;
	if(constant != 0)
	{
		u_prime	= (4 * X) / constant; 
		v_prime = (9 * Y) / constant;
	}else
	{
		u_prime	= 4.0;
		v_prime	= 9.0/15.0;
	}

	/**compute u* and v* */
	U = (13 * L * (u_prime - Un_prime));
	V = (13 * L * (v_prime - Vn_prime));
}



// RGB -> CIELab
inline void ColorConversion::RGBtoLab(double R, double G, double B, double &L, double &a, double &b)
{
	double X, Y, Z;
	RGBtoXYZ(R, G, B, X, Y, Z);
	XYZtoLab(X, Y, Z, L, a, b);
}


// RGB -> CIELUV
inline void ColorConversion::RGBtoLUV(double R, double G, double B, double &L, double &U, double &V)
{
	double X, Y, Z;
	RGBtoXYZ(R, G, B, X, Y, Z);
	XYZtoLUV(X, Y, Z, L, U, V);
}



inline double ColorConversion::inverseF(double input, double factor)
{
	const double thresh = 6.0/29;
	if(input>thresh)
		return (pow(input, 3) * factor);
	else
		return (input - 16.0/116) * 3 * pow(thresh, 2) * factor;
}



// CIELab -> XYZ
inline void ColorConversion::LabtoXYZ(double L, double a, double b, double &X, double &Y, double &Z)
{
	const double Xo=244.66128;
	const double Yo=255.0;
	const double Zo=277.63227;

	double fY = (L + 16)/166;
	double fX = fY + a/500;
	double fZ = fY - b/200;

	X = inverseF(fX,Xo)*2.9;
	Y = inverseF(fY,Yo)*2.9;
	Z = inverseF(fZ,Zo)*2.9;
}



// CIELUV -> XYZ
inline void ColorConversion::LUVtoXYZ(double L, double U, double V, double &X, double &Y, double &Z)
{
	const double Yn	= 1.00000;
	const double Un_prime = 0.19784977571475;
	const double Vn_prime = 0.46834507665248;

	if(L < 8.0)
		Y = Yn * L / 903.3;
	else
	{
		Y = (L + 16.0) / 116.0;
		Y *= Yn * Y * Y;
	}

	double u_prime	= U / (13 * L) + Un_prime;
	double v_prime	= V / (13 * L) + Vn_prime;

	X = 9 * u_prime * Y / (4 * v_prime);
	Z = (12 - 3 * u_prime - 20 * v_prime) * Y / (4 * v_prime);
}



// XYZ -> RGB
inline void ColorConversion::XYZtoRGB(double X, double Y, double Z, double &R, double &G, double &B)
{

	R = RoundC(3.2389*X  - 1.5312*Y - 0.5294*Z);
	G = RoundC(-0.9688*X + 1.8742*Y + 0.0508*Z);
	B = RoundC(0.0556*X - 0.2039*Y + 1.0568*Z);

	if(R < 0)	R = 0; if(R > 255)	R = 255;
	if(G < 0)	G = 0; if(G > 255)	G = 255;
	if(B < 0)	B = 0; if(B > 255)	B = 255;
}



// CIELab -> RGB
inline void ColorConversion::LabtoRGB(double L, double a, double b, double &R, double &G, double &B)
{
	double X, Y, Z;
	LabtoXYZ(L, a, b, X, Y, Z);
	XYZtoRGB(X, Y, Z, R, G, B);
}


// CIELUV -> RGB
inline void ColorConversion::LUVtoRGB(double L, double U, double V, double &R, double &G, double &B)
{
	if(L < 0.1)
	{
		R = 0;
		G = 0;
		B = 0;
		return;
	}
	double X, Y, Z;
	LUVtoXYZ(L, U, V, X, Y, Z);
	XYZtoRGB(X, Y, Z, R, G, B);
}