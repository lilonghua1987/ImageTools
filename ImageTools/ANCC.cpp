#include "StereoMatch.h"


void StereoMatch::ANCC(float theta, int w)
{
	Matrix<float> logL = ImageExpandTools::logImage(imgL);
	Matrix<float> logR = ImageExpandTools::logImage(imgR);

	for (ulong r = 0; r < logL.row; r++)
	{
		for (ulong c = 0; c < logL.column; c++)
		{
			float meanL = 0.f;
			float meanR = 0.f;
			for (ulong i = 0; i < logL.channel; i++)
			{
				meanL += logL.at(r,c,i);
				meanR += logR.at(r,c,i);
			}
			meanL /= 3.f;
			meanR /= 3.f;

			for (ulong i = 0; i < logL.channel; i++)
			{
				logL.at(r,c,i) = (meanL - logL.at(r,c,i))*100;
				logR.at(r,c,i) = (meanR - logR.at(r,c,i))*100;
			}
		}
	}

	Matrix<float> labL = ColorConversion::RGBToLab(imgL);
	Matrix<float> labR = ColorConversion::RGBToLab(imgR);

	for (int r = 0; r < (int)imgL.height; r++)
	{
		for (int c = 0; c < (int)imgL.width; c++)
		{

		}
	}
}