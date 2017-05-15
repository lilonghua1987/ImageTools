#include "Reconstruction.h"

Reconstruction::Reconstruction(void)
{
}


Reconstruction::~Reconstruction(void)
{
}


void Reconstruction::shading(const Image& img, int iter, float Ps, float Qs)
{
	int width = img.width,height = img.height;
	float p,q,pq,PQs,fZ,dfZ,Eij,Wn=0.0001f*0.0001f,Y,K;

	Matrix<float> Zn(height,width),Zn1(height,width),Si(height,width),Si1(height,width);

	Matrix<float> result(height,width);

	for(int i=0;i<height;i++)
	{
       for(int j=0;j<width;j++)
	   {
           Si1.at(i,j) = 0.01f; 
	   }
	}

	/************************************************************************/
	for(int I=1;I<=iter;I++){
		for(int i=0;i<height;i++)
		{
			for(int j=0;j<width;j++){ /* calculate -f(Zij) & df(Zij) */
				if(j-1 < 0 || i-1 < 0) /* take care boundary */
				{
					p = q = 0.0;
				}else {
					p = Zn1.at(i,j) - Zn1.at(i,(j-1));
					q = Zn1.at(i,j) - Zn1.at((i-1),j);
				}
				pq = 1.0f + p*p + q*q;
				PQs = 1.0f + Ps*Ps + Qs*Qs;
				Eij = img.at<BYTE>(i,j)/255.0f;
				fZ = -1.0f*(Eij - tools::Max(0.0f,(1+p*Ps+q*Qs)/(sqrt(pq)*sqrt(PQs))));
				dfZ = static_cast<float>(-1.0f*((Ps+Qs)/(sqrt(pq)*sqrt(PQs))-(p+q)*(1.0+p*Ps+q*Qs)/(sqrt(pq*pq*pq)*sqrt(PQs)))) ;
				Y = fZ + dfZ*Zn1.at(i,j);
				K = Si1.at(i,j)*dfZ/(Wn+dfZ*Si1.at(i,j)*dfZ);
				Si.at(i,j) = (1.0f - K*dfZ)*Si1.at(i,j); 
				Zn.at(i,j) = Zn1.at(i,j) + K*(Y-dfZ*Zn1.at(i,j));
			}
		}

		for(int i=0;i<height;i++)
		{
			for(int j=0;j<width;j++)
			{
				result.at(i,j) = Zn.at(i,j);
				Zn1.at(i,j) = Zn.at(i,j);
				Si1.at(i,j) = Si.at(i,j);
			}
		}
	}

	/*Log log("shadImg.txt");
	log.systemLog(result);*/

	ImageTools::SaveImage(img,"temp/shadimg.jpg");
}