#include "BPStereo.h"


const float BPStereo::discK = 1.7f;
const float BPStereo::dataK = 15.f;
const float BPStereo::lambda = 0.07f;

BPStereo::BPStereo(int disp,int iter,int levels,float sigma)
	:disp(disp)
	,iter(iter)
	,levels(levels)
	,sigma(sigma)
{
}


BPStereo::~BPStereo(void)
{
}

// dt of 1d function
inline void BPStereo::dt(float f[]) {
	for (auto q = 1; q < disp; q++) {
		float prev = f[q-1] + 1.0F;
		if (prev < f[q])
			f[q] = prev;
	}
	for (auto q = disp-2; q >= 0; q--) {
		float prev = f[q+1] + 1.0F;
		if (prev < f[q])
			f[q] = prev;
	}
}

void BPStereo::msg(float s1[], float s2[], float s3[], float s4[], float dst[])
{
	float val;

	// aggregate and find min
	float minimum = (float)infVal;
	for (auto value = 0; value < disp; value++) {
		dst[value] = s1[value] + s2[value] + s3[value] + s4[value];
		if (dst[value] < minimum)
			minimum = dst[value];
	}

	// dt
	dt(dst);

	// truncate 
	minimum += discK;
	for (auto value = 0; value < disp; value++)
		if (minimum < dst[value])
			dst[value] = minimum;

	// normalize
	val = 0;
	for (auto value = 0; value < disp; value++) 
		val += dst[value];

	val /= disp;
	for (auto value = 0; value < disp; value++) 
		dst[value] -= val;
}

Matrix<float> BPStereo::dataCosts(const Image& imL,const Image& imR)
{
	int width = imL.width;
	int height = imL.height;
	Matrix<float> costs(height, width, disp);

	Matrix<float> sm1 = ImageExpandTools::Image2FloatMat(imL), sm2 = ImageExpandTools::Image2FloatMat(imR);
	if (sigma >= 0.1) {
		sm1 = ImageProcess::GausscianSeparateBlur(sm1, sigma);
		sm2 = ImageProcess::GausscianSeparateBlur(sm2, sigma);
	}

	for (int y = 0; y < height; y++) 
	{
		for (int x = disp-1; x < width; x++)
		{
			for (auto value = 0; value < disp; value++)
			{
				float val = abs(sm1.at(y, x)-sm2.at(y, x-value));	
				costs.at(y, x, value) = lambda * tools::Min(val, dataK);
			}
		}
	}

	return costs;
}

Matrix<float> BPStereo::dataCosts(const Image& img)
{
	int width = img.width;
	int height = img.height;
	Matrix<float> costs(height, width, disp);

	for (int y = 0; y < height; y++) {
		for (auto x = disp-1; x < width; x++) {
			for (auto value = 0; value < disp; value++) {
				float val = powf((img.getPixel<float>(y, x).red-value),2.0);
				costs.at(y, x, value) = lambda * tools::Min(val, dataK);
			}
		}
	}

	return costs;
}

void BPStereo::bp(Matrix<float>& u, Matrix<float>& d,Matrix<float>& l, Matrix<float>& r,Matrix<float>& data,int iter)
{
	int width = data.column;
	int height = data.row;

	for (auto t = 0; t < iter; t++) {
		std::cout << "iter " << t << "\n";

		for (int y = 1; y < height-1; y++) {
			for (int x = ((y+t) % 2) + 1; x < width-1; x+=2) {

				msg(u.ptr(y+1, x),l.ptr(y, x+1),r.ptr(y, x-1),data.ptr(y, x), u.ptr(y, x));

				msg(d.ptr(y-1, x),l.ptr(y, x+1),r.ptr(y, x-1),data.ptr(y, x), d.ptr(y, x));

				msg(u.ptr(y+1, x),d.ptr(y-1, x),r.ptr(y, x-1),data.ptr(y, x), r.ptr(y, x));

				msg(u.ptr(y+1, x),d.ptr(y-1, x),l.ptr(y, x+1),data.ptr(y, x), l.ptr(y, x));

			}
		}
	}
}

Image BPStereo::generate(Matrix<float>& u, Matrix<float>& d, Matrix<float>& l, Matrix<float>& r, Matrix<float>& data)
{
	int width = data.column;
	int height = data.row;
	Image out(width, height);

	for (int y = 1; y < height-1; y++) {
		for (int x = 1; x < width-1; x++) {
			// keep track of best value for current pixel
			int best = 0;
			float best_val = (float)infVal;
			for (auto value = 0; value < disp; value++) 
			{
				float val = u.at(y+1, x, value) + d.at(y-1, x, value) + l.at(y, x+1, value) + r.at(y, x-1, value) + data.at(y, x, value);
				if (val < best_val) 
				{
					best_val = val;
					best = value;
				}
			}
			out.at<uchar>(y,x) = best;
		}
	}

	ImageTools::NormalImage(out,0,255);
	return out;
}


Image BPStereo::bpMatch(const Image& imL,const Image& imR)
{
	std::vector<Matrix<float>> u(levels);
	std::vector<Matrix<float>> d(levels);
	std::vector<Matrix<float>> l(levels);
	std::vector<Matrix<float>> r(levels);
	std::vector<Matrix<float>> data(levels);

	// data costs
	data[0] = dataCosts(imL,imR);

	// data pyramid
	for (int i = 1; i < levels; i++) {
		int old_width = data[i-1].column;
		int old_height = data[i-1].row;
		int new_width = (int)ceil(old_width/2.0);
		int new_height = (int)ceil(old_height/2.0);

		assert(new_width >= 1);
		assert(new_height >= 1);

		data[i] = Matrix<float>(new_height, new_width, disp);
		for (int y = 0; y < old_height; y++) 
		{
			for (int x = 0; x < old_width; x++) 
			{
				for (auto value = 0; value < disp; value++) 
				{
					data[i].at(y/2, x/2, value) += data[i-1].at(y, x, value);
				}
			}
		}
	}

	// run bp from coarse to fine
	for (int i = levels-1; i >= 0; i--) 
	{
		int width = data[i].column;
		int height = data[i].row;

		// allocate & init memory for messages
		if (i == levels-1) {
			// in the coarsest level messages are initialized to zero
			u[i] =  Matrix<float>(height, width, disp);
			d[i] =  Matrix<float>(height, width, disp);
			l[i] =  Matrix<float>(height, width, disp);
			r[i] =  Matrix<float>(height, width, disp);
		} else {
			// initialize messages from values of previous level
			u[i] =  Matrix<float>(height, width, disp);
			d[i] =  Matrix<float>(height, width, disp);
			l[i] =  Matrix<float>(height, width, disp);
			r[i] =  Matrix<float>(height, width, disp);

			for (int y = 0; y < height; y++) {
				for (int x = 0; x < width; x++) {
					for (auto value = 0; value < disp; value++) {
						u[i].at(y, x, value) = u[i+1].at(y/2, x/2, value);
						d[i].at(y, x, value) = d[i+1].at(y/2, x/2, value);
						l[i].at(y, x, value) = l[i+1].at(y/2, x/2, value);
						r[i].at(y, x, value) = r[i+1].at(y/2, x/2, value);
					}
				}
			}
		} 

		// BP
		bp(u[i], d[i], l[i], r[i], data[i], iter);    
	}

	return generate(u[0], d[0], l[0], r[0], data[0]);
}


Image BPStereo::restore(const Image& img)
{
	std::vector<Matrix<float>> u(levels);
	std::vector<Matrix<float>> d(levels);
	std::vector<Matrix<float>> l(levels);
	std::vector<Matrix<float>> r(levels);
	std::vector<Matrix<float>> data(levels);

	// data costs
	data[0] = dataCosts(img);

	// data pyramid
	for (int i = 1; i < levels; i++) {
		int old_width = data[i-1].column;
		int old_height = data[i-1].row;
		int new_width = (int)ceil(old_width/2.0);
		int new_height = (int)ceil(old_height/2.0);

		assert(new_width >= 1);
		assert(new_height >= 1);

		data[i] = Matrix<float>(new_height, new_width, disp);
		for (int y = 0; y < old_height; y++) 
		{
			for (int x = 0; x < old_width; x++) 
			{
				for (int value = 0; value < disp; value++)
				{
					data[i].at(y/2, x/2, value) += data[i-1].at(y, x, value);
				}
			}
		}
	}

	// run bp from coarse to fine
	for (int i = levels-1; i >= 0; i--) {
		int width = data[i].column;
		int height = data[i].row;

		// allocate & init memory for messages
		if (i == levels-1) {
			// in the coarsest level messages are initialized to zero
			u[i] = Matrix<float>(height, width, disp);
			d[i] = Matrix<float>(height, width, disp);
			l[i] = Matrix<float>(height, width, disp);
			r[i] = Matrix<float>(height, width, disp);
		} else {
			// initialize messages from values of previous level
			u[i] = Matrix<float>(height, width, disp);
			d[i] = Matrix<float>(height, width, disp);
			l[i] = Matrix<float>(height, width, disp);
			r[i] = Matrix<float>(height, width, disp);

			for (int y = 0; y < height; y++) {
				for (int x = 0; x < width; x++) {
					for (auto value = 0; value < disp; value++) {
						u[i].at(y, x, value) = u[i+1].at(y/2, x/2, value);
						d[i].at(y, x, value) = d[i+1].at(y/2, x/2, value);
						l[i].at(y, x, value) = l[i+1].at(y/2, x/2, value);
						r[i].at(y, x, value) = r[i+1].at(y/2, x/2, value);
					}
				}
			}
		} 

		// BP
		bp(u[i], d[i], l[i], r[i], data[i], iter);    
	}

	Image out = generate(u[0], d[0], l[0], r[0], data[0]);

	return out;
}