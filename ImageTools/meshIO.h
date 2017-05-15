#pragma once

#include "Image.h"
#include "Matrix.h"

template<class T>
class MeshIO
{
public:
	MeshIO();
	~MeshIO();

	static bool exportPLY(Image& d, const Image& i, float sf, const std::string& path)
	{
		assert(d.channel == 1);
		assert(d.width != 0 && d.height != 0);
		if (d.width == 0 || d.height == 0)
		{
			return false;
		}
		assert(i.width != 0 && i.height != 0);
		if (i.width == 0 || i.height == 0)
		{
			return false;
		}
		assert(d.height == i.height && d.width == i.width);
		if (d.height != i.height || d.width != i.width)
		{
			return false;
		}
		assert(sf > 0);
		if (sf <= 0)
		{
			return false;
		}
		assert(!path.empty());

		unsigned int vertex = 0;
		for (ulong y = 0; y < d.height; y++)
		{
			for (ulong x = 0; x < d.width; x++)
			{
				if (d.at<uchar>(y, x)) vertex++;
			}
		}

		ofstream out;
		out.open(path);
		out << "ply" << endl;
		out << "format ascii 1.0" << endl;
		out << "comment created by DisparityMap2Mesh" << endl;
		out << "element vertex " << vertex << endl;
		out << "property float x" << endl;
		out << "property float y" << endl;
		out << "property float z" << endl;
		out << "property uchar red" << endl;
		out << "property uchar green" << endl;
		out << "property uchar blue" << endl;
		out << "end_header" << endl;

		const float maxV = (float)d.maxVal<uchar>();
		const float minV = (float)d.minVal<uchar>();
		const float d0 = (float)d.width;
		const float deltaZ = (float)d.width / 4.0f; // /10.0f;
		const float a = deltaZ / (1.0f / d0 - 1.0f / (maxV - minV + d0));
		const float z0 = deltaZ - a / d0;

		for (ulong y = 0; y < d.height; y++)
		{
			for (ulong x = 0; x < d.width; x++)
			{
				Pixel<BYTE> point = i.getPixel<BYTE>(y, x);
				if (d.at<BYTE>(y, x))
				{
					float z = a / (d.at<BYTE>(y, x) - minV + d0) + z0;
					out << x << "\t" << (float)(d.height - 1 - y) << "\t" << z << "\t" << (int)point.red << "\t" << (int)point.green << "\t" << (int)point.blue << std::endl;
				}
			}
		}
		out.close();
		return true;
	};

	static bool exportPLY(Matrix<T>& d, const Image& img, float sf, const std::string& path)
	{
		assert(d.channel == 1);
		assert(d.column != 0 && d.row != 0);
		assert(img.width != 0 && img.height != 0);
		assert(d.row == img.height && d.column == img.width);
		assert(!path.empty());

		unsigned int vertex = 0;
		for (ulong y = 0; y < d.row; y++)
		{
			for (ulong x = 0; x < d.column; x++)
			{
				if (d.at(y, x)) vertex++;
			}
		}

		ofstream out;
		out.open(path);
		out << "ply" << endl;
		out << "format ascii 1.0" << endl;
		out << "comment created by DisparityMap2Mesh" << endl;
		out << "element vertex " << vertex << endl;
		out << "property float x" << endl;
		out << "property float y" << endl;
		out << "property float z" << endl;
		out << "property uchar red" << endl;
		out << "property uchar green" << endl;
		out << "property uchar blue" << endl;
#ifdef TRIANGULATION
		out << "element face " << triangles.size() << endl;
		out << "property list uchar int vertex_index" << endl;
#endif
		out << "end_header" << endl;

		const float maxV = (float)d.maxValue();
		const float minV = (float)d.minValue();
		const float d0 = (float)d.column;
		const float deltaZ = (float)tools::Max(d.column, d.row); // /10.0f;
		const float a = deltaZ / (1.0f / d0 - 1.0f / (maxV - minV + d0));
		const float z0 = deltaZ - a / d0;

		for (ulong y = 0; y < d.row; y++)
		{
			for (ulong x = 0; x < d.column; x++)
			{
				Pixel<BYTE> point = img.getPixel<BYTE>(y, x);
				if (d.at(y, x) != 0.f)
				{
					float z = a / (d.at(y, x) - minV + d0) + z0;
					out << x << "\t" << (d.row - 1 - y) << "\t" << z << "\t" << (int)point.red << "\t" << (int)point.green << "\t" << (int)point.blue << endl;
				}
			}
		}
		out.close();
		return true;
	};
};