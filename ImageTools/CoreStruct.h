#pragma once

#ifndef CORESCTRUCT_H

#ifndef uchar
typedef unsigned char uchar;
#endif // !uchar

#ifndef ulong
typedef unsigned long ulong;
#endif

#ifndef uint
typedef unsigned int uint;
#endif

#ifdef _DEBUG
#  define RUNTIMECHECK true
#else
#  define RUNTIMECHECK false
#endif // DEBUG
#define TYPENAME(t1,t2) decltype(t1()*t2())

#include <iostream>
#include <assert.h>

template <class _BTy>
class Pixel
{
public:
	_BTy red;
	_BTy green;
	_BTy blue;
	_BTy alpha;

public:

	Pixel()
		:red(0)
		,green(0)
		,blue(0)
		,alpha(0)
	{
	};

	Pixel(_BTy red, _BTy green, _BTy blue , _BTy alpha = 0)
		:red(red),
		green(green)
		,blue(blue)
		,alpha(alpha)
	{
		static_assert(std::_Is_numeric<_BTy>::value, "_BTy should be numeric");
	};

	template <typename T>
	Pixel(const Pixel<T> &p)
	{
		this->alpha = (_BTy)p.alpha;
		this->blue = (_BTy)p.blue;
		this->green = (_BTy)p.green;
		this->red = (_BTy)p.red;
	};

	template <typename T>
	auto operator *(const T &rate)const ->Pixel<TYPENAME(_BTy,T)>
	{
		static_assert(std::_Is_numeric<T>::value, "T should be numeric");
		return Pixel<TYPENAME(_BTy,T)>(red * rate, green * rate, blue * rate, alpha * rate);
	};

	template <typename T>
	auto operator *(const Pixel<T> &pixel)const ->Pixel<TYPENAME(_BTy,T)>
	{
		return Pixel<TYPENAME(_BTy,T)>(red * pixel.red, green * pixel.green, blue * pixel.blue, alpha * pixel.alpha);
	};

	template <typename T>
	Pixel& operator = (const Pixel<T> &p)
	{
		this->alpha = (_BTy)p.alpha;
		this->blue = (_BTy)p.blue;
		this->green = (_BTy)p.green;
		this->red = (_BTy)p.red;
		return *this;
	};

	template <typename T>
	Pixel& operator += (const Pixel<T> &p)
	{
		this->alpha= static_cast<_BTy>(this->alpha +p.alpha);
		this->blue = static_cast<_BTy>(this->blue + p.blue);
		this->green = static_cast<_BTy>(this->green + p.green);
		this->red = static_cast<_BTy>(this->red + p.red);
		return *this;
	};

	template <typename T>
	auto operator/(const T &t)const ->Pixel<TYPENAME(_BTy,T)>
	{
		static_assert(std::_Is_numeric<T>::value, "T should be numeric");
		assert(t != 0);
		return Pixel<TYPENAME(_BTy,T)>(red / t, green / t, blue / t, alpha / t);
	};

	Pixel operator-(const Pixel& p)const
	{
		return Pixel(red - p.red, green - p.green, blue - p.blue, alpha - p.alpha);
	}

	template<typename T>
	friend auto operator -  (const T& t, const Pixel& p) ->Pixel<TYPENAME(_BTy,T)>
	{
		static_assert(std::_Is_numeric<T>::value, "T should be numeric");
		return Pixel<TYPENAME(_BTy,T)>(t-p.red,t-p.green,t-p.blue,t-p.alpha);
	}

	template<class T>
	auto operator + (const Pixel<T>& p)const ->Pixel<TYPENAME(_BTy,T)>
	{
		static_assert(std::_Is_numeric<T>::value, "T should be numeric");
		return Pixel<TYPENAME(_BTy,T)>(red+p.red,green+p.green,blue+p.blue,alpha+p.alpha);
	}

	template<typename T>
	friend auto operator +  (const T& t, const Pixel& p) ->Pixel<TYPENAME(_BTy,T)>
	{
		static_assert(std::_Is_numeric<T>::value, "T should be numeric");
		return Pixel<TYPENAME(_BTy,T)>(t+p.red,t+p.green,t+p.blue,t+p.alpha);
	}

	bool operator == (const Pixel &p1)const;

	bool operator < (const _BTy& p)const
	{
		static_assert(std::_Is_numeric<_BTy>::value, "_BTy should be numeric");
		return ((this->red < p) || (this->green < p) || (this->blue < p) || (this->alpha < p));
	}

	bool operator > (const _BTy& p)const
	{
		static_assert(std::_Is_numeric<_BTy>::value, "_BTy should be numeric");
		return ((this->red > p) || (this->green > p) || (this->blue > p) || (this->alpha > p));
	}

	Pixel& randColor();

	double distanceL1(const Pixel& p)const
	{
		return abs(red - p.red) + abs(green - p.green) + abs(blue - p.blue) + abs(alpha - p.alpha);
	}

	double distanceL2(const Pixel& p)const
	{
		return (red - p.red) * (red - p.red) + (green - p.green) * (green - p.green) + (blue - p.blue) * (blue - p.blue) + (alpha - p.alpha) * (alpha - p.alpha);
	}
};

template<class _BTy>
inline bool Pixel<_BTy>::operator==(const Pixel<_BTy> &p1)const
{
	static_assert(std::_Is_numeric<_BTy>::value, "_BTy should be numeric");
	return ((p1.red == red) && (p1.green == green) && (p1.blue == blue) && (p1.alpha == alpha));
};

template <class _BTy>
inline Pixel<_BTy>& Pixel<_BTy>::randColor()
{
	static_assert(std::_Is_numeric<_BTy>::value, "_BTy should be numeric");
	this->red = rand()/256;
	this->green = rand()/256;
	this->blue = rand()/256;
	this->alpha = rand()/256;
	return *this;
}



//
template <typename T>
struct Point2D
{
	T x; //表示横坐标，在图片中也可以表示宽度
	T y; //表示纵坐标，在图片中也可以表示高度

	Point2D()
		:x(0)
		,y(0)
	{};

	Point2D(T x, T y)
		:x(x)
		,y(y)
	{};

	template<typename T2>
	Point2D(const Point2D<T2>& point)
		:x(static_cast<T>(point.x))
		,y(static_cast<T>(point.y))
	{}

	template<typename T2>
	auto operator-(const Point2D<T2>& point)const ->Point2D<TYPENAME(T,T2)>
	{
		return Point2D<TYPENAME(T,T2)>(x-point.x,y-point.y);
	}

	Point2D operator-()const
	{
		return Point2D(-x,-y);
	}

	template<typename T2>
	auto operator+(const Point2D<T2>& point)const ->Point2D<TYPENAME(T,T2)>
	{
		return Point2D<TYPENAME(T,T2)>(x+point.x,y+point.y);
	}

	template<typename T2>
	auto operator / (const T2& t)const ->Point2D<TYPENAME(T,T2)>
	{
		static_assert(std::_Is_numeric<T>::value, "T should be numeric");
		static_assert(std::_Is_numeric<T2>::value, "t should be numeric");
		assert(t != 0);
		return Point2D<TYPENAME(T,T2)>(x/t,y/t);
	}

	template<typename T2>
	Point2D& operator=(const Point2D<T2>& point)
	{
		this->x = static_cast<T>(point.x);
		this->y = static_cast<T>(point.y);
		return *this;
	}

	template<typename T2>
	bool operator==(const Point2D<T2>& point)const
	{
		return (point.x == x && point.y == y);
	}

	template<typename T2>
	bool operator!=(const Point2D<T2>& point)const
	{
		return (point.x != x || point.y != y);
	}

	bool operator>(const Point2D& point)const
	{
		return (point.x < x && point.y < y);
	}

	bool operator>=(const Point2D& point)const
	{
		return (point.x <= x && point.y <= y);
	}

	bool operator<(const Point2D& point)const
	{
		return (point.x > x && point.y > y);
	}

	bool operator<=(const Point2D& point)const
	{
		return (point.x >= x && point.y >= y);
	}

	template<typename T2>
	friend std::ostream& operator << (std::ostream &outs, const Point2D<T2>& point);

	template<typename T2>
	friend std::istream& operator >> (std::istream &ins, Point2D<T2>& point);
};


template<typename T2>
inline std::ostream& operator << (std::ostream& outs, const Point2D<T2>& point)
{
	outs << "[" << point.x << "," << point.y << "]";
	return outs;
};

template<typename T2>
inline std::istream& operator >> (std::istream& ins, Point2D<T2>& point)
{
	std::cout << "Please input the point [x y]" << std::endl;
	ins >> point.x >> point.y;
	return ins;
};


template <typename T>
struct Point3D:public Point2D<T>
{
	T z;

	Point3D()
		:Point2D()
		,z(0)
	{};

	Point3D(T x, T y, T z)
		:Point2D(x,y)
		,z(z)
	{};

	Point3D(const Point3D& point)
		:Point2D(point.x,point.y)
		,z(point.z)
	{}

	Point3D(const Point2D& point, T z)
		:Point2D(point)
		,z(z)
	{}

	template<class T1>
	auto operator-(const Point3D<T1>& point)const ->Point3D<TYPENAME(T,T1)>
	{
		return Point3D<TYPENAME(T,T1)>(x - point.x,y - point.y,z - point.z);
	}

	Point3D operator-()const
	{
		return Point3D(-x,-y,-z);
	}

	template<class T1>
	auto operator+(const Point3D<T1>& point)const ->Point3D<TYPENAME(T,T1)>
	{
		return Point3D<TYPENAME(T,T1)>(x + point.x,y + point.y,z + point.z);
	}

	template<typename T1>
	auto operator / (const T1& t)const ->Point3D<TYPENAME(T,T1)>
	{
		static_assert(std::_Is_numeric<T>::value, "T should be numeric");
		static_assert(std::_Is_numeric<T1>::value, "t should be numeric");
		assert(t != 0);
		return Point3D<TYPENAME(T,T1)>(x/t,y/t,z/t);
	}

	template<typename T1>
	auto operator * (const T1& t)const ->Point3D<TYPENAME(T,T1)>
	{
		static_assert(std::_Is_numeric<T>::value, "T should be numeric");
		static_assert(std::_Is_numeric<T1>::value, "t should be numeric");
		return Point3D<TYPENAME(T,T1)>(x*t,y*t,z*t);
	}

	Point3D& operator=(const Point3D& point)
	{
		this->x = point.x;
		this->y = point.y;
		this->z = point.z;
		return *this;
	}

	bool operator==(const Point3D& point)const
	{
		return (point.x == x && point.y == y && point.z == z);
	}

	bool operator!=(const Point3D& point)const
	{
		return (point.x != x || point.y != y || point.z != z);
	}

	bool operator>(const Point3D& point)const
	{
		return (point.x < x && point.y < y && point.z < z);
	}

	bool operator>=(const Point3D& point)const
	{
		return (point.x <= x && point.y <= y && point.z <= z);
	}

	bool operator<(const Point3D& point)const
	{
		return (point.x > x && point.y > y && point.z > z);
	}

	bool operator<=(const Point3D& point)const
	{
		return (point.x >= x && point.y >= y && point.z >= z);
	}

	virtual double distanceL1(const Point3D& p)const
	{
		return abs(x - p.x) + abs(y - p.y) + abs(z - p.z);
	}

	virtual double distanceL2(const Point3D& p)const
	{
		return (x - p.x) * (x - p.x) + (y - p.y) * (y - p.y) + (z - p.z) * (z - p.z);
	}

	friend std::ostream& operator << (std::ostream& outs, const Point3D& point)
	{
		outs << "[" << point.x << "," << point.y << "," << point.z << "]" << std::endl;
		return outs;
	}

	friend std::istream& operator >> (std::istream& ins, Point3D& point)
	{
		ins >> point.x >> point.y >> point.z;
		return ins;
	}
};

#endif // !CORESCTRUCT_H
