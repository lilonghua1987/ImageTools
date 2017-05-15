/*
 * Representation of a plane in 3D.
 * Written by Simon Fuhrmann.
 */
#pragma  once
#ifndef MATH_PLANE_HEADER
#define MATH_PLANE_HEADER

#include "defines.h"
#include "vector.h"

MATH_NAMESPACE_BEGIN

template <class T> class Plane3;
typedef Plane3<float> Plane3f;
typedef Plane3<double> Plane3d;

/**
 * Class that represents a plane in hesse form.
 * This type of plane allows efficient calculation of orthogonal
 * distances. The normal is expected to have unit length.
 */
template <class T>
class Plane3
{
public:
    typedef Vector<T, 3> Vec3T;

public:
    Vec3T n;
    T d;

public:
    /** Creates an uninitialized plane. */
    Plane3 (void);

    /** Creates a plane with normal n and distance d from the origin. */
    Plane3 (Vec3T const& n, T const& d);

    /** Creates a plane containing p with normal n. */
    Plane3 (Vec3T const& n, Vec3T const& p);

    /** Creates the plane from three points. */
    Plane3 (Vec3T const& p1, Vec3T const& p2, Vec3T const& p3);

    /** Returns the signed distance from a point to the plane. */
    T point_dist (Vec3T const& p) const;

    /** Flips the orientation of the plane. */
    Plane3<T> invert (void) const;

	// Determine whether point P in triangle ABC
	static bool PointinTriangle(Vec3T const& A, Vec3T const& B, Vec3T const& C, Vec3T const& P);

	Vec3T projection(Vec3T const& p);
};

/* ---------------------------------------------------------------- */

template <class T>
inline
Plane3<T>::Plane3 (void)
{
}

template <class T>
inline
Plane3<T>::Plane3 (Vec3T const& n, T const& d)
    : n(n), d(d)
{
}

template <class T>
inline
Plane3<T>::Plane3 (Vec3T const& n, Vec3T const& p)
    : n(n), d(p.dot(n))
{
}

template <class T>
inline
Plane3<T>::Plane3 (Vec3T const& p1, Vec3T const& p2, Vec3T const& p3)
{
    n = (p2 - p1).cross(p3 - p1).normalize();
    d = p1.dot(n);
}

template <class T>
inline T
Plane3<T>::point_dist (Vec3T const& p) const
{
    return p.dot(n) - d;
}

template <class T>
inline Plane3<T>
Plane3<T>::invert (void) const
{
    return Plane3<T>(-n, -d);
}

template <class T>
inline bool
Plane3<T>::PointinTriangle(Vec3T const& A, Vec3T const& B, Vec3T const& C, Vec3T const& P)
{
	Vec3T v0 = C - A ;
	Vec3T v1 = B - A ;
	Vec3T v2 = P - A ;

	T a = v0.dot(v0) ;
	T c = v0.dot(v1) ;
	T d = v0.dot(v2) ;
	T b = v1.dot(v1) ;
	T e = v1.dot(v2) ;

	T f = a * b - b * b;
	T x = b * d - c * e;
	T y = a * e - c * d;

	if (x < 0 || y < 0) // if u out of range, return directly
	{
		return false ;
	}

	return x + y <= f ;
}


template<class T>
inline Vector<T,3> Plane3<T>::projection(Vec3T const& p)
{
	T tM = n.dot(n);
	assert(tM != 0);
	T tZ = n.dot(p) - d;
	return p - n*tZ/tM;
}

MATH_NAMESPACE_END

#endif /* MATH_PLANE_HEADER */
