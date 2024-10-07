#ifndef VECTORUTILITIES_H
#define VECTORUTILITIES_H

#include <array>
#include <cmath> // sqrt
#include "Definitions.hpp"

// real_number dotProduct(const Point<DIM, real_number> &v, const Point<DIM, real_number> &w);

// std::array<Point<DIM, real_number>, DIM> dyadicProduct(const Point<DIM, real_number> &v, const Point<DIM, real_number> &w);

// Point<DIM, real_number> matVec(const std::array<Point<DIM, real_number>, DIM> &m, const Point<DIM, real_number> &v);

// real_number getVectorNorm(const Point<DIM, real_number> &v);

// void normalizeVector(Point<DIM, real_number> &v);

// Point<DIM, real_number> getPerpendicularUnit2D(const Point<DIM, real_number> &v);

// std::array<Point<DIM, real_number>, 3> getBoundaryPositions(const Point<DIM, real_number> &r, const Point<DIM, real_number> &normal, real_number dp);

// Point<DIM, real_number> crossProduct(const Point<DIM, real_number> &v, const Point<DIM, real_number> &w);

// void ApplyRotation(Point<DIM, real_number> &x, const real_number theta, const Point<DIM, real_number> centre);

inline __device__ __host__ real_number dotProduct(const Point<DIM, real_number> &v, const Point<DIM, real_number> &w)
{
	real_number result = 0.0f;
	if constexpr (DIM == 2)
		result = v.get(0) * w.get(0) + v.get(1) * w.get(1);
	else if constexpr (DIM == 3)
		result = v.get(0) * w.get(0) + v.get(1) * w.get(1) + v.get(2) * w.get(2);

	return result;
}

inline __device__ __host__ std::array<Point<DIM, real_number>, DIM> dyadicProduct(const Point<DIM, real_number> &v, const Point<DIM, real_number> &w)
{
	std::array<Point<DIM, real_number>, DIM> dyad;
	for (size_t i = 0; i < DIM; i++)
	{
		dyad[i] = v.get(i) * w;
	}
	return dyad;
}

inline __device__ __host__ Point<DIM, real_number> matVec(const std::array<Point<DIM, real_number>, DIM> &m, const Point<DIM, real_number> &v)
{
	Point<DIM, real_number> res;
	for (size_t i = 0; i < DIM; i++)
	{
		res.get(i) = dotProduct(m[i], v);
	}
	return res;
}

inline __device__ __host__ real_number getVectorNorm(const Point<DIM, real_number> &v)
{
	return sqrt(norm2(v));
}

inline __device__ __host__ void normalizeVector(Point<DIM, real_number> &v)
{
	const real_number norm = getVectorNorm(v);
	if (norm > 0.0f)
		v = v / norm;
}

inline __device__ __host__ Point<DIM, real_number> getPerpendicularUnit2D(const Point<DIM, real_number> &v)
{
	Point<DIM, real_number> perp;
	perp.get(0) = -v.get(1);
	perp.get(1) = v.get(0);
	normalizeVector(perp);
	return perp;
}

inline __device__ __host__ std::array<Point<DIM, real_number>, 3> getBoundaryPositions(const Point<DIM, real_number> &r, const Point<DIM, real_number> &normal, real_number dp)
{
	// aplies offset to a vector and returns the vectors pointing at the virtual wall particles
	// r needs to be pointing to the wall particle

	Point<DIM, real_number> offset_1 = -0.5f * normal * dp; // offset from wall to first particle
	Point<DIM, real_number> offset_2 = -1.0f * normal * dp; // offset from first to second particle, and from second to third
	std::array<Point<DIM, real_number>, 3> r_virtual;
	r_virtual[0] = r + offset_1;
	r_virtual[1] = r_virtual[0] + offset_2;
	r_virtual[2] = r_virtual[1] + offset_2;

	return r_virtual;
}

inline __device__ __host__ real_number crossProduct2D(const Point<DIM, real_number> &v, const Point<DIM, real_number> &w)
{
	return v.get(0) * w.get(1) - v.get(1) * w.get(0);
}

inline __device__ __host__ Point<DIM, real_number> crossProduct3D(const Point<DIM, real_number> &v, const Point<DIM, real_number> &w)
{
	Point<DIM, real_number> res;
	res.get(0) = v.get(1) * w.get(2) - v.get(2) * w.get(1);
	res.get(1) = v.get(2) * w.get(0) - v.get(0) * w.get(2);
	res.get(2) = v.get(0) * w.get(1) - v.get(1) * w.get(0);
	return res;
}

inline __device__ __host__ Point<DIM, real_number> ApplyRotation(const Point<DIM, real_number> x, const real_number theta, const Point<DIM, real_number> centre)
{

	Point<DIM, real_number> x_rotated;
	x_rotated.get(0) = cos(theta) * (x.get(0) - centre.get(0)) - sin(theta) * (x.get(1) - centre.get(1)) + centre.get(0);
	x_rotated.get(1) = sin(theta) * (x.get(0) - centre.get(0)) + cos(theta) * (x.get(1) - centre.get(1)) + centre.get(1);

	return x_rotated;
}

// function to know the type of a variable
// used to know types of auto variables
// template <class T>
// std::string type_name()
// {
// 	typedef typename std::remove_reference<T>::type TR;
// 	std::unique_ptr<char, void (*)(void *)> own(
// #ifndef _MSC_VER
// 		abi::__cxa_demangle(typeid(TR).name(), nullptr,
// 							nullptr, nullptr),
// #else
// 		nullptr,
// #endif
// 		std::free);
// 	std::string r = own != nullptr ? own.get() : typeid(TR).name();
// 	if (std::is_const<TR>::value)
// 		r += " const";
// 	if (std::is_volatile<TR>::value)
// 		r += " volatile";
// 	if (std::is_lvalue_reference<T>::value)
// 		r += "&";
// 	else if (std::is_rvalue_reference<T>::value)
// 		r += "&&";
// 	return r;
// }

#endif // VECTORUTILITIES_H
