#ifndef VECTORUTILITIES_H
#define VECTORUTILITIES_H

#include <array>
#include <cmath> // sqrt
#include "Definitions.hpp"

double dotProduct(const Point<DIM, double> &v, const Point<DIM, double> &w);

std::array<Point<DIM, double>, DIM> dyadicProduct(const Point<DIM, double> &v, const Point<DIM, double> &w);

Point<DIM, double> matVec(const std::array<Point<DIM, double>, DIM> &m, const Point<DIM, double> &v);

double getVectorNorm(const Point<DIM, double> &v);

void normalizeVector(Point<DIM, double> &v);

Point<DIM, double> getPerpendicularUnit2D(const Point<DIM, double> &v);

std::array<Point<DIM, double>, 3> getBoundaryPositions(const Point<DIM, double> &r, const Point<DIM, double> &normal, double dp);

Point<DIM, double> crossProduct(const Point<DIM, double> &v, const Point<DIM, double> &w);

void ApplyRotation(Point<DIM, double> &x, const double theta, const Point<DIM, double> centre);

// function to know the type of a variable
// used to know types of auto variables
template <class T>
std::string type_name()
{
	typedef typename std::remove_reference<T>::type TR;
	std::unique_ptr<char, void (*)(void *)> own(
#ifndef _MSC_VER
		abi::__cxa_demangle(typeid(TR).name(), nullptr,
							nullptr, nullptr),
#else
		nullptr,
#endif
		std::free);
	std::string r = own != nullptr ? own.get() : typeid(TR).name();
	if (std::is_const<TR>::value)
		r += " const";
	if (std::is_volatile<TR>::value)
		r += " volatile";
	if (std::is_lvalue_reference<T>::value)
		r += "&";
	else if (std::is_rvalue_reference<T>::value)
		r += "&&";
	return r;
}

#endif // VECTORUTILITIES_H
