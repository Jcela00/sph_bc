#ifndef KERNEL_HPP
#define KERNEL_HPP
#include "Definitions.hpp"

double Wab(double r, const double H, const double Kquintic);
Point<DIM, double> DWab(const Point<DIM, double> &dx, const double r, const double H, const double Kquintic);

#endif // KERNEL_HPP