#ifndef PHYSICS_HPP
#define PHYSICS_HPP

#include "Definitions.hpp"
#include "VectorUtilities.hpp"

// Declaration
void EqState(particles &vd, const double rho0, const double B, const double gamma, const double xi);

// Declaration And Definition because of inline
inline double InvEqState_particle(const double p, const double rho0, const double B, const double gamma, const double xi)
{
    return rho0 * std::pow(((p - xi) / B + 1.0), 1.0 / gamma);
}

inline double EqState_particle(const double rho, const double rho0, const double B, const double gamma, const double xi)
{
    return B * (std::pow(rho / rho0, gamma) - 1.0) + xi;
}

inline Point<DIM, double> Pi_physical(const Point<DIM, double> &dr, const double &r, const Point<DIM, double> &dv, const Point<DIM, double> &dW, const double eta)
{
    return eta * (dv * dotProduct(dr, dW)) / (r * r);
}

inline double PressureForce(const double &rhoa, const double &rhob, const double &prsa, const double &prsb)
{
    return -1.0 * (rhob * prsa + rhoa * prsb) / (rhoa + rhob);
}

#endif // PHYSICS_HPP
