#ifndef PHYSICS_HPP
#define PHYSICS_HPP

#include "Definitions.hpp"
#include "VectorUtilities.hpp"

// Declaration
void EqState(particles &vd, const real_number rho0, const real_number B, const real_number gamma, const real_number xi);

// Declaration And Definition because of inline
inline real_number InvEqState_particle(const real_number p, const real_number rho0, const real_number B, const real_number gamma, const real_number xi)
{
    return rho0 * std::pow(((p - xi) / B + 1.0), 1.0 / gamma);
}

inline real_number EqState_particle(const real_number rho, const real_number rho0, const real_number B, const real_number gamma, const real_number xi)
{
    return B * (std::pow(rho / rho0, gamma) - 1.0) + xi;
}

inline Point<DIM, real_number> Pi_physical(const Point<DIM, real_number> &dr, const real_number &r, const Point<DIM, real_number> &dv, const Point<DIM, real_number> &dW, const real_number eta)
{
    real_number denominator = r * r;
    Point<DIM, real_number> result = (eta * (dv * dotProduct(dr, dW)));
    result = result / denominator;
    // return eta * (dv * dotProduct(dr, dW)) / (r * r); // USED to be like this but when changing to real_number the quotient cant be done in one line

    return result;
}

inline real_number PressureForce(const real_number &rhoa, const real_number &rhob, const real_number &prsa, const real_number &prsb)
{
    return -1.0 * (rhob * prsa + rhoa * prsb) / (rhoa + rhob);
}

#endif // PHYSICS_HPP
