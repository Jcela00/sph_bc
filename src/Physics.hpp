#ifndef PHYSICS_HPP
#define PHYSICS_HPP

#include "Definitions.hpp"
#include "VectorUtilities.hpp"

void EqState(particles &vd, const double rho0, const double B, const double gamma, const double xi);

double InvEqState_particle(const double p, const double rho0, const double B, const double gamma, const double xi);

double EqState_particle(const double rho, const double rho0, const double B, const double gamma, const double xi);

Point<DIM, double> Pi_physical(const Point<DIM, double> &dr, const double &r, const Point<DIM, double> &dv, const Point<DIM, double> &dW, const double eta);

double PressureForce(const double &rhoa, const double &rhob, const double &prsa, const double &prsb);

#endif // PHYSICS_HPP
