#include "Physics.hpp"

void EqState(particles &vd, const double rho0, const double B, const double gamma, const double xi)
{
    auto it = vd.getDomainIterator();

    while (it.isNext())
    {
        auto a = it.get();

        if (vd.getProp<type>(a) == FLUID)
        {
            double rho_a = vd.template getProp<rho>(a);
            double rho_frac = rho_a / rho0;

            vd.template getProp<pressure>(a) = B * (std::pow(rho_frac, gamma) - 1.0) + xi;
        }

        ++it;
    }
}

double InvEqState_particle(const double p, const double rho0, const double B, const double gamma, const double xi)
{
    return rho0 * std::pow(((p - xi) / B + 1.0), 1.0 / gamma);
}

double EqState_particle(const double rho, const double rho0, const double B, const double gamma, const double xi)
{
    return B * (std::pow(rho / rho0, gamma) - 1.0) + xi;
}

Point<DIM, double> Pi_physical(const Point<DIM, double> &dr, const double &r, const Point<DIM, double> &dv, const Point<DIM, double> &dW, const double eta)
{
    return eta * (dv * dotProduct(dr, dW)) / (r * r);
}

double PressureForce(const double &rhoa, const double &rhob, const double &prsa, const double &prsb)
{
    return -1.0 * (rhob * prsa + rhoa * prsb) / (rhoa + rhob);
}