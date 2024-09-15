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
