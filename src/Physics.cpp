#include "Physics.hpp"

void EqState(particles &vd, const real_number rho0, const real_number B, const real_number gamma, const real_number xi)
{
    auto it = vd.getDomainIterator();

    while (it.isNext())
    {
        auto a = it.get();

        if (vd.getProp<type>(a) == FLUID)
        {
            real_number rho_a = vd.template getProp<rho>(a);
            real_number rho_frac = rho_a / rho0;

            vd.template getProp<pressure>(a) = B * (std::pow(rho_frac, gamma) - 1.0) + xi;
        }

        ++it;
    }
}
