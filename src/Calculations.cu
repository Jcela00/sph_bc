#include "Calculations.hpp"

void CalcVolume(particles &vd, real_number dp)
{
    // This function computes the volume of the virtual particles for the new no-slip BC

    auto part = vd.getDomainIterator();

    // For each particle ...
    while (part.isNext())
    {
        // Key of the particle a
        vect_dist_key_dx a = part.get();

        // if particle BOUNDARY
        if (vd.getProp<vd0_type>(a) != FLUID)
        {
            real_number dxwall = vd.getProp<vd9_volume>(a)[0];
            real_number kappa = vd.getProp<vd9_volume>(a)[1];

            for (int i = 0; i < 3; i++)
            {

                // if constexpr (DIM == 2)
                //     vd.template getProp<vd9_volume>(a)[i] = dp * dp;
                // else if constexpr (DIM == 3)
                //     vd.template getProp<vd9_volume>(a)[i] = dp * dp * dp;

                // n=i+1
                real_number n = static_cast<real_number>(i + 1);
                if constexpr (DIM == 2)
                    vd.template getProp<vd9_volume>(a)[i] = std::max(0.0, 0.5 * (2.0 * dp + dp * dp * kappa - 2.0 * n * dp * dp * kappa) * dxwall);
                else if constexpr (DIM == 3)
                    vd.template getProp<vd9_volume>(a)[i] = std::max(0.0, (1.0 / 3.0) * (3.0 * dp + kappa * dp * dp * (3.0 - 6.0 * n) + kappa * kappa * dp * dp * dp * (3 * n * n - 3 * n + 1)) * dxwall);
            }
        }
        ++part;
    }
}
