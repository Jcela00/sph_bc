#ifndef KERNEL_HPP
#define KERNEL_HPP
#include "Definitions.hpp"

// real_number Wab(real_number r, const real_number H, const real_number Kquintic);
// Point<DIM, real_number> DWab(const Point<DIM, real_number> &dx, const real_number r, const real_number H, const real_number Kquintic);

inline real_number Wab(real_number r, const real_number H, const real_number Kquintic)
{
    const real_number q = r / H;
    const real_number tmp3 = (3.0 - q) * (3.0 - q) * (3.0 - q) * (3.0 - q) * (3.0 - q);
    const real_number tmp2 = (2.0 - q) * (2.0 - q) * (2.0 - q) * (2.0 - q) * (2.0 - q);
    const real_number tmp1 = (1.0 - q) * (1.0 - q) * (1.0 - q) * (1.0 - q) * (1.0 - q);

    if (q < 0.0)
        return 0.0;
    else if (q < 1.0)
        return Kquintic * (tmp3 - 6.0 * tmp2 + 15.0 * tmp1);
    else if (q < 2.0)
        return Kquintic * (tmp3 - 6.0 * tmp2);
    else if (q < 3.0)
        return Kquintic * tmp3;
    else
        return 0.0;
}

inline Point<DIM, real_number> DWab(const Point<DIM, real_number> &dx, const real_number r, const real_number H, const real_number Kquintic)
{
    Point<DIM, real_number> DW;
    const real_number q = r / H;

    const real_number tmp3 = (3.0 - q) * (3.0 - q) * (3.0 - q) * (3.0 - q);
    const real_number tmp2 = (2.0 - q) * (2.0 - q) * (2.0 - q) * (2.0 - q);
    const real_number tmp1 = (1.0 - q) * (1.0 - q) * (1.0 - q) * (1.0 - q);

    const real_number c1 = -5.0 * Kquintic / H;
    real_number factor;

    if (q < 1.0)
    {
        if (q > 0.0)
            factor = c1 * (tmp3 - 6.0 * tmp2 + 15.0 * tmp1) / r;
        else
            factor = 0.0;
    }
    else if (q < 2.0)
        factor = c1 * (tmp3 - 6.0 * tmp2) / r;
    else if (q < 3.0)
        factor = c1 * tmp3 / r;
    else
        factor = 0.0;

    DW = factor * dx;

    return DW;
}

#endif // KERNEL_HPP