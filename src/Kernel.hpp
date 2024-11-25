#ifndef KERNEL_HPP
#define KERNEL_HPP
#include "Definitions.hpp"

inline __device__ __host__ real_number Wab(const real_number r, const real_number H, const real_number Kquintic)
{
    const real_number q = r / H;
    const real_number tmp3 = (3.0 - q) * (3.0 - q) * (3.0 - q) * (3.0 - q) * (3.0 - q);
    const real_number tmp2 = (2.0 - q) * (2.0 - q) * (2.0 - q) * (2.0 - q) * (2.0 - q);
    const real_number tmp1 = (1.0 - q) * (1.0 - q) * (1.0 - q) * (1.0 - q) * (1.0 - q);

    const real_number switch3 = (q < 3.0) ? 1.0 : 0.0;
    const real_number switch2 = (q < 2.0) ? 1.0 : 0.0;
    const real_number switch1 = (q < 1.0) ? 1.0 : 0.0;

    return Kquintic * (switch3 * tmp3 - 6.0 * switch2 * tmp2 + 15.0 * switch1 * tmp1);
}

inline __device__ __host__ Point<DIM, real_number> DWab(const Point<DIM, real_number> &dx, const real_number r, const real_number H, const real_number Kquintic)
{
    Point<DIM, real_number> DW;
    const real_number q = r / H;

    const real_number tmp3 = (3.0 - q) * (3.0 - q) * (3.0 - q) * (3.0 - q);
    const real_number tmp2 = (2.0 - q) * (2.0 - q) * (2.0 - q) * (2.0 - q);
    const real_number tmp1 = (1.0 - q) * (1.0 - q) * (1.0 - q) * (1.0 - q);

    const real_number switch3 = (q < 3.0) ? 1.0 : 0.0;
    const real_number switch2 = (q < 2.0) ? 1.0 : 0.0;
    const real_number switch1 = (q < 1.0) ? 1.0 : 0.0;

    real_number factor = (q != 0.0) ? (-5.0 * Kquintic / H) * (switch3 * tmp3 - 6.0 * switch2 * tmp2 + 15.0 * switch1 * tmp1) / r : 0.0;

    DW = factor * dx;

    return DW;
}

#endif // KERNEL_HPP