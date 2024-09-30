#ifndef KERNEL_HPP
#define KERNEL_HPP
#include "Definitions.hpp"

// real_number Wab(real_number r, const real_number H, const real_number Kquintic);
// Point<DIM, real_number> DWab(const Point<DIM, real_number> &dx, const real_number r, const real_number H, const real_number Kquintic);

// inline __device__ __host__ real_number Wab(real_number r, const real_number H, const real_number Kquintic)
// {
//     const real_number q = r / H;
//     const real_number tmp3 = (3.0 - q) * (3.0 - q) * (3.0 - q) * (3.0 - q) * (3.0 - q);
//     const real_number tmp2 = (2.0 - q) * (2.0 - q) * (2.0 - q) * (2.0 - q) * (2.0 - q);
//     const real_number tmp1 = (1.0 - q) * (1.0 - q) * (1.0 - q) * (1.0 - q) * (1.0 - q);

//     if (q < 0.0)
//         return 0.0;
//     else if (q < 1.0)
//         return Kquintic * (tmp3 - 6.0 * tmp2 + 15.0 * tmp1);
//     else if (q < 2.0)
//         return Kquintic * (tmp3 - 6.0 * tmp2);
//     else if (q < 3.0)
//         return Kquintic * tmp3;
//     else
//         return 0.0;
// }

// inline __device__ __host__ Point<DIM, real_number> DWab(const Point<DIM, real_number> &dx, const real_number r, const real_number H, const real_number Kquintic)
// {
//     Point<DIM, real_number> DW;
//     const real_number q = r / H;

//     const real_number tmp3 = (3.0 - q) * (3.0 - q) * (3.0 - q) * (3.0 - q);
//     const real_number tmp2 = (2.0 - q) * (2.0 - q) * (2.0 - q) * (2.0 - q);
//     const real_number tmp1 = (1.0 - q) * (1.0 - q) * (1.0 - q) * (1.0 - q);

//     const real_number c1 = -5.0 * Kquintic / H;
//     real_number factor;

//     if (q < 1.0)
//     {
//         if (q > 0.0)
//             factor = c1 * (tmp3 - 6.0 * tmp2 + 15.0 * tmp1) / r;
//         else
//             factor = 0.0;
//     }
//     else if (q < 2.0)
//         factor = c1 * (tmp3 - 6.0 * tmp2) / r;
//     else if (q < 3.0)
//         factor = c1 * tmp3 / r;
//     else
//         factor = 0.0;

//     DW = factor * dx;

//     return DW;
// }

inline __device__ __host__ real_number Wab(const real_number r, const real_number H, const real_number Kquintic)
{
    const real_number q = r / H;
    const real_number tmp3 = (3.0 - q) * (3.0 - q) * (3.0 - q) * (3.0 - q) * (3.0 - q);
    const real_number tmp2 = (2.0 - q) * (2.0 - q) * (2.0 - q) * (2.0 - q) * (2.0 - q);
    const real_number tmp1 = (1.0 - q) * (1.0 - q) * (1.0 - q) * (1.0 - q) * (1.0 - q);

    const real_number switch3 = (q < 3.0f) ? 1.0f : 0.0f;
    const real_number switch2 = (q < 2.0f) ? 1.0f : 0.0f;
    const real_number switch1 = (q < 1.0f) ? 1.0f : 0.0f;

    return Kquintic * (switch3 * tmp3 - 6.0f * switch2 * tmp2 + 15.0f * switch1 * tmp1);
}

inline __device__ __host__ Point<DIM, real_number> DWab(const Point<DIM, real_number> &dx, const real_number r, const real_number H, const real_number Kquintic)
{
    Point<DIM, real_number> DW;
    const real_number q = r / H;

    const real_number tmp3 = (3.0f - q) * (3.0f - q) * (3.0f - q) * (3.0f - q);
    const real_number tmp2 = (2.0f - q) * (2.0f - q) * (2.0f - q) * (2.0f - q);
    const real_number tmp1 = (1.0f - q) * (1.0f - q) * (1.0f - q) * (1.0f - q);

    const real_number switch3 = (q < 3.0f) ? 1.0f : 0.0f;
    const real_number switch2 = (q < 2.0f) ? 1.0f : 0.0f;
    const real_number switch1 = (q < 1.0f) ? 1.0f : 0.0f;

    real_number factor = (q != 0.0f) ? (-5.0 * Kquintic / H) * (switch3 * tmp3 - 6.0f * switch2 * tmp2 + 15.0f * switch1 * tmp1) / r : 0.0f;

    DW = factor * dx;

    return DW;
}

#endif // KERNEL_HPP