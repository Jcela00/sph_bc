#ifndef ASSIGN_NORMALS_HPP
#define ASSIGN_NORMALS_HPP

#include "Definitions.hpp"
#include "VectorUtilities.hpp"
#include "Kernel.hpp"
#include <limits>

// template <typename vd_type, typename NN_type>
// __global__ void AssignNormalsGPU(vd_type vd, NN_type NN)
// {
//     // Key of the particle a
//     unsigned int a;
//     GET_PARTICLE_SORT(a, NN);

//     // if particle FLUID
//     if (vd.template getProp<vd0_type>(a) == FLUID)
//     {
//         Point<DIM, real_number> xa = vd.getPos(a);

//         Point<DIM, real_number> normal_averaged = {0.0};
//         real_number sumW = 0.0;

//         auto Np = NN.getNNIteratorBox(NN.getCell(xa));

//         while (Np.isNext() == true)
//         {
//             auto b = Np.get_sort();

//             const unsigned int typeb = vd.template getProp<vd0_type>(b);

//             if (typeb != FLUID)
//             {
//                 const Point<DIM, real_number> xb = vd.getPos(b);
//                 const Point<DIM, real_number> dr = xa - xb;
//                 const real_number r2 = norm2(dr);

//                 if (r2 < _params_gpu_.r_cut2)
//                 {
//                     real_number W = Wab(sqrt(r2), _params_gpu_.dp, _params_gpu_.Kquintic);
//                     Point<DIM, real_number> normal_b = vd.template getProp<vd8_normal>(b);
//                     normal_averaged += normal_b * W;
//                     sumW += W;
//                 }
//             }

//             ++Np;
//         }

//         if (sumW != 0.0)
//         {
//             normal_averaged = normal_averaged / sumW;
//             real_number norm = getVectorNorm(normal_averaged);
//             if (norm > 0.1)
//             {
//                 for (unsigned int xyz = 0; xyz < DIM; xyz++)
//                 {
//                     vd.template getProp<vd8_normal>(a)[xyz] = normal_averaged.get(xyz) / norm;
//                 }
//             }
//             else
//             {
//                 for (unsigned int xyz = 0; xyz < DIM; xyz++)
//                 {
//                     vd.template getProp<vd8_normal>(a)[xyz] = 0.0;
//                 }
//             }
//         }
//         else
//         {
//             for (unsigned int xyz = 0; xyz < DIM; xyz++)
//             {
//                 vd.template getProp<vd8_normal>(a)[xyz] = 0.0;
//             }
//         }
//     }
// }

// CLOSEST VERSION
template <typename vd_type, typename NN_type>
__global__ void AssignNormalsGPU(vd_type vd, NN_type NN)
{
    // Key of the particle a
    unsigned int a;
    GET_PARTICLE_SORT(a, NN);

    // if particle FLUID
    if (vd.template getProp<vd0_type>(a) == FLUID)
    {
        Point<DIM, real_number> xa = vd.getPos(a);

        constexpr real_number maxreal = std::numeric_limits<real_number>::max();
        real_number DistMin = maxreal;

        auto Np = NN.getNNIteratorBox(NN.getCell(xa));

        while (Np.isNext() == true)
        {
            auto b = Np.get_sort();

            const unsigned int typeb = vd.template getProp<vd0_type>(b);

            if (typeb != FLUID)
            {
                const Point<DIM, real_number> xb = vd.getPos(b);
                const Point<DIM, real_number> dr = xa - xb;
                const real_number r2 = norm2(dr);

                if (r2 < _params_gpu_.r_cut2)
                {
                    // real_number normalDist = std::max(0.25 * _params_gpu_.dp, abs(dotProduct(dr, normal_b)));

                    // real_number normalDist = r2

                    // MIN VERSION
                    if (r2 < DistMin)
                    {
                        vd.template getProp<vd8_normal>(a) = vd.template getProp<vd8_normal>(b);
                        DistMin = r2;
                    }
                }
            }

            ++Np;
        }
        if (DistMin == maxreal) // fluid is not close to any boundary
        {
            for (unsigned int xyz = 0; xyz < DIM; xyz++)
            {
                vd.template getProp<vd8_normal>(a)[xyz] = 0.0;
            }
        }
    }
}

template <typename vd_type, typename NN_type>
inline void AssignNormals(vd_type &vd, NN_type &NN)
{
    // This function computes the density of particles from the summation formulation

    auto it = vd.getDomainIteratorGPU();

    vd.template updateCellListGPU<vd0_type, vd8_normal>(NN);
    CUDA_LAUNCH(AssignNormalsGPU, it, vd.toKernel(), NN.toKernel());
    vd.template restoreOrder<vd0_type, vd8_normal>(NN);
}

// inline __device__ __host__ bool NormalsInteract(Point<DIM, real_number> na, Point<DIM, real_number> nb, Point<DIM, real_number> rab, real_number cos_angle)
// {
//     // Check if either vector is zero
//     if (isZeroVector(na) || isZeroVector(nb))
//         return true;

//     // if (dotProduct(na, nb) > cos_angle)
//     // {
//     //     return true;
//     // }
//     // else
//     // {
//     //     if (dotProduct(na, rab) > 0)
//     //         return true;
//     //     else
//     //         return false;s
//     // }

//     return (dotProduct(na, nb) > cos_angle) || (dotProduct(na, rab) > 0); // simplified version
// }

inline __device__ __host__ bool NormalsInteract(Point<DIM, real_number> na, Point<DIM, real_number> nb, unsigned int typeb, Point<DIM, real_number> rba, real_number cos_angle)
{
    // if (typeb == FLUID)
    // {
    //     // rba = rba / getVectorNorm(rba);
    //     // return isZeroVector(na) || isZeroVector(nb) || (dotProduct(na, nb) > cos_angle) || (dotProduct(na, rba) > cos_angle) || (dotProduct(nb, -rba) > cos_angle);
    //     return true;
    // }
    // else
    //     return isZeroVector(na) || isZeroVector(nb) || (dotProduct(na, nb) > cos_angle);

    return true;
}
inline __device__ __host__ bool NormalsInteractCurv(Point<DIM, real_number> na, Point<DIM, real_number> nb)
{
    // return (dotProduct(na, nb) > 0.0);
    return true;
}
#endif // ASSIGN_NORMALS_HPP
