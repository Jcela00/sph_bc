#ifndef ASSIGN_NORMALS_HPP
#define ASSIGN_NORMALS_HPP

#include "Definitions.hpp"
#include "VectorUtilities.hpp"
#include "Kernel.hpp"

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

        real_number r2min = 1e10;

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
                    if (r2 < r2min)
                    {
                        vd.template getProp<vd8_normal>(a) = vd.template getProp<vd8_normal>(b);
                        r2min = r2;
                    }
                }
            }
            ++Np;
        }
        if (r2min == 1e10) // fluid is not close to any boundary
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

inline __device__ __host__ bool NormalsInteract(Point<DIM, real_number> na, Point<DIM, real_number> nb, real_number cos_angle)
{
    // Check if either vector is zero
    if (isZeroVector(na) || isZeroVector(nb))
        return true;

    // Check if the angle between the vectors is within the threshold
    return dotProduct(na, nb) > cos_angle;
}

#endif // ASSIGN_NORMALS_HPP
