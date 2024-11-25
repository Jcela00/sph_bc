#ifndef CALC_DENSITY_HPP
#define CALC_DENSITY_HPP

#include "Definitions.hpp"
#include "VectorUtilities.hpp"
#include "Kernel.hpp"

template <typename vd_type, typename NN_type>
__global__ void CalcDensityGPU_new(vd_type vd, NN_type NN)
{
    // Key of the particle a
    unsigned int a;
    GET_PARTICLE_SORT(a, NN);

    // if particle FLUID
    if (vd.template getProp<vd0_type>(a) == FLUID)
    {
        Point<DIM, real_number> xa = vd.getPos(a);

        real_number rho_sum = 0.0;

        auto Np = NN.getNNIteratorBox(NN.getCell(xa));
        while (Np.isNext() == true)
        {
            auto b = Np.get_sort();

            const Point<DIM, real_number> xb = vd.getPos(b);
            const Point<DIM, real_number> dr = xa - xb;
            const real_number r2 = norm2(dr);

            if (r2 < _params_gpu_.r_cut2)
            {
                const real_number r = sqrt(r2);

                if (vd.template getProp<vd0_type>(b) == FLUID)
                {
                    // evaluate kernel
                    const real_number w = Wab(r, _params_gpu_.H, _params_gpu_.Kquintic);
                    rho_sum += w * _params_gpu_.MassFluid;
                }
                else // if boundary particle
                {
                    const Point<DIM, real_number> normal = vd.template getProp<vd8_normal>(b);
                    const Point<3, real_number> vol = vd.template getProp<vd9_volume>(b);

                    // Apply offsets to dr to get 3 vectrors pointing to virtual particles
                    const std::array<Point<DIM, real_number>, 3> R_virtual = getBoundaryPositions(-1.0 * dr, normal, _params_gpu_.dp);

                    // iterate the 3 virtual particles
                    for (int i = 0; i < 3; i++)
                    {
                        const real_number W = Wab(getVectorNorm(R_virtual[i]), _params_gpu_.H, _params_gpu_.Kquintic);
                        rho_sum += W * vol[i] * _params_gpu_.rho0; // this is equivalent to +=W*mass
                    }
                }
            }

            ++Np;
        }

        vd.template getProp<vd1_rho>(a) = rho_sum;
    }
}

template <typename vd_type, typename NN_type>
__global__ void CalcDensityGPU_old(vd_type vd, NN_type NN)
{
    // Key of the particle a
    unsigned int a;
    GET_PARTICLE_SORT(a, NN);

    // if particle FLUID
    if (vd.template getProp<vd0_type>(a) == FLUID)
    {
        Point<DIM, real_number> xa = vd.getPos(a);

        real_number rho_sum = 0.0;

        auto Np = NN.getNNIteratorBox(NN.getCell(xa));
        while (Np.isNext() == true)
        {
            auto b = Np.get_sort();

            const Point<DIM, real_number> xb = vd.getPos(b);
            const Point<DIM, real_number> dr = xa - xb;
            const real_number r2 = norm2(dr);

            if (r2 < _params_gpu_.r_cut2)
            {
                const real_number r = sqrt(r2);

                const real_number w = Wab(r, _params_gpu_.H, _params_gpu_.Kquintic);
                if (vd.template getProp<vd0_type>(b) == FLUID)
                {
                    rho_sum += w * _params_gpu_.MassFluid;
                }
                else // if boundary particle
                {
                    rho_sum += w * _params_gpu_.MassBound;
                }
            }

            ++Np;
        }

        vd.template getProp<vd1_rho>(a) = rho_sum;
    }
}

template <typename vd_type, typename NN_type>
inline void CalcDensity(vd_type &vd, NN_type &NN, const Parameters &params)
{
    // This function computes the density of particles from the summation formulation

    auto it = vd.getDomainIteratorGPU();

    if (params.BC_TYPE == NO_SLIP)
    {
        vd.template updateCellListGPU<vd0_type, vd1_rho, vd2_pressure, vd4_velocity>(NN);
        CUDA_LAUNCH(CalcDensityGPU_old, it, vd.toKernel(), NN.toKernel());
        vd.template restoreOrder<vd0_type, vd1_rho, vd2_pressure, vd4_velocity>(NN);
    }
    else if (params.BC_TYPE == NEW_NO_SLIP)
    {
        vd.template updateCellListGPU<vd0_type, vd1_rho, vd2_pressure, vd4_velocity, vd8_normal, vd9_volume>(NN);
        CUDA_LAUNCH(CalcDensityGPU_new, it, vd.toKernel(), NN.toKernel());
        vd.template restoreOrder<vd0_type, vd1_rho, vd2_pressure, vd4_velocity, vd8_normal, vd9_volume>(NN);
    }
}

#endif
