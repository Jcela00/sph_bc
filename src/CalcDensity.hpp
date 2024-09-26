#ifndef CALC_DENSITY_HPP
#define CALC_DENSITY_HPP

#include "Definitions.hpp"
#include "VectorUtilities.hpp"
#include "Kernel.hpp"

template <typename vd_type, typename NN_type>
__global__ void CalcDensityGPU(vd_type vd, NN_type NN, const Parameters params)
{
    // Key of the particle a
    auto a = GET_PARTICLE(vd);

    // if particle FLUID
    if (vd.template getProp<type>(a) == FLUID)
    {

        // Get the position xb of the particle a
        Point<DIM, real_number> xa = vd.getPos(a);

        // initialize density sum
        real_number rho_sum = 0.0;
        auto Np = NN.getNNIteratorBox(NN.getCell(xa));

        // iterate the neighborhood particles
        while (Np.isNext() == true)
        {
            // Key of b particle
            auto b = Np.get();

            // Get the position xb of the particle b
            const Point<DIM, real_number> xb = vd.getPos(b);

            // Get the vector pointing at xa from xb
            const Point<DIM, real_number> dr = xa - xb;

            // take the norm squared of this vector
            const real_number r2 = norm2(dr);
            const real_number r = sqrt(r2);

            // If the particles interact ...
            if (r < params.r_threshold)
            {

                if (vd.template getProp<type>(b) == FLUID)
                {

                    // evaluate kernel
                    const real_number w = Wab(r, params.H, params.Kquintic);
                    rho_sum += w * params.MassFluid;
                }
                else // if boundary particle
                {
                    if (params.BC_TYPE == NO_SLIP)
                    {
                        const real_number r = sqrt(r2);
                        // evaluate kernel
                        const real_number w = Wab(r, params.H, params.Kquintic);
                        rho_sum += w * params.MassBound;
                    }
                    else if (params.BC_TYPE == NEW_NO_SLIP) // need to evaluate kernel at virtual particles
                    {
                        // get normal vector of b
                        const Point<DIM, real_number> normal = vd.template getProp<normal_vector>(b);
                        // get volumes, and curvature of virtual particles
                        const Point<3, real_number> vol = vd.template getProp<vd_volume>(b);
                        // const real_number kappa = vd.template getProp<curvature_boundary>(b);

                        // Apply offsets to dr to get 3 vectrors pointing to virtual particles
                        const std::array<Point<DIM, real_number>, 3> R_virtual = getBoundaryPositions(-1.0 * dr, normal, params.dp);

                        // iterate the 3 virtual particles
                        for (int i = 0; i < 3; i++)
                        {
                            const real_number W = Wab(getVectorNorm(R_virtual[i]), params.H, params.Kquintic);
                            rho_sum += W * vol[i] * params.rho_zero; // this is equivalent to +=W*mass
                        }
                    }
                }
            }

            ++Np;
        }
        if (rho_sum != 0.0)
        {
            vd.template getProp<rho>(a) = rho_sum;
        }
        else
        {
            vd.template getProp<rho>(a) = params.rho_zero;
        }
    }
}

template <typename vd_type, typename CellList>
inline void CalcDensity(vd_type &vd, CellList &NN, const Parameters &params)
{
    // This function computes the density of particles from the summation formulation

    auto it = vd.getDomainIteratorGPU();

    // Update the cell-list
    vd.updateCellListGPU(NN);

    CUDA_LAUNCH(CalcDensityGPU, it, vd.toKernel(), NN.toKernel(), params);
}

#endif
