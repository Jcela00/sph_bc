#ifndef EXTRAPOLATEVELOCITY_HPP
#define EXTRAPOLATEVELOCITY_HPP

#include "Definitions.hpp"
#include "VectorUtilities.hpp"
#include "Kernel.hpp"

template <typename vd_type, typename NN_type>
__global__ void ExtrapolateVelocityGPU(vd_type vd, NN_type NN)
{

    // Key of the b particle
    unsigned int b;
    GET_PARTICLE_SORT(b, NN);

    // if particle boundary
    if (vd.template getProp<vd0_type>(b) != FLUID)
    {

        // Get the position xb of the boundary particle
        Point<DIM, real_number> xb = vd.getPos(b);

        // initialize sums
        Point<DIM, real_number> sum_vW = (DIM == 2) ? Point<DIM, real_number>{0.0, 0.0} : Point<DIM, real_number>{0.0, 0.0, 0.0};

        real_number sum_pW = 0.0;
        real_number sum_W = 0.0;

        auto Np = NN.getNNIteratorBox(NN.getCell(xb));

        // iterate the neighborhood fluid particles
        while (Np.isNext() == true)
        {
            // Key of fluid particle
            auto f = Np.get_sort();

            // if (b == f) skip this particle
            if (b == f)
            {
                ++Np;
                continue;
            }

            // Skip other boundary particles
            if (vd.template getProp<vd0_type>(f) != FLUID)
            {
                ++Np;
                continue;
            }

            // Get the position xf of the fluid particle
            Point<DIM, real_number> xf = vd.getPos(f);

            // Get the velocity of the fluid particle
            Point<DIM, real_number> vf = vd.template getProp<vd4_velocity>(f);

            // Get the density of the fluid particle
            real_number rhof = vd.template getProp<vd1_rho>(f);

            // Get the pressure of the fluid particle
            real_number Pf = vd.template getProp<vd2_pressure>(f);

            // Get the vector pointing at xb from xf rwf
            Point<DIM, real_number> dr = xb - xf;

            // take the norm squared of this vector
            real_number r2 = norm2(dr);

            // If the particles interact ...
            if (r2 < _params_gpu_.r_cut2)
            {
                // calculate distance
                real_number r = sqrt(r2);

                // evaluate kernel
                real_number w = Wab(r, _params_gpu_.H, _params_gpu_.Kquintic);
                // compute v*W
                sum_vW += w * vf;

                // compute Pf +rhof*(g-a) dot rwf
                // at the moment a = 0
                const real_number dot = dotProduct(dr, _params_gpu_.gravity_vector);
                sum_pW += w * (Pf + rhof * dot);
                sum_W += w;
            }

            ++Np;
        }
        if (sum_W != 0.0)
        {
            // Set the velocity of the boundary particle b ( to be used in the momentum equation to impose BC)
            // We use the v_transport field because boundary particles dont use this array

            // solid particles store the centre of rotation in the force transport field since it is unused
            // vector pointing from centre of rotation to marker particle
            const Point<DIM, real_number> radial_vec = {xb.get(0) - vd.template getProp<vd7_force_t>(b)[0], xb.get(1) - vd.template getProp<vd7_force_t>(b)[1]};
            const real_number radial = getVectorNorm(radial_vec);
            // vector tangential to the radial vector, rotation velocity is in this direction
            const Point<DIM, real_number> rotation_tangential = getPerpendicularUnit2D(radial_vec);
            Point<DIM, real_number> vw = vd.template getProp<vd4_velocity>(b);
            vw.get(0) += radial * vd.template getProp<vd10_omega>(b) * rotation_tangential.get(0);
            vw.get(1) += radial * vd.template getProp<vd10_omega>(b) * rotation_tangential.get(1);

            for (int xyz = 0; xyz < DIM; ++xyz)
            {
                vd.template getProp<vd5_velocity_t>(b)[xyz] = 2.0 * vw[xyz] - sum_vW.get(xyz) / sum_W;
            }

            // Set the pressure of the boundary particle b
            vd.template getProp<vd2_pressure>(b) = sum_pW / sum_W;
            // Compute density from inverted Eq of state
            vd.template getProp<vd1_rho>(b) = InvEqState_particle(vd.template getProp<vd2_pressure>(b), _params_gpu_.rho0, _params_gpu_.B, _params_gpu_.gamma, _params_gpu_.xi);
        }
        else
        {

            for (int xyz = 0; xyz < DIM; ++xyz)
            {
                vd.template getProp<vd5_velocity_t>(b)[xyz] = 2.0 * vd.template getProp<vd4_velocity>(b)[xyz];
            }

            vd.template getProp<vd2_pressure>(b) = 0.0;
            vd.template getProp<vd1_rho>(b) = _params_gpu_.rho0;
        }
    }
}

template <typename vd_type, typename NN_type>
inline void ExtrapolateVelocity(vd_type &vd, NN_type &NN)
{
    // This function fills the value of v_transport for the boundary particles,
    // with the extrapolated velocity for the old no slip BC

    auto it = vd.getDomainIteratorGPU();

    vd.template updateCellListGPU<vd0_type, vd1_rho, vd2_pressure, vd4_velocity, vd7_force_t, vd10_omega>(NN);

    CUDA_LAUNCH(ExtrapolateVelocityGPU, it, vd.toKernel(), NN.toKernel());

    vd.template restoreOrder<vd0_type, vd1_rho, vd2_pressure, vd4_velocity, vd7_force_t, vd10_omega, vd5_velocity_t>(NN);
}

#endif // EXTRAPOLATEVEVELOCITY_HPP