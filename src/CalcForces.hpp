#ifndef CALC_FORCES_HPP
#define CALC_FORCES_HPP

#include "Definitions.hpp"
#include "Interactions.hpp"

template <typename CellList>
void calc_forces(particles &vd,
                 CellList &NN,
                 double &obstacle_force_x,
                 double &obstacle_force_y,
                 bool calc_drag,
                 const Parameters &params)
{
    vector_dist_iterator part = vd.getDomainIterator();

    const double dp = params.dp;
    // const double max_curvature = 1.0 / (3.0 * dp);

    // Update the cell-list
    vd.updateCellList(NN);

    // For each particle ...
    while (part.isNext())
    {
        // get key of particle a

        auto a = part.get();

        // We threat FLUID particle differently from BOUNDARY PARTICLES
        // Boundary particles dont need any force computation

        if (vd.getProp<type>(a) == FLUID) // INTERACTION OF A FLUID PARTICLE WITH ITS NEIGHBORHOOD
        {

            // Get the position xp of the particle
            Point<DIM, double> xa = vd.getPos(a);

            // Take the mass of the particle dependently if it is FLUID or BOUNDARY
            double massa = params.MassFluid;

            // Get the density and pressure of the particle a
            double rhoa = vd.getProp<rho>(a);
            double Pa = vd.getProp<pressure>(a);

            // Get the Velocity of the particle a
            Point<DIM, double> va = vd.getProp<velocity>(a);

            // Reset the force counter (0 + gravity)
            for (int xyz = 0; xyz < DIM; xyz++)
            {
                vd.template getProp<force>(a)[xyz] = params.gravity_vector.get(xyz);
                vd.template getProp<force_transport>(a)[xyz] = 0.0;
            }
            vd.template getProp<drho>(a) = 0.0;

            // Get an iterator over the neighborhood particles of p
            auto Np = NN.getNNIterator(NN.getCell(xa));

            // For each neighborhood particle
            while (Np.isNext() == true)
            {
                // get particle b
                unsigned long b = Np.get();

                // if (a == b) skip this particle
                if (a.getKey() == b)
                {
                    ++Np;
                    continue;
                };

                // Get the position xb of the particle
                Point<DIM, double> xb = vd.getPos(b);

                // Get the distance between a and b
                // in fluid - boundary its xf-xb i.e. vector pointing at fluid from boundary
                Point<DIM, double> dr = xa - xb;

                // take the norm (squared) of this vector
                double r2 = norm2(dr);
                bool accumulate_force = true;
                // if they interact
                if (r2 < params.r_threshold * params.r_threshold)
                {
                    if (vd.getProp<type>(b) == BOUNDARY)
                    {
                        if (calc_drag && (xb.get(1) > 2.0 * dp || xb.get(1) < params.length[1] - 2.0 * dp)) // to exclude solid walls
                            accumulate_force = true;
                        else
                            accumulate_force = false;

                        if (params.BC_TYPE == NO_SLIP)
                            interact_fluid_boundary_old(vd, a, massa, rhoa, Pa, xa, va, xb, dr, r2, b, accumulate_force, obstacle_force_x, obstacle_force_y, params);
                        else if (params.BC_TYPE == NEW_NO_SLIP)
                        {
                            // const Point<DIM, double> normal_b = vd.getProp<normal_vector>(b);
                            // const double kappa = vd.getProp<curvature_boundary>(b);
                            // if (dotProduct(normal_b, dr) > (max_cos / max_curvature) * kappa)
                            // if (dotProduct(normal_b, dr) > 0.0)
                                interact_fluid_boundary_new(vd, a, massa, rhoa, Pa, xa, va, xb, dr, b, accumulate_force, obstacle_force_x, obstacle_force_y, params);
                        }
                    }
                    else
                    {
                        interact_fluid_fluid(vd, a, massa, rhoa, Pa, xa, va, xb, dr, r2, b, params);
                    }
                }
                ++Np;
            }
        }

        ++part;
    }
}

#endif // CALC_FORCES_HPP
