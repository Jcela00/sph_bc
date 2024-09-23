#ifndef CALC_FORCES_HPP
#define CALC_FORCES_HPP

#include "Definitions.hpp"
#include "Interactions.hpp"

template <typename CellList>
void calc_forces(particles &vd,
                 CellList &NN,
                 real_number &obstacle_force_x,
                 real_number &obstacle_force_y,
                 bool calc_drag,
                 const Parameters &params)
{
    // get domain iterator
    vector_dist_iterator part = vd.getDomainIterator();

    const real_number dp = params.dp;

    // Update the cell-list
    vd.updateCellList(NN);

    // For each particle ...
    while (part.isNext())
    {
        // get key of particle a
        auto a = part.get();

        // if the particle is FLUID
        if (vd.getProp<type>(a) == FLUID)
        {

            const Point<DIM, real_number> xa = vd.getPos(a);
            const real_number massa = params.MassFluid;
            const real_number rhoa = vd.getProp<rho>(a);
            const real_number Pa = vd.getProp<pressure>(a);
            const Point<DIM, real_number> va = vd.getProp<velocity>(a);

            // Reset the force counter (0 + gravity)
            for (int xyz = 0; xyz < DIM; xyz++)
            {
                vd.template getProp<force>(a)[xyz] = params.gravity_vector.get(xyz);
                vd.template getProp<force_transport>(a)[xyz] = 0.0;
            }

            // Reset the density difference
            vd.template getProp<drho>(a) = 0.0;

            // Get an iterator over the neighborhood particles of p
            auto Np = NN.getNNIteratorBox(NN.getCell(xa));

            // For each neighborhood particle
            while (Np.isNext() == true)
            {
                // get key of particle b
                unsigned long b = Np.get();

                // if (a == b) skip this particle
                if (a.getKey() == b)
                {
                    ++Np;
                    continue;
                }

                // Get the position xb of the particle
                Point<DIM, real_number> xb = vd.getPos(b);

                // Get the distance between a and b
                // in fluid - boundary its xf-xb i.e. vector pointing at fluid from boundary
                Point<DIM, real_number> dr = xa - xb;

                // take the norm (squared) of this vector
                real_number r2 = norm2(dr);

                // in general we do not accumulate force, only in the case of interacting with an obstacle
                bool accumulate_force = false;

                // if they interact
                if (r2 < params.r_threshold * params.r_threshold)
                {
                    if (vd.getProp<type>(b) == BOUNDARY)
                    {
                        if (calc_drag && ((xb.get(1) > 0.0 && xb.get(1) < params.length[1])))
                            accumulate_force = true;
                        else
                            accumulate_force = false;

                        if (params.BC_TYPE == NO_SLIP)
                            interact_fluid_boundary_old(vd, a, massa, rhoa, Pa, xa, va, xb, dr, r2, b, accumulate_force, obstacle_force_x, obstacle_force_y, params);
                        else if (params.BC_TYPE == NEW_NO_SLIP)
                            interact_fluid_boundary_new(vd, a, massa, rhoa, Pa, xa, va, xb, dr, b, accumulate_force, obstacle_force_x, obstacle_force_y, params);
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

template <typename CellList>
void calc_forces_regularize(particles &vd,
                            CellList &NN,
                            const Parameters &params)
{
    // get domain iterator
    vector_dist_iterator part = vd.getDomainIterator();

    // Update the cell-list
    vd.updateCellList(NN);

    // For each particle ...
    while (part.isNext())
    {
        // get key of particle a
        auto a = part.get();

        // if the particle is FLUID
        if (vd.getProp<type>(a) == FLUID)
        {
            const Point<DIM, real_number> xa = vd.getPos(a);
            const real_number massa = params.MassFluid;
            const real_number rhoa = vd.getProp<rho>(a);

            // Reset the force counter (0 + gravity)
            for (int xyz = 0; xyz < DIM; xyz++)
            {
                vd.template getProp<force_transport>(a)[xyz] = 0.0;
            }

            // Reset the density difference
            vd.template getProp<drho>(a) = 0.0;

            // Get an iterator over the neighborhood particles of p
            auto Np = NN.getNNIteratorBox(NN.getCell(xa));

            // For each neighborhood particle
            while (Np.isNext() == true)
            {
                // get key of particle b
                unsigned long b = Np.get();

                // if (a == b) skip this particle
                if (a.getKey() == b)
                {
                    ++Np;
                    continue;
                }

                // Get the position xb of the particle
                Point<DIM, real_number> xb = vd.getPos(b);

                // Get the distance between a and b
                // in fluid - boundary its xf-xb i.e. vector pointing at fluid from boundary
                Point<DIM, real_number> dr = xa - xb;

                // take the norm (squared) of this vector
                real_number r2 = norm2(dr);

                // if they interact
                if (r2 < params.r_threshold * params.r_threshold)
                {
                    if (vd.getProp<type>(b) == BOUNDARY)
                    {
                        if (params.BC_TYPE == NO_SLIP)
                        {
                            // interact_fluid_boundary_old(vd, a, massa, rhoa, Pa, xa, va, xb, dr, r2, b, accumulate_force, obstacle_force_x, obstacle_force_y, params);
                            // NOT READY YET
                        }
                        else if (params.BC_TYPE == NEW_NO_SLIP)
                            interact_fluid_boundary_new_regularize(vd, a, rhoa, dr, b, params);
                    }
                    else
                    {
                        interact_fluid_fluid_regularize(vd, a, massa, rhoa, dr, r2, b, params);
                    }
                }
                ++Np;
            }
        }

        ++part;
    }
}

#endif // CALC_FORCES_HPP
