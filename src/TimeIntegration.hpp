#ifndef TIMEINTEGRATION_HPP
#define TIMEINTEGRATION_HPP

#include "Definitions.hpp"
#include "VectorUtilities.hpp"
#include "Physics.hpp"

void max_velocity(particles &vd, Vcluster<> &v_cl, real_number &max_vel);
real_number calc_deltaT(particles &vd, Vcluster<> &v_cl, const Parameters &params);
Point<DIM, real_number> SolidBodyAcceleration(real_number t, const Parameters &params);

template <typename CellList>
void kick_drift_int(particles &vd,
                    CellList &NN,
                    const real_number dt,
                    Vcluster<> &v_cl,
                    real_number &obstalce_force_x,
                    real_number &obstacle_force_y,
                    bool calc_drag,
                    real_number t,
                    Parameters &params)
{
    // particle iterator
    auto part = vd.getDomainIterator();

    const real_number dt_2 = dt * 0.5;
    // vd.ghost_get<type, rho, pressure, velocity, v_transport, normal_vector, curvature_boundary, arc_length, vd_volume, vd_omega>();

    Point<DIM, real_number> aSolid_minus = SolidBodyAcceleration(t - dt_2, params);
    Point<DIM, real_number> aSolid_plus = SolidBodyAcceleration(t + dt_2, params);

    // For each particle ...
    while (part.isNext())
    {
        // get particle a key
        vect_dist_key_dx a = part.get();

        if (vd.template getProp<type>(a) == BOUNDARY) // BOUNDARY
        {
            // move boundary particles with linear velocity and acceleration
            for (int xyz = 0; xyz < DIM; xyz++)
            {
                // trick to avoid accelerating the solid walls in the moving obstacle scenario
                if (vd.template getProp<vd_omega>(a) != 0.0)
                    vd.template getProp<force>(a)[xyz] = aSolid_minus.get(xyz);
                else
                    vd.template getProp<force>(a)[xyz] = 0.0;

                vd.template getProp<velocity>(a)[xyz] += dt_2 * vd.template getProp<force>(a)[xyz];
                vd.getPos(a)[xyz] += dt * vd.template getProp<velocity>(a)[xyz];

                // centre also moves
                vd.template getProp<force_transport>(a)[xyz] += dt * vd.template getProp<velocity>(a)[xyz];
            }

            // also rotate
            if (vd.template getProp<vd_omega>(a) != 0.0) // if rotating
            {
                real_number theta = vd.template getProp<vd_omega>(a) * dt;
                Point<DIM, real_number> normal = vd.template getProp<normal_vector>(a);

                // apply rotation
                Point<DIM, real_number> centre = vd.getProp<force_transport>(a); // stored in force_transport
                Point<DIM, real_number> x = vd.getPos(a);
                ApplyRotation(x, theta, centre);
                vd.getPos(a)[0] = x.get(0);
                vd.getPos(a)[1] = x.get(1);
                // update normals ( rotated wrt origin )
                ApplyRotation(normal, theta, {0.0, 0.0});
                vd.template getProp<normal_vector>(a)[0] = normal.get(0);
                vd.template getProp<normal_vector>(a)[1] = normal.get(1);
            }
        }
        else // FLUID
        {
            for (int xyz = 0; xyz < DIM; xyz++)
            {
                vd.template getProp<velocity>(a)[xyz] = vd.template getProp<velocity>(a)[xyz] + dt_2 * vd.template getProp<force>(a)[xyz];
                vd.template getProp<v_transport>(a)[xyz] = vd.template getProp<velocity>(a)[xyz] + dt_2 * vd.template getProp<force_transport>(a)[xyz];
                vd.getPos(a)[xyz] += dt * vd.template getProp<v_transport>(a)[xyz];
            }

            if (params.DENSITY_TYPE == DENSITY_DIFFERENTIAL)
                vd.template getProp<rho>(a) = vd.template getProp<rho>(a) + dt * vd.template getProp<drho>(a);
        }
        ++part;
    }
    // map particles if they have gone outside the domain
    vd.map();
    vd.ghost_get<type, rho, pressure, velocity, v_transport, normal_vector, curvature_boundary, arc_length, vd_volume, vd_omega>();

    // in density summation we need to update the density after moving the particles
    if (params.DENSITY_TYPE == DENSITY_SUMMATION)
    {
        CalcDensity(vd, NN, params);
        vd.ghost_get<rho>(KEEP_PROPERTIES);
    }

    // Calculate pressure from the density
    EqState(vd, params.rho_zero, params.B, params.gamma_, params.xi);
    vd.ghost_get<pressure>(KEEP_PROPERTIES);

    if (params.BC_TYPE == NO_SLIP)
    {
        ExtrapolateVelocity(vd, NN, params);
        vd.ghost_get<rho, pressure, v_transport>(KEEP_PROPERTIES);
    }

    obstalce_force_x = 0.0;
    obstacle_force_y = 0.0;
    calc_forces(vd, NN, obstalce_force_x, obstacle_force_y, calc_drag, params);
    if (calc_drag)
    {
        v_cl.sum(obstalce_force_x);
        v_cl.sum(obstacle_force_y);
        v_cl.execute();

        obstalce_force_x = obstalce_force_x * params.MassFluid;
        obstacle_force_y = obstacle_force_y * params.MassFluid;
    }

    // vd.ghost_get<type, rho, pressure, force, velocity, force_transport, v_transport, normal_vector, curvature_boundary, arc_length, vd_volume, vd_omega>();

    // particle iterator
    auto part2 = vd.getDomainIterator();

    // For each particle ...
    while (part2.isNext())
    {
        // get particle a key
        vect_dist_key_dx a = part2.get();

        if (vd.template getProp<type>(a) == BOUNDARY)
        {
            for (int xyz = 0; xyz < DIM; xyz++)
            {
                // trick to avoid accelerating the solid walls in the moving obstacle scenario
                if (vd.template getProp<vd_omega>(a) != 0.0)
                    vd.template getProp<force>(a)[xyz] = aSolid_plus.get(xyz);
                else
                    vd.template getProp<force>(a)[xyz] = 0.0;

                vd.template getProp<velocity>(a)[xyz] += dt_2 * vd.template getProp<force>(a)[xyz];
            }
        }
        else
        {
            for (int xyz = 0; xyz < DIM; xyz++)
            {
                vd.template getProp<velocity>(a)[xyz] = vd.template getProp<velocity>(a)[xyz] + dt_2 * vd.template getProp<force>(a)[xyz];
                vd.template getProp<v_transport>(a)[xyz] = vd.template getProp<velocity>(a)[xyz] + dt_2 * vd.template getProp<force_transport>(a)[xyz];
            }
        }

        ++part2;
    }

    // increment the iteration counter
    params.cnt++;
}

template <typename CellList>
void kick_drift_int_regularize(particles &vd,
                               CellList &NN,
                               const real_number dt,
                               Vcluster<> &v_cl,
                               size_t &niter,
                               Parameters &params)
{
    // particle iterator
    auto part = vd.getDomainIterator();

    const real_number dt_2 = dt * 0.5;
    // vd.ghost_get<type, rho, pressure, velocity, v_transport, normal_vector, curvature_boundary, arc_length, vd_volume, vd_omega>();

    // For each particle ...
    while (part.isNext())
    {
        // get particle a key
        vect_dist_key_dx a = part.get();

        if (vd.template getProp<type>(a) == FLUID)
        {
            for (int xyz = 0; xyz < DIM; xyz++)
            {
                vd.template getProp<v_transport>(a)[xyz] = vd.template getProp<velocity>(a)[xyz] + dt_2 * vd.template getProp<force_transport>(a)[xyz];
                vd.getPos(a)[xyz] += dt * vd.template getProp<v_transport>(a)[xyz];
            }
        }
        ++part;
    }

    // map particles if they have gone outside the domain
    vd.map();
    vd.ghost_get<type, rho, pressure, velocity, v_transport, normal_vector, curvature_boundary, arc_length, vd_volume, vd_omega>();

    // in density summation we need to update the density after moving the particles
    if (params.DENSITY_TYPE == DENSITY_SUMMATION)
    {
        CalcDensity(vd, NN, params);
        vd.ghost_get<rho>(KEEP_PROPERTIES);
    }

    // Calculate pressure from the density
    EqState(vd, params.rho_zero, params.B, params.gamma_, params.xi);
    vd.ghost_get<pressure>(KEEP_PROPERTIES);

    calc_forces_regularize(vd, NN, params);

    // particle iterator
    auto part2 = vd.getDomainIterator();

    // For each particle ...
    while (part2.isNext())
    {
        // get particle a key
        vect_dist_key_dx a = part2.get();

        if (vd.template getProp<type>(a) == FLUID)
        {
            for (int xyz = 0; xyz < DIM; xyz++)
            {
                vd.template getProp<v_transport>(a)[xyz] = vd.template getProp<velocity>(a)[xyz] + dt_2 * vd.template getProp<force_transport>(a)[xyz];
            }
        }

        ++part2;
    }

    // increment the iteration counter
    niter++;
}

#endif // TIMEINTEGRATION_HPP
