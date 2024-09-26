#ifndef TIMEINTEGRATION_HPP
#define TIMEINTEGRATION_HPP

#include "Definitions.hpp"
#include "VectorUtilities.hpp"
#include "Physics.hpp"
#include "EqState.hpp"
#include "CalcDensity.hpp"
#include "CalcForces.hpp"

// void max_velocity(particles &vd, Vcluster<> &v_cl, real_number &max_vel)
// {
//     // Calculate the maximum acceleration
//     auto part = vd.getDomainIterator();

//     while (part.isNext())
//     {
//         auto a = part.get();

//         Point<DIM, real_number> vel(vd.getProp<velocity>(a));
//         real_number vel2 = norm2(vel);

//         if (vel2 >= max_vel)
//             max_vel = vel2;

//         ++part;
//     }
//     max_vel = std::sqrt(max_vel);

//     // Vcluster<> &v_cl = create_vcluster();
//     v_cl.max(max_vel);
//     v_cl.execute();
// }

// inline __device__ __host__ real_number calc_deltaT(particles &vd, Vcluster<> &v_cl, const Parameters params)

inline __device__ __host__ real_number calc_deltaT(const Parameters params)
{
    real_number Maxvel = 0.0;
    // max_velocity(vd, v_cl, Maxvel);
    Maxvel = params.cbar / 10.0;

    real_number dt_u = 0.25 * params.H / (params.cbar + abs(Maxvel));
    real_number dt_visc = 0.25 * params.H * params.H / (params.nu);
    real_number dt_g = 0.25 * sqrt(params.H / getVectorNorm(params.gravity_vector));
    real_number dt = params.CFLnumber * std::min({dt_u, dt_visc, dt_g});
    // if (dt < DtMin)
    // 	dt = DtMin;

    return dt;
}

inline __device__ __host__ Point<DIM, real_number> SolidBodyAcceleration(real_number t, const Parameters &params)
{
    Point<DIM, real_number> a;
    if (params.SCENARIO == MOVING_OBSTACLE)
    {
        real_number period = 10.0;
        real_number amplitude = 0.25;
        a.get(0) = amplitude * (4.0 * M_PI * M_PI / (period * period)) * sin(2.0 * M_PI * t / period);
        a.get(1) = 0.0;
    }
    else
    {
        a = {0.0, 0.0};
    }

    return a;
}

template <typename vd_type>
__global__ void kdi_1_gpu(vd_type vd,
                          const real_number dt,
                          real_number t,
                          Parameters params)
{
    const real_number dt_2 = dt * 0.5;
    Point<DIM, real_number> aSolid_minus = 0.0; // SolidBodyAcceleration(t - dt_2, params);

    // get particle a key
    auto a = GET_PARTICLE(vd);

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
            Point<DIM, real_number> centre = vd.template getProp<force_transport>(a); // stored in force_transport
            Point<DIM, real_number> x = vd.getPos(a);
            // ApplyRotation(x, theta, centre);
            vd.getPos(a)[0] = x.get(0);
            vd.getPos(a)[1] = x.get(1);
            // update normals ( rotated wrt origin )
            // ApplyRotation(normal, theta, {0.0, 0.0});
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
}
template <typename vd_type>
__global__ void kdi_2_gpu(vd_type vd,
                          const real_number dt,
                          real_number t,
                          Parameters params)
{

    const real_number dt_2 = dt * 0.5;
    Point<DIM, real_number> aSolid_plus = 0.0; // SolidBodyAcceleration(t + dt_2, params);

    // get particle a key
    auto a = GET_PARTICLE(vd);

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
}

template <typename CellList>
void kick_drift_int(particles &vd,
                    CellList &NN,
                    const real_number dt,
                    real_number t,
                    Parameters &params)
{
    vd.ghost_get<type, rho, pressure, velocity, v_transport, normal_vector, curvature_boundary, arc_length, vd_volume, vd_omega, vd_vorticity>(RUN_ON_DEVICE);

    // particle iterator
    auto it = vd.getDomainIteratorGPU();

    CUDA_LAUNCH(kdi_1_gpu, it, vd.toKernel(), dt, t, params);
    cudaError_t err = cudaDeviceSynchronize(); // Wait for the kernel to finish
    if (err != cudaSuccess)
    {
        printf("CUDA error: %s\n", cudaGetErrorString(err));
    }

    // map particles if they have gone outside the domain
    vd.map(RUN_ON_DEVICE);

    vd.ghost_get<type, rho, pressure, velocity, v_transport, normal_vector, curvature_boundary, arc_length, vd_volume, vd_omega, vd_vorticity>(RUN_ON_DEVICE);

    // in density summation we need to update the density after moving the particles
    if (params.DENSITY_TYPE == DENSITY_SUMMATION)
    {
        CalcDensity(vd, NN, params);
        vd.ghost_get<type, rho, pressure, velocity, v_transport, normal_vector, curvature_boundary, arc_length, vd_volume, vd_omega>(RUN_ON_DEVICE);
    }

    // Calculate pressure from the density
    EqState(vd, params.rho_zero, params.B, params.gamma_, params.xi);
    vd.ghost_get<type, rho, pressure, velocity, v_transport, normal_vector, curvature_boundary, arc_length, vd_volume, vd_omega, vd_vorticity>(RUN_ON_DEVICE);

    // if (params.BC_TYPE == NO_SLIP)
    // {
    //     ExtrapolateVelocity(vd, NN, params);
    //     vd.ghost_get<rho, pressure, v_transport>(KEEP_PROPERTIES);
    // }

    // vd.deviceToHostPos();
    // vd.deviceToHostProp<type, rho, pressure, force, velocity, force_transport, v_transport, normal_vector, curvature_boundary, arc_length, vd_volume, vd_omega, vd_vorticity>();
    // vd.write_frame("BEFORE force", 0, params.WRITER);
    calc_forces(vd, NN, params);

    // vd.deviceToHostPos();
    // vd.deviceToHostProp<type, rho, pressure, force, velocity, force_transport, v_transport, normal_vector, curvature_boundary, arc_length, vd_volume, vd_omega, vd_vorticity>();
    // vd.write_frame("AFTER force", 0, params.WRITER);

    // if (calc_drag)
    // {
    //     v_cl.sum(obstalce_force_x);
    //     v_cl.sum(obstacle_force_y);
    //     v_cl.execute();

    //     obstalce_force_x = obstalce_force_x * params.MassFluid;
    //     obstacle_force_y = obstacle_force_y * params.MassFluid;
    // }

    // BEFORE GPU GHOST GET
    vd.ghost_get<type, rho, pressure, velocity, v_transport, normal_vector, curvature_boundary, arc_length, vd_volume, vd_omega, vd_vorticity>(RUN_ON_DEVICE);

    // particle iterator
    auto it2 = vd.getDomainIteratorGPU();

    CUDA_LAUNCH(kdi_2_gpu, it2, vd.toKernel(), dt, t, params);
    err = cudaDeviceSynchronize(); // Wait for the kernel to finish
    if (err != cudaSuccess)
    {
        printf("CUDA error: %s\n", cudaGetErrorString(err));
    }

    // vd.deviceToHostPos();
    // vd.deviceToHostProp<type, rho, pressure, force, velocity, force_transport, v_transport, normal_vector, curvature_boundary, arc_length, vd_volume, vd_omega, vd_vorticity>();
    // vd.write_frame("AFTER second kdi", 0, params.WRITER);
    // increment the iteration counter
    params.cnt++;
}

// template <typename CellList>
// void kick_drift_int_regularize(particles &vd,
//                                CellList &NN,
//                                const real_number dt,
//                                Vcluster<> &v_cl,
//                                size_t &niter,
//                                Parameters &params)
// {
//     // particle iterator
//     auto part = vd.getDomainIterator();

//     const real_number dt_2 = dt * 0.5;
//     // vd.ghost_get<type, rho, pressure, velocity, v_transport, normal_vector, curvature_boundary, arc_length, vd_volume, vd_omega>();

//     // For each particle ...
//     while (part.isNext())
//     {
//         // get particle a key
//         vect_dist_key_dx a = part.get();

//         if (vd.template getProp<type>(a) == FLUID)
//         {
//             for (int xyz = 0; xyz < DIM; xyz++)
//             {
//                 vd.template getProp<v_transport>(a)[xyz] = vd.template getProp<velocity>(a)[xyz] + dt_2 * vd.template getProp<force_transport>(a)[xyz];
//                 vd.getPos(a)[xyz] += dt * vd.template getProp<v_transport>(a)[xyz];
//             }
//         }
//         ++part;
//     }

//     // map particles if they have gone outside the domain
//     vd.map();
//     vd.ghost_get<type, rho, pressure, velocity, v_transport, normal_vector, curvature_boundary, arc_length, vd_volume, vd_omega>();

//     // in density summation we need to update the density after moving the particles
//     if (params.DENSITY_TYPE == DENSITY_SUMMATION)
//     {
//         CalcDensity(vd, NN, params);
//         vd.ghost_get<rho>(KEEP_PROPERTIES);
//     }

//     // Calculate pressure from the density
//     EqState(vd, params.rho_zero, params.B, params.gamma_, params.xi);
//     vd.ghost_get<pressure>(KEEP_PROPERTIES);

//     calc_forces_regularize(vd, NN, params);

//     // particle iterator
//     auto part2 = vd.getDomainIterator();

//     // For each particle ...
//     while (part2.isNext())
//     {
//         // get particle a key
//         vect_dist_key_dx a = part2.get();

//         if (vd.template getProp<type>(a) == FLUID)
//         {
//             for (int xyz = 0; xyz < DIM; xyz++)
//             {
//                 vd.template getProp<v_transport>(a)[xyz] = vd.template getProp<velocity>(a)[xyz] + dt_2 * vd.template getProp<force_transport>(a)[xyz];
//             }
//         }

//         ++part2;
//     }

//     // increment the iteration counter
//     niter++;
// }

#endif // TIMEINTEGRATION_HPP
