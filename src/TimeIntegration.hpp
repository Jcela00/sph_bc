#ifndef TIMEINTEGRATION_HPP
#define TIMEINTEGRATION_HPP

#include "Definitions.hpp"
#include "VectorUtilities.hpp"
#include "Physics.hpp"
#include "EqState.hpp"
#include "CalcDensity.hpp"
#include "ExtrapolateVelocity.hpp"
#include "CalcForces.hpp"

template <typename vd_type>
__global__ void max_velocity_gpu(vd_type vd)
{
    auto a = GET_PARTICLE(vd);

    Point<DIM, real_number> vel(vd.template getProp<vd4_velocity>(a));
    vd.template getProp<vd12_vel_red>(a) = norm(vel);
}

inline __device__ __host__ void max_velocity(particles &vd, real_number &max_vel)
{
    Vcluster<> &v_cl = create_vcluster();

    auto part = vd.getDomainIteratorGPU();

    CUDA_LAUNCH(max_velocity_gpu, part, vd.toKernel());

    max_vel = reduce_local<vd12_vel_red, _max_>(vd);

    v_cl.max(max_vel);
    v_cl.execute();
}

inline __device__ __host__ real_number calc_deltaT(particles &vd, const Parameters &params)
{
    real_number Maxvel = 0.0;
    max_velocity(vd, Maxvel);

    real_number dt_u = 0.25 * params.H / (params.cbar + abs(Maxvel));
    real_number dt_visc = 0.25 * params.H * params.H / (params.nu);
    real_number dt_g = 0.25 * sqrt(params.H / getVectorNorm(params.gravity_vector));
    real_number dt = params.CFLnumber * std::min({dt_u, dt_visc, dt_g});

    return dt;
}

template <typename vd_type>
__global__ void avg_vel_x_GPU(vd_type vd)
{
    auto a = GET_PARTICLE(vd);

    vd.template getProp<vd12_vel_red>(a) = vd.template getProp<vd4_velocity>(a)[0];
}

inline __device__ __host__ void avg_vel_x(particles &vd, real_number &avg_vel)
{
    Vcluster<> &v_cl = create_vcluster();

    auto part = vd.getDomainIteratorGPU();

    CUDA_LAUNCH(avg_vel_x_GPU, part, vd.toKernel());

    avg_vel = reduce_local<vd12_vel_red, _add_>(vd);

    v_cl.sum(avg_vel);
    v_cl.execute();
}

inline __device__ __host__ Point<3, real_number> ComputeDragLift(particles &vd, const Parameters &params)
{
    Vcluster<> &v_cl = create_vcluster();

    Point<3, real_number> VxDragLift{0.0};
    real_number vx = 0.0;
    avg_vel_x(vd, vx);

    real_number fx = reduce_local<vd13_force_red_x, _add_>(vd);
    real_number fy = reduce_local<vd14_force_red_y, _add_>(vd);

    v_cl.sum(fx);
    v_cl.sum(fy);
    v_cl.execute();

    VxDragLift[0] = vx;
    VxDragLift[1] = fx * params.MassFluid;
    VxDragLift[2] = fy * params.MassFluid;

    return VxDragLift;
}
// inline __device__ __host__ Point<DIM, real_number> SolidBodyAcceleration(real_number t, const Parameters &params)
// {
//     Point<DIM, real_number> a;
//     if (params.SCENARIO == MOVING_OBSTACLE)
//     {
//         real_number period = 10.0;
//         real_number amplitude = 0.25;
//         a.get(0) = amplitude * (4.0 * M_PI * M_PI / (period * period)) * sin(2.0 * M_PI * t / period);
//         a.get(1) = 0.0;
//     }
//     else
//     {
//         a = {0.0, 0.0};
//     }

//     return a;
// }

template <typename vd_type>
__global__ void kdi_1_gpu(vd_type vd, const real_number dt)
{
    const real_number dt_2 = dt * 0.5;

    // get particle a key
    auto a = GET_PARTICLE(vd);

    if (vd.template getProp<vd0_type>(a) == BOUNDARY) // BOUNDARY
    {
        // move boundary particles with their linear velocity
        for (int xyz = 0; xyz < DIM; xyz++)
        {
            vd.getPos(a)[xyz] += dt * vd.template getProp<vd4_velocity>(a)[xyz];
        }
    }
    else if (vd.template getProp<vd0_type>(a) == OBSTACLE)
    {
        if constexpr (DIM == 2)
        {
            // move obstacle particles with their linear velocity
            for (int xyz = 0; xyz < DIM; xyz++)
            {
                vd.getPos(a)[xyz] += dt * vd.template getProp<vd4_velocity>(a)[xyz];
                // centre also moves
                vd.template getProp<vd7_force_t>(a)[xyz] += dt * vd.template getProp<vd4_velocity>(a)[xyz];
            }

            // also rotate ( if omega == 0 this does nothing )
            real_number omega = vd.template getProp<vd10_omega>(a);
            real_number theta = omega * dt;
            Point<DIM, real_number> normal = vd.template getProp<vd8_normal>(a);

            // apply rotation
            Point<DIM, real_number> centre = vd.template getProp<vd7_force_t>(a); // stored in force_transport
            Point<DIM, real_number> x = vd.getPos(a);

            real_number xnorm = norm(x - centre);
            x = ApplyRotation(x, theta, centre);
            real_number xnorm_rot = norm(x - centre);
            // rotation should preserve the distance to the centre
            vd.getPos(a)[0] = centre.get(0) + (x.get(0) - centre.get(0)) * xnorm / xnorm_rot;
            vd.getPos(a)[1] = centre.get(1) + (x.get(1) - centre.get(1)) * xnorm / xnorm_rot;
            // update normals ( rotated wrt origin )
            normal = ApplyRotation(normal, theta, {0.0, 0.0});
            normalizeVector(normal); // renormalize to avoid numerical errors
            vd.template getProp<vd8_normal>(a)[0] = normal.get(0);
            vd.template getProp<vd8_normal>(a)[1] = normal.get(1);
        }
    }
    else if (vd.template getProp<vd0_type>(a) == FLUID)
    {
        // Otherwise it is a fluid particle
        for (int xyz = 0; xyz < DIM; xyz++)
        {
            vd.template getProp<vd4_velocity>(a)[xyz] = vd.template getProp<vd4_velocity>(a)[xyz] + dt_2 * vd.template getProp<vd6_force>(a)[xyz];
            vd.template getProp<vd5_velocity_t>(a)[xyz] = vd.template getProp<vd4_velocity>(a)[xyz] + dt_2 * vd.template getProp<vd7_force_t>(a)[xyz];
            vd.getPos(a)[xyz] += dt * vd.template getProp<vd5_velocity_t>(a)[xyz];
        }

        if (_params_gpu_.DENSITY_TYPE == DENSITY_DIFFERENTIAL)
            vd.template getProp<vd1_rho>(a) = vd.template getProp<vd1_rho>(a) + dt_2 * vd.template getProp<vd3_drho>(a);
    }

    return;
}

template <typename vd_type>
__global__ void kdi_2_gpu(vd_type vd, const real_number dt)
{
    const real_number dt_2 = dt * 0.5;

    // get particle a key
    auto a = GET_PARTICLE(vd);

    if (vd.template getProp<vd0_type>(a) == FLUID)
    {
        for (int xyz = 0; xyz < DIM; xyz++)
        {
            vd.template getProp<vd4_velocity>(a)[xyz] = vd.template getProp<vd4_velocity>(a)[xyz] + dt_2 * vd.template getProp<vd6_force>(a)[xyz];
            vd.template getProp<vd5_velocity_t>(a)[xyz] = vd.template getProp<vd4_velocity>(a)[xyz] + dt_2 * vd.template getProp<vd7_force_t>(a)[xyz];
        }
        if (_params_gpu_.DENSITY_TYPE == DENSITY_DIFFERENTIAL)
            vd.template getProp<vd1_rho>(a) = vd.template getProp<vd1_rho>(a) + dt_2 * vd.template getProp<vd3_drho>(a);
    }
}

template <typename CellList>
void kick_drift_int(particles &vd,
                    CellList &NN,
                    const real_number dt,
                    real_number t,
                    Parameters &params,
                    AuxiliarParameters &auxParams)
{

    // particle iterator
    auto it = vd.getDomainIteratorGPU();
    CUDA_LAUNCH(kdi_1_gpu, it, vd.toKernel(), dt);

    // map particles if they have gone outside the domain
    vd.map(RUN_ON_DEVICE);
    vd.ghost_get<vd0_type, vd1_rho, vd2_pressure, vd4_velocity, vd5_velocity_t, vd8_normal, vd9_volume, vd10_omega>(RUN_ON_DEVICE);

    // in density summation we need to update the density after moving the particles
    if (params.DENSITY_TYPE == DENSITY_SUMMATION)
    {
        CalcDensity(vd, NN, params);
        vd.ghost_get<vd1_rho>(KEEP_PROPERTIES | RUN_ON_DEVICE);
    }

    // Calculate pressure from the density
    EqState(vd, params.rho0, params.B, params.gamma, params.xi);
    vd.ghost_get<vd2_pressure>(KEEP_PROPERTIES | RUN_ON_DEVICE);

    if (params.BC_TYPE == NO_SLIP)
    {
        ExtrapolateVelocity(vd, NN);
        vd.ghost_get<vd1_rho, vd2_pressure, vd5_velocity_t>(KEEP_PROPERTIES | RUN_ON_DEVICE);
    }

    calc_forces(vd, NN, params);

    // particle iterator
    it = vd.getDomainIteratorGPU();

    CUDA_LAUNCH(kdi_2_gpu, it, vd.toKernel(), dt);

    if (params.DENSITY_TYPE == DENSITY_DIFFERENTIAL) // If density differential, density has been updated in kdi_2
    {
        // Calculate pressure from the density
        EqState(vd, params.rho0, params.B, params.gamma, params.xi);
        vd.ghost_get<vd2_pressure>(KEEP_PROPERTIES | RUN_ON_DEVICE);
    }

    // increment the iteration counter
    auxParams.cnt++;
}

#endif // TIMEINTEGRATION_HPP
