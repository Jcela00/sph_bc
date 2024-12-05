#ifndef TIMEINTEGRATION_HPP
#define TIMEINTEGRATION_HPP

#include "Definitions.hpp"
#include "VectorUtilities.hpp"
#include "Physics.hpp"
#include "EqState.hpp"
#include "CalcDensity.hpp"
#include "ExtrapolateVelocity.hpp"
#include "CalcForces.hpp"
#include "AssignNormals.hpp"

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

__device__ real_number SquareAcceleration(const real_number t)
{
    const real_number A = 2.82095118865606;
    const real_number mu = 0.525652150817692;
    const real_number sigma = 0.141421507636543;

    return A * std::exp(-0.5 * (t - mu) * (t - mu) / (sigma * sigma));
}

__device__ real_number VelocityRamp(const real_number t)
{
    const real_number timetosteady = 4.0;
    if (t <= timetosteady)
    {
        return (t / timetosteady) * (t / timetosteady) * (3.0 - 2.0 * (t / timetosteady));
    }
    else
    {
        return 1.0;
    }

    // return 1.0;
}

// template <typename vd_type, typename NN_type>
// __global__ void SensorInflow_gpu(vd_type vd, NN_type NN, real_number t)
// {
//     // Key of the particle a
//     unsigned int a;
//     GET_PARTICLE_SORT(a, NN);

//     // if particle FLUID
//     if (vd.template getProp<vd0_type>(a) == FLUID)
//     {
//         Point<DIM, real_number> xa = vd.getPos(a);

//         if (xa.get(0) < 10.0 * _params_gpu_.H && xa.get(0) > 8.5 * _params_gpu_.H)
//         {
//             real_number rho_sum = 0.0;
//             real_number p_sum = 0.0;
//             real_number u_sum = 0.0;
//             real_number w_sum = 0.0;

//             auto Np = NN.getNNIteratorBox(NN.getCell(xa));
//             while (Np.isNext() == true)
//             {
//                 auto b = Np.get_sort();

//                 const Point<DIM, real_number> xb = vd.getPos(b);

//                 if (xb.get(0) > 10.0 * _params_gpu_.H && vd.template getProp<vd0_type>(b) == FLUID)
//                 {
//                     const Point<DIM, real_number> dr = xa - xb;
//                     const real_number r2 = norm2(dr);

//                     if (r2 < _params_gpu_.r_cut2)
//                     {
//                         const real_number r = sqrt(r2);

//                         const real_number w = Wab(r, _params_gpu_.H, _params_gpu_.Kquintic);

//                         w_sum += w;
//                         rho_sum += vd.template getProp<vd1_rho>(b);
//                         p_sum += vd.template getProp<vd2_pressure>(b);
//                         u_sum += vd.template getProp<vd4_velocity>(b)[0];
//                     }
//                 }
//                 ++Np;
//             }

//             if (w_sum != 0.0)
//             {
//                 vd.template getProp<vd1_rho>(a) = rho_sum / w_sum;
//                 vd.template getProp<vd2_pressure>(a) = p_sum / w_sum;
//                 vd.template getProp<vd4_velocity>(a)[0] = u_sum / w_sum;
//             }
//             else
//             {
//                 vd.template getProp<vd1_rho>(a) = _params_gpu_.rho0;
//                 vd.template getProp<vd2_pressure>(a) = 0.0;
//                 const real_number fac = VelocityRamp(t);
//                 vd.template getProp<vd4_velocity>(a)[0] = fac * _params_gpu_.Vinflow[0];
//             }
//         }
//     }
// }

// template <typename vd_type, typename NN_type>
// inline void SensorInflow(vd_type &vd, NN_type &NN, const Parameters &params, const real_number t)
// {
//     // This function sensors velocity and pressure from domain interior

//     auto it = vd.getDomainIteratorGPU();

//     vd.template updateCellListGPU<vd0_type, vd1_rho, vd2_pressure, vd4_velocity>(NN);
//     CUDA_LAUNCH(SensorInflow_gpu, it, vd.toKernel(), NN.toKernel(), t);
//     vd.template restoreOrder<vd0_type, vd1_rho, vd2_pressure, vd4_velocity>(NN);
// }

// template <typename vd_type>
// __global__ void kdi_1_gpu(vd_type vd, const real_number dt, const real_number t)
// {
//     const real_number dt_2 = dt * 0.5;

//     // get particle a key
//     auto a = GET_PARTICLE(vd);

//     if (vd.template getProp<vd0_type>(a) == BOUNDARY) // BOUNDARY
//     {
//         // move boundary particles with their linear velocity
//         for (int xyz = 0; xyz < DIM; xyz++)
//         {
//             vd.getPos(a)[xyz] += dt * vd.template getProp<vd4_velocity>(a)[xyz];
//         }
//     }
//     else if (vd.template getProp<vd0_type>(a) == OBSTACLE)
//     {
//         if constexpr (DIM == 2)
//         {
//             // move obstacle particles with their linear velocity
//             for (int xyz = 0; xyz < DIM; xyz++)
//             {
//                 if (_params_gpu_.SCENARIO == MOVING_OBSTACLE)
//                 {
//                     bool is_x = (xyz == 0);
//                     vd.template getProp<vd6_force>(a)[xyz] = SquareAcceleration(t) * is_x;
//                     vd.template getProp<vd4_velocity>(a)[xyz] += dt * vd.template getProp<vd6_force>(a)[xyz];
//                 }

//                 vd.getPos(a)[xyz] += dt * vd.template getProp<vd4_velocity>(a)[xyz];
//                 // centre also moves
//                 vd.template getProp<vd7_force_t>(a)[xyz] += dt * vd.template getProp<vd4_velocity>(a)[xyz];
//             }

//             // also rotate ( if omega == 0 this does nothing )
//             real_number omega = vd.template getProp<vd10_omega>(a);
//             real_number theta = omega * dt;
//             Point<DIM, real_number> normal = vd.template getProp<vd8_normal>(a);

//             // apply rotation
//             Point<DIM, real_number> centre = vd.template getProp<vd7_force_t>(a); // stored in force_transport
//             Point<DIM, real_number> x = vd.getPos(a);

//             real_number xnorm = norm(x - centre);
//             x = ApplyRotation(x, theta, centre);
//             real_number xnorm_rot = norm(x - centre);
//             // rotation should preserve the distance to the centre
//             vd.getPos(a)[0] = centre.get(0) + (x.get(0) - centre.get(0)) * xnorm / xnorm_rot;
//             vd.getPos(a)[1] = centre.get(1) + (x.get(1) - centre.get(1)) * xnorm / xnorm_rot;
//             // update normals ( rotated wrt origin )
//             normal = ApplyRotation(normal, theta, {0.0, 0.0});
//             normalizeVector(normal); // renormalize to avoid numerical errors
//             vd.template getProp<vd8_normal>(a)[0] = normal.get(0);
//             vd.template getProp<vd8_normal>(a)[1] = normal.get(1);
//         }
//     }
//     else if (vd.template getProp<vd0_type>(a) == FLUID) // Otherwise it is a fluid particle
//     {

//         for (int xyz = 0; xyz < DIM; xyz++)
//         {
//             vd.template getProp<vd4_velocity>(a)[xyz] = vd.template getProp<vd4_velocity>(a)[xyz] + dt_2 * vd.template getProp<vd6_force>(a)[xyz];
//             vd.template getProp<vd5_velocity_t>(a)[xyz] = vd.template getProp<vd4_velocity>(a)[xyz] + dt_2 * vd.template getProp<vd7_force_t>(a)[xyz];
//             vd.getPos(a)[xyz] += dt * vd.template getProp<vd5_velocity_t>(a)[xyz];
//         }

//         if (_params_gpu_.DENSITY_TYPE == DENSITY_DIFFERENTIAL)
//             vd.template getProp<vd1_rho>(a) = vd.template getProp<vd1_rho>(a) + dt_2 * vd.template getProp<vd3_drho>(a);

//         if (_params_gpu_.SCENARIO == CUSTOM)
//         {
//             Point<DIM, real_number> pos = vd.getPos(a);
//             if (pos.get(0) < 10.0 * _params_gpu_.H)
//             {
//                 if (pos.get(0) > 8.5 * _params_gpu_.H)
//                 {
//                     const real_number fac = VelocityRamp(t);

//                     const real_number a_speed = _params_gpu_.cbar;
//                     const real_number rhoref = _params_gpu_.rho0;
//                     const real_number pref = 0;
//                     const real_number uref = fac * _params_gpu_.Vinflow[0];

//                     const real_number rho = vd.template getProp<vd1_rho>(a);
//                     const real_number p = vd.template getProp<vd2_pressure>(a);
//                     const real_number u = vd.template getProp<vd4_velocity>(a)[0];

//                     const real_number J1 = 0.0;
//                     const real_number J2 = 0.0;
//                     const real_number J3 = -rho * a_speed * (u - uref) + (p - pref);

//                     // vd.template getProp<vd1_rho>(a) = _params_gpu_.rho0;
//                     // vd.template getProp<vd2_pressure>(a) = 0.0;
//                     vd.template getProp<vd1_rho>(a) = rhoref + (1.0 / (a_speed * a_speed)) * (-J1 + 0.5 * J2 + 0.5 * J3);
//                     vd.template getProp<vd2_pressure>(a) = pref + 0.5 * (J2 + J3);
//                     vd.template getProp<vd4_velocity>(a)[0] = uref + (1.0 / (2 * rho * a_speed)) * (J2 - J3);
//                     vd.template getProp<vd4_velocity>(a)[1] = fac * _params_gpu_.Vinflow[1];
//                     if constexpr (DIM == 3)
//                         vd.template getProp<vd4_velocity>(a)[2] = fac * _params_gpu_.Vinflow[2];
//                 }
//                 else
//                 {
//                     const real_number fac = VelocityRamp(t);

//                     const real_number rhoref = _params_gpu_.rho0;
//                     const real_number pref = 0;
//                     const real_number uref = fac * _params_gpu_.Vinflow[0];

//                     vd.template getProp<vd1_rho>(a) = rhoref;
//                     vd.template getProp<vd2_pressure>(a) = pref;
//                     vd.template getProp<vd4_velocity>(a)[0] = uref;
//                     vd.template getProp<vd4_velocity>(a)[1] = fac * _params_gpu_.Vinflow[1];
//                     if constexpr (DIM == 3)
//                         vd.template getProp<vd4_velocity>(a)[2] = fac * _params_gpu_.Vinflow[2];
//                 }
//             }
//             else if (pos.get(0) > _params_gpu_.length[0] - 10.0 * _params_gpu_.H)
//             {
//                 const real_number fac = VelocityRamp(t);

//                 const real_number a_speed = _params_gpu_.cbar;
//                 const real_number rhoref = _params_gpu_.rho0;
//                 const real_number pref = 0;

//                 const real_number rho = vd.template getProp<vd1_rho>(a);
//                 const real_number p = vd.template getProp<vd2_pressure>(a);
//                 const real_number u = vd.template getProp<vd4_velocity>(a)[0];
//                 const real_number uref = u;

//                 const real_number J1 = -a_speed * a_speed * (rho - rhoref) + (p - pref);
//                 const real_number J2 = rho * a_speed * (u - uref) + (p - pref);
//                 const real_number J3 = 0.0;

//                 vd.template getProp<vd1_rho>(a) = rhoref + (1.0 / (a_speed * a_speed)) * (-J1 + 0.5 * J2 + 0.5 * J3);
//                 vd.template getProp<vd2_pressure>(a) = pref + 0.5 * (J2 + J3);
//                 vd.template getProp<vd4_velocity>(a)[0] = uref + (1.0 / (2 * rho * a_speed)) * (J2 - J3);
//             }
//         }
//     }
// }

// template <typename vd_type>
// __global__ void kdi_2_gpu(vd_type vd, const real_number dt, const real_number t)
// {
//     const real_number dt_2 = dt * 0.5;

//     // get particle a key
//     auto a = GET_PARTICLE(vd);

//     if (vd.template getProp<vd0_type>(a) == FLUID)
//     {
//         for (int xyz = 0; xyz < DIM; xyz++)
//         {
//             vd.template getProp<vd4_velocity>(a)[xyz] = vd.template getProp<vd4_velocity>(a)[xyz] + dt_2 * vd.template getProp<vd6_force>(a)[xyz];
//             vd.template getProp<vd5_velocity_t>(a)[xyz] = vd.template getProp<vd4_velocity>(a)[xyz] + dt_2 * vd.template getProp<vd7_force_t>(a)[xyz];
//         }
//         if (_params_gpu_.DENSITY_TYPE == DENSITY_DIFFERENTIAL)
//             vd.template getProp<vd1_rho>(a) = vd.template getProp<vd1_rho>(a) + dt_2 * vd.template getProp<vd3_drho>(a);

//         if (_params_gpu_.SCENARIO == CUSTOM)
//         {
//             Point<DIM, real_number> pos = vd.getPos(a);
//             if (pos.get(0) < 10.0 * _params_gpu_.H)
//             {
//                 // if (pos.get(0) > 8.5 * _params_gpu_.H)
//                 // {
//                 //     const real_number fac = VelocityRamp(t);

//                 //     const real_number a_speed = _params_gpu_.cbar;
//                 //     const real_number rhoref = _params_gpu_.rho0;
//                 //     const real_number pref = 0;
//                 //     const real_number uref = fac * _params_gpu_.Vinflow[0];

//                 //     const real_number rho = vd.template getProp<vd1_rho>(a);
//                 //     const real_number p = vd.template getProp<vd2_pressure>(a);
//                 //     const real_number u = vd.template getProp<vd4_velocity>(a)[0];

//                 //     const real_number J1 = 0.0;
//                 //     const real_number J2 = 0.0;
//                 //     const real_number J3 = -rho * a_speed * (u - uref) + (p - pref);

//                 //     // vd.template getProp<vd1_rho>(a) = _params_gpu_.rho0;
//                 //     // vd.template getProp<vd2_pressure>(a) = 0.0;
//                 //     vd.template getProp<vd1_rho>(a) = rhoref + (1.0 / (a_speed * a_speed)) * (-J1 + 0.5 * J2 + 0.5 * J3);
//                 //     vd.template getProp<vd2_pressure>(a) = pref + 0.5 * (J2 + J3);
//                 //     vd.template getProp<vd4_velocity>(a)[0] = uref + (1.0 / (2 * rho * a_speed)) * (J2 - J3);
//                 //     vd.template getProp<vd4_velocity>(a)[1] = fac * _params_gpu_.Vinflow[1];
//                 //     if constexpr (DIM == 3)
//                 //         vd.template getProp<vd4_velocity>(a)[2] = fac * _params_gpu_.Vinflow[2];
//                 // }
//                 // else
//                 // {
//                 const real_number fac = VelocityRamp(t);

//                 const real_number rhoref = _params_gpu_.rho0;
//                 const real_number pref = 0;
//                 const real_number uref = fac * _params_gpu_.Vinflow[0];

//                 vd.template getProp<vd1_rho>(a) = rhoref;
//                 vd.template getProp<vd2_pressure>(a) = pref;
//                 vd.template getProp<vd4_velocity>(a)[0] = uref;
//                 vd.template getProp<vd4_velocity>(a)[1] = fac * _params_gpu_.Vinflow[1];
//                 if constexpr (DIM == 3)
//                     vd.template getProp<vd4_velocity>(a)[2] = fac * _params_gpu_.Vinflow[2];
//                 // }
//             }
//             else if (pos.get(0) > _params_gpu_.length[0] - 10.0 * _params_gpu_.H)
//             {
//                 const real_number fac = VelocityRamp(t);

//                 const real_number a_speed = _params_gpu_.cbar;
//                 const real_number rhoref = _params_gpu_.rho0;
//                 const real_number pref = 0;

//                 const real_number rho = vd.template getProp<vd1_rho>(a);
//                 const real_number p = vd.template getProp<vd2_pressure>(a);
//                 const real_number u = vd.template getProp<vd4_velocity>(a)[0];
//                 const real_number uref = u;

//                 const real_number J1 = -a_speed * a_speed * (rho - rhoref) + (p - pref);
//                 const real_number J2 = rho * a_speed * (u - uref) + (p - pref);
//                 const real_number J3 = 0.0;

//                 vd.template getProp<vd1_rho>(a) = rhoref + (1.0 / (a_speed * a_speed)) * (-J1 + 0.5 * J2 + 0.5 * J3);
//                 vd.template getProp<vd2_pressure>(a) = pref + 0.5 * (J2 + J3);
//                 vd.template getProp<vd4_velocity>(a)[0] = uref + (1.0 / (2 * rho * a_speed)) * (J2 - J3);
//             }
//         }
//     }
// }

template <typename vd_type>
__global__ void kdi_1_gpu(vd_type vd, const real_number dt, const real_number t)
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
                if (_params_gpu_.SCENARIO == MOVING_OBSTACLE)
                {
                    bool is_x = (xyz == 0);
                    vd.template getProp<vd6_force>(a)[xyz] = SquareAcceleration(t) * is_x;
                    vd.template getProp<vd4_velocity>(a)[xyz] += dt * vd.template getProp<vd6_force>(a)[xyz];
                }

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
    else if (vd.template getProp<vd0_type>(a) == FLUID) // Otherwise it is a fluid particle
    {
        if (_params_gpu_.SCENARIO == CUSTOM)
        {
            Point<DIM, real_number> pos = vd.getPos(a);
            if (pos.get(0) < 5.0 * _params_gpu_.H)
            {
                const real_number fac = VelocityRamp(t);
                // vd.template getProp<vd1_rho>(a) = _params_gpu_.rho0;
                for (int xyz = 0; xyz < DIM; xyz++)
                {
                    vd.template getProp<vd4_velocity>(a)[xyz] = fac * _params_gpu_.Vinflow[xyz];
                }
            }
            // else if (pos.get(0) > _params_gpu_.length[0] - 10.0 * _params_gpu_.H)
            // {
            // }
        }

        for (int xyz = 0; xyz < DIM; xyz++)
        {
            vd.template getProp<vd4_velocity>(a)[xyz] = vd.template getProp<vd4_velocity>(a)[xyz] + dt_2 * vd.template getProp<vd6_force>(a)[xyz];
            vd.template getProp<vd5_velocity_t>(a)[xyz] = vd.template getProp<vd4_velocity>(a)[xyz] + dt_2 * vd.template getProp<vd7_force_t>(a)[xyz];
            vd.getPos(a)[xyz] += dt * vd.template getProp<vd5_velocity_t>(a)[xyz];
        }

        if (_params_gpu_.DENSITY_TYPE == DENSITY_DIFFERENTIAL)
            vd.template getProp<vd1_rho>(a) = vd.template getProp<vd1_rho>(a) + dt_2 * vd.template getProp<vd3_drho>(a);
    }
}

template <typename vd_type>
__global__ void kdi_2_gpu(vd_type vd, const real_number dt, const real_number t)
{
    const real_number dt_2 = dt * 0.5;

    // get particle a key
    auto a = GET_PARTICLE(vd);

    if (vd.template getProp<vd0_type>(a) == FLUID)
    {
        if (_params_gpu_.SCENARIO == CUSTOM)
        {
            Point<DIM, real_number> pos = vd.getPos(a);

            if (pos.get(0) < 5.0 * _params_gpu_.H)
            {
                // const real_number fac = 1.0;
                const real_number fac = VelocityRamp(t);
                // vd.template getProp<vd1_rho>(a) = _params_gpu_.rho0;

                for (int xyz = 0; xyz < DIM; xyz++)
                {
                    vd.template getProp<vd4_velocity>(a)[xyz] = fac * _params_gpu_.Vinflow[xyz];
                }
            }
            // else if (pos.get(0) > _params_gpu_.length[0] - 10.0 * _params_gpu_.H)
            // {

            // }
        }
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

    // Particle iterator
    auto it = vd.getDomainIteratorGPU();

    // Kick drift step
    CUDA_LAUNCH(kdi_1_gpu, it, vd.toKernel(), dt, t);

    // map particles if they have gone outside the domain
    vd.map(RUN_ON_DEVICE);
    vd.ghost_get<vd0_type, vd1_rho, vd2_pressure, vd4_velocity, vd5_velocity_t, vd8_normal, vd9_volume, vd10_omega>(RUN_ON_DEVICE);

    // In new bc assign normals to fluid particles near to boundaries
    if (params.BC_TYPE == NEW_NO_SLIP)
    {
        AssignNormals(vd, NN);
        vd.ghost_get<vd8_normal>(KEEP_PROPERTIES | RUN_ON_DEVICE);
    }

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

    // Compute forces
    calc_forces(vd, NN, params);

    // particle iterator
    it = vd.getDomainIteratorGPU();

    // Second Kick step
    CUDA_LAUNCH(kdi_2_gpu, it, vd.toKernel(), dt, t);

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
