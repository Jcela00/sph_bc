#ifndef CALC_FORCES_HPP
#define CALC_FORCES_HPP

#include "Definitions.hpp"
#include "VectorUtilities.hpp"
#include "Kernel.hpp"
#include "Physics.hpp"
#include <cstdio>
#include <stdio.h>
#include <unistd.h>

template <typename vd_ker_type>
inline __device__ __host__ void interact_fluid_boundary_new(vd_ker_type &vd,
                                                            unsigned int &fluid_key,
                                                            const real_number massf,
                                                            const real_number rhof,
                                                            const real_number Pf,
                                                            const Point<DIM, real_number> xf,
                                                            const Point<DIM, real_number> vf,
                                                            const Point<DIM, real_number> xw,
                                                            const Point<DIM, real_number> r_wall_to_fluid,
                                                            unsigned long boundary_key,
                                                            const Parameters params)
{

    real_number dp = params.dp;
    // Points from fluid to wall
    Point<DIM, real_number> r_fluid_to_wall = -1.0 * r_wall_to_fluid;
    // real_number dist2marker = getVectorNorm(r_fluid_to_wall);
    Point<DIM, real_number> vw = vd.template getProp<velocity>(boundary_key);

    real_number ang_vel = vd.template getProp<vd_omega>(boundary_key);

    if (ang_vel != 0.0) // if the solid is rotating, we need to add the tangential velocity of the rotation to vw
    {
        // marker particles store centre of solid body in force_transport since it is unused
        // vector pointing from centre of rotation to marker particle
        const Point<DIM, real_number> radial_vec = {xw.get(0) - vd.template getProp<force_transport>(boundary_key)[0],
                                                    xw.get(1) - vd.template getProp<force_transport>(boundary_key)[1]};
        const real_number radius = getVectorNorm(radial_vec);
        // get vector tangential to the radial vector, rotation velocity is in this direction
        const Point<DIM, real_number> tangential_rotation = getPerpendicularUnit2D(radial_vec);

        // Wall velocity is linear velocity + w*R*tangential
        vw.get(0) += radius * ang_vel * tangential_rotation.get(0);
        vw.get(1) += radius * ang_vel * tangential_rotation.get(1);
    }

    // Get normal and tangential vectors for velocity mirroring
    const Point<DIM, real_number> normal = vd.template getProp<normal_vector>(boundary_key);
    const Point<DIM, real_number> tangential = getPerpendicularUnit2D(normal);

    // wall acceleration
    const Point<DIM, real_number> aw = vd.template getProp<force>(boundary_key);

    // Difference between fluid transport and momentum velocity
    const Point<DIM, real_number> vtf = vd.template getProp<v_transport>(fluid_key);
    const Point<DIM, real_number> vdiff_f = vtf - vf;

    // Project vf and vw on tangential and normal directions
    real_number vt = dotProduct(vf, tangential);
    real_number vn = dotProduct(vf, normal);
    real_number vwt = dotProduct(vw, tangential);

    // vertical distance from fluid particle to wall
    real_number lf = dotProduct(r_fluid_to_wall, normal);
    lf = (lf < 0.0 ? -1.0 * lf : lf); // absolute value

    // Get array of vectors from fluid to 3 boundary particles and its norms
    std::array<Point<DIM, real_number>, 3> r_boundary = getBoundaryPositions(r_fluid_to_wall, normal, dp);
    std::array<real_number, 3> r_boundary_norm = {getVectorNorm(r_boundary[0]), getVectorNorm(r_boundary[1]), getVectorNorm(r_boundary[2])};

    // const real_number dist2third = r_boundary_norm[2];
    // distance from 3 boundary particles to marker
    std::array<real_number, 3> lwall = {0.5f * dp, 1.5f * dp, 2.5f * dp};

    // project forces on normal direction
    real_number g_normal = dotProduct(params.gravity_vector, normal);
    real_number a_normal = dotProduct(aw, normal);

    // For gradient of Af tensor
    const std::array<Point<DIM, real_number>, DIM> Af = dyadicProduct(rhof * vf, vdiff_f);
    // to avoid division by zero
    lf = std::max(lf, 0.25f * dp);

    const Point<3, real_number> vol = vd.template getProp<vd_volume>(boundary_key);

    for (int i = 0; i < 3; i++) // for the 3 boundary particles
    {
        // real_number rmax = sqrt(3.0 * 3.0 - (0.5 + (real_number)i) * (0.5 + (real_number)i)) * dp;
        // real_number rmin = (3.0 - (0.5 + (real_number)i)) * dp;
        // real_number kappa_max = 1.0 / (3.0 * dp);

        // kappa 0 gets rmax, kappa = kappa_max gets rmin
        // real_number r_interp = (rmin - rmax) / kappa_max * kappa + rmax;
        // real_number r_interp = rmin;

        // if (dist2marker < r_interp)
        // {
        const real_number Mass_boundary = vol[i] * params.rho_zero;
        const Point<DIM, real_number> v_boundary = ((vwt - vt) * (lwall[i] / lf) + vwt) * tangential + vn * normal;
        const real_number p_boundary = Pf - rhof * (g_normal - a_normal) * (lf + lwall[i]); //  dot(r_boundary,normal) = -(lf+lw)
        const real_number rho_boundary = InvEqState_particle(p_boundary, params.rho_zero, params.B, params.gamma_, params.xi);

        // flip sign of r_boundary to get vector pointing from boundary to fluid (Force routines use the vector pointing from b to a)
        r_boundary[i] = -1.0 * r_boundary[i];

        // Evaluate kernel gradient
        const Point<DIM, real_number> DW = DWab(r_boundary[i], r_boundary_norm[i], params.H, params.Kquintic);

        // Compute forces
        const Point<DIM, real_number> v_rel = vf - v_boundary;
        const real_number Vb = (Mass_boundary / rho_boundary); // vb is mass/rho instead of directly vol[i] because it allows to variate with density

        const Point<DIM, real_number> ViscosityTerm = Pi_physical(r_boundary[i], r_boundary_norm[i], v_rel, DW, params.eta);
        const real_number PressureTerm = PressureForce(rhof, rho_boundary, Pf, p_boundary); //-p_boundary - Pf;
        const Point<DIM, real_number> DivATerm = 0.5 * matVec(Af, DW);
        // printf("Hello from interact_fluid_boundary_new\n");
        for (int xyz = 0; xyz < DIM; ++xyz)
        {
            // write to particles
            vd.template getProp<force>(fluid_key)[xyz] += 2.0 * (Vb / rhof) * (PressureTerm * DW.get(xyz) + ViscosityTerm.get(xyz) + DivATerm.get(xyz));
            vd.template getProp<force_transport>(fluid_key)[xyz] += -2.0 * (Vb / rhof) * (params.Pbackground) * DW.get(xyz);
        }
        // if (accumulate_force) // we accumulate x & y force on cylinder ( this is just to compute drag coefficient)
        // {
        //     obstacle_force_x += -2.0 * (Vb / rhof) * (PressureTerm * DW.get(0) + ViscosityTerm.get(0) + DivATerm.get(0));
        //     obstacle_force_y += -2.0 * (Vb / rhof) * (PressureTerm * DW.get(1) + ViscosityTerm.get(1) + DivATerm.get(1));
        // }

        if (params.DENSITY_TYPE == DENSITY_DIFFERENTIAL) // this doesnt work well I havent touched in a long time and I have made may changes
        {
            vd.template getProp<drho>(fluid_key) += Mass_boundary * dotProduct(v_rel, DW);
        }
        // }
    }
}

template <typename vd_ker_type>
inline __device__ __host__ void interact_fluid_boundary_old(vd_ker_type &vd,
                                                            unsigned int &fluid_key,
                                                            const real_number massf,
                                                            const real_number rhof,
                                                            const real_number Pf,
                                                            const Point<DIM, real_number> xf,
                                                            const Point<DIM, real_number> vf,
                                                            const Point<DIM, real_number> xb,
                                                            const Point<DIM, real_number> r_ab,
                                                            const real_number r2,
                                                            unsigned long boundary_key,
                                                            const Parameters params)
{
    const real_number massb = params.MassBound;
    const real_number rhob = vd.template getProp<rho>(boundary_key);
    const real_number Pb = vd.template getProp<pressure>(boundary_key);
    const Point<DIM, real_number> vb = vd.template getProp<velocity>(boundary_key);
    const Point<DIM, real_number> vb_noslip = vd.template getProp<v_transport>(boundary_key); // here we store the extrapolated velocity for no slip BC

    const real_number r = sqrt(r2);

    const Point<DIM, real_number> v_rel = vf - vb;
    const Point<DIM, real_number> v_rel_aux = vf - vb_noslip;

    const Point<DIM, real_number> DW = DWab(r_ab, r, params.H, params.Kquintic);

    const real_number Va2 = (massf / rhof) * (massf / rhof);
    const real_number Vb2 = (massb / rhob) * (massb / rhob);

    const Point<DIM, real_number> ViscosityTerm = Pi_physical(r_ab, r, v_rel_aux, DW, params.eta);
    const real_number PressureTerm = PressureForce(rhof, rhob, Pf, Pb);

    const Point<DIM, real_number> vtf = vd.template getProp<v_transport>(fluid_key);
    const Point<DIM, real_number> vdiff_f = vtf - vf;
    // Boundary particles have no transport velocity difference
    std::array<Point<DIM, real_number>, DIM> Af = dyadicProduct(rhof * vf, vdiff_f);
    Point<DIM, real_number> DivATerm = 0.5 * matVec(Af, DW);

    for (int xyz = 0; xyz < DIM; ++xyz)
    {
        vd.template getProp<force>(fluid_key)[xyz] += (Va2 + Vb2) * (PressureTerm * DW.get(xyz) + ViscosityTerm.get(xyz) + DivATerm.get(xyz)) / massf;
        vd.template getProp<force_transport>(fluid_key)[xyz] += -1.0 * (Va2 + Vb2) * params.Pbackground * DW.get(xyz) / massf;
    }
    // if (accumulate_force) // we accumulate x and y force on obstacle for drag and lift coefficient
    // {
    //     obstacle_force_x += -(Va2 + Vb2) * (PressureTerm * DW.get(0) + ViscosityTerm.get(0) + DivATerm.get(0)) / massf; //+ 1.0 * (Va2 + Vb2) * params.Pbackground * DW.get(0) / massf;
    //     obstacle_force_y += -(Va2 + Vb2) * (PressureTerm * DW.get(1) + ViscosityTerm.get(1) + DivATerm.get(1)) / massf; //+ 1.0 * (Va2 + Vb2) * params.Pbackground * DW.get(1) / massf;
    //     // obstacle_force_y += sqrt(std::pow((Va2 + Vb2) * (PressureTerm * DW.get(0) + ViscosityTerm.get(0) + DivATerm.get(0)) / massf, 2) + std::pow((Va2 + Vb2) * (PressureTerm * DW.get(1) + ViscosityTerm.get(1) + DivATerm.get(1)) / massf, 2));
    // }

    if (params.DENSITY_TYPE == DENSITY_DIFFERENTIAL)
    {
        vd.template getProp<drho>(fluid_key) += massb * dotProduct(v_rel, DW);
    }
}

template <typename vd_ker_type>
inline __device__ __host__ void interact_fluid_fluid(vd_ker_type &vd,
                                                     unsigned int &a_key,
                                                     const real_number massa,
                                                     const real_number rhoa,
                                                     const real_number Pa,
                                                     const Point<DIM, real_number> xa,
                                                     const Point<DIM, real_number> va,
                                                     const Point<DIM, real_number> xb,
                                                     const Point<DIM, real_number> r_ab,
                                                     const real_number r2,
                                                     const unsigned long b_key,
                                                     const Parameters params)
{

    const real_number massb = params.MassFluid;
    const real_number rhob = vd.template getProp<rho>(b_key);
    const real_number Pb = vd.template getProp<pressure>(b_key);
    const Point<DIM, real_number> vb = vd.template getProp<velocity>(b_key);

    const real_number r = sqrt(r2);

    const Point<DIM, real_number> v_rel = va - vb;

    const Point<DIM, real_number> DW = DWab(r_ab, r, params.H, params.Kquintic);

    const real_number Va2 = (massa / rhoa) * (massa / rhoa);
    const real_number Vb2 = (massb / rhob) * (massb / rhob);

    const Point<DIM, real_number> ViscosityTerm = Pi_physical(r_ab, r, v_rel, DW, params.eta);
    const real_number PressureTerm = PressureForce(rhoa, rhob, Pa, Pb);

    const Point<DIM, real_number> vta = vd.template getProp<v_transport>(a_key);
    const Point<DIM, real_number> vtb = vd.template getProp<v_transport>(b_key);
    const Point<DIM, real_number> vdiff_a = vta - va;
    const Point<DIM, real_number> vdiff_b = vtb - vb;

    // real_number DivATerm = 0.5 * (rhoa * dotProduct(va, vdiff_a) + rhob * dotProduct(vb, vdiff_b));
    std::array<Point<DIM, real_number>, DIM> Aa = dyadicProduct(rhoa * va, vdiff_a);
    std::array<Point<DIM, real_number>, DIM> Ab = dyadicProduct(rhob * vb, vdiff_b);
    std::array<Point<DIM, real_number>, DIM> SumA;
    SumA[0] = Aa[0] + Ab[0];
    SumA[1] = Aa[1] + Ab[1];
    if constexpr (DIM == 3)
        SumA[2] = Aa[2] + Ab[2];

    Point<DIM, real_number> DivATerm = 0.5 * matVec(SumA, DW);

    for (int xyz = 0; xyz < DIM; ++xyz)
    {
        vd.template getProp<force>(a_key)[xyz] += (Va2 + Vb2) * (PressureTerm * DW.get(xyz) + ViscosityTerm.get(xyz) + DivATerm.get(xyz)) / massa;
        vd.template getProp<force_transport>(a_key)[xyz] += -1.0 * (Va2 + Vb2) * params.Pbackground * DW.get(xyz) / massa;
    }
    if (params.DENSITY_TYPE == DENSITY_DIFFERENTIAL)
    {
        vd.template getProp<drho>(a_key) += massb * dotProduct(v_rel, DW);
    }
}

template <typename vd_type, typename NN_type>
__global__ void calc_forcesGPU(vd_type vd,
                               NN_type NN,
                               const Parameters params)
{
    printf("Hello from inside calc forces gpu\n");
    const real_number dp = params.dp;

    auto a = GET_PARTICLE(vd);

    // if the particle is FLUID
    if (vd.template getProp<type>(a) == FLUID)
    {
        printf("PArticle is fluid\n");

        const Point<DIM, real_number> xa = vd.getPos(a);
        const real_number massa = params.MassFluid;
        const real_number rhoa = vd.template getProp<rho>(a);
        const real_number Pa = vd.template getProp<pressure>(a);
        Point<DIM, real_number> va = vd.template getProp<velocity>(a);

        // Reset the force counter (0 + gravity)
        for (int xyz = 0; xyz < DIM; xyz++)
        {
            vd.template getProp<force>(a)[xyz] = params.gravity_vector.get(xyz);
            vd.template getProp<force_transport>(a)[xyz] = 0.0f;
        }

        // Reset the density difference
        vd.template getProp<drho>(a) = 0.0;

        // Get an iterator over the neighborhood particles of p
        auto Np = NN.getNNIteratorBox(NN.getCell(xa));

        // For each neighborhood particle
        while (Np.isNext() == true)
        {
            printf("Hello from inside loop\n");

            // get key of particle b
            auto b = Np.get();

            // if (a == b) skip this particle
            if (a == b)
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

            // if they interact
            printf("r2: %f\n", r2);
            printf("r_threshold^2: %f\n", params.r_threshold * params.r_threshold);
            if (r2 < params.r_threshold * params.r_threshold && r2 >= 1e-16)
            {
                if (vd.template getProp<type>(b) == BOUNDARY)
                {
                    if (params.BC_TYPE == NO_SLIP)
                    {
                        interact_fluid_boundary_old(vd, a, massa, rhoa, Pa, xa, va, xb, dr, r2, b, params);
                    }
                    else if (params.BC_TYPE == NEW_NO_SLIP)
                    {
                        printf("Interacting fluid - boundary\n");
                        interact_fluid_boundary_new(vd, a, massa, rhoa, Pa, xa, va, xb, dr, b, params);
                    }
                }
                else
                {
                    printf("Interacting fluid - fluid\n");
                    interact_fluid_fluid(vd, a, massa, rhoa, Pa, xa, va, xb, dr, r2, b, params);
                }
            }
            ++Np;
        }
    }
}

template <typename NN_type>
inline void calc_forces(particles &vd,
                        NN_type &NN,
                        const Parameters &params)
{
    // get domain iterator
    auto it = vd.getDomainIteratorGPU();
    // // Update the cell-list
    vd.updateCellListGPU(NN);

    CUDA_LAUNCH(calc_forcesGPU, it, vd.toKernel(), NN.toKernel(), params);
    cudaError_t err = cudaDeviceSynchronize(); // Wait for the kernel to finish
    if (err != cudaSuccess)
    {
        printf("CUDA error: %s\n", cudaGetErrorString(err));
    }
}

// template <typename CellList>
// void calc_forces_regularize(particles &vd,
//                             CellList &NN,
//                             const Parameters &params)
// {
//     // get domain iterator
//     vector_dist_iterator part = vd.getDomainIterator();

//     // Update the cell-list
//     vd.updateCellList(NN);

//     // For each particle ...
//     while (part.isNext())
//     {
//         // get key of particle a
//         auto a = part.get();

//         // if the particle is FLUID
//         if (vd.getProp<type>(a) == FLUID)
//         {
//             const Point<DIM, real_number> xa = vd.getPos(a);
//             const real_number massa = params.MassFluid;
//             const real_number rhoa = vd.getProp<rho>(a);

//             // Reset the force counter (0 + gravity)
//             for (int xyz = 0; xyz < DIM; xyz++)
//             {
//                 vd.template getProp<force_transport>(a)[xyz] = 0.0;
//             }

//             // Reset the density difference
//             vd.template getProp<drho>(a) = 0.0;

//             // Get an iterator over the neighborhood particles of p
//             auto Np = NN.getNNIteratorBox(NN.getCell(xa));

//             // For each neighborhood particle
//             while (Np.isNext() == true)
//             {
//                 // get key of particle b
//                 unsigned long b = Np.get();

//                 // if (a == b) skip this particle
//                 if (a.getKey() == b)
//                 {
//                     ++Np;
//                     continue;
//                 }

//                 // Get the position xb of the particle
//                 Point<DIM, real_number> xb = vd.getPos(b);

//                 // Get the distance between a and b
//                 // in fluid - boundary its xf-xb i.e. vector pointing at fluid from boundary
//                 Point<DIM, real_number> dr = xa - xb;

//                 // take the norm (squared) of this vector
//                 real_number r2 = norm2(dr);

//                 // if they interact
//                 if (r2 < params.r_threshold * params.r_threshold)
//                 {
//                     if (vd.getProp<type>(b) == BOUNDARY)
//                     {
//                         if (params.BC_TYPE == NO_SLIP)
//                         {
//                             // interact_fluid_boundary_old(vd, a, massa, rhoa, Pa, xa, va, xb, dr, r2, b, accumulate_force, obstacle_force_x, obstacle_force_y, params);
//                             // NOT READY YET
//                         }
//                         else if (params.BC_TYPE == NEW_NO_SLIP)
//                             interact_fluid_boundary_new_regularize(vd, a, rhoa, dr, b, params);
//                     }
//                     else
//                     {
//                         interact_fluid_fluid_regularize(vd, a, massa, rhoa, dr, r2, b, params);
//                     }
//                 }
//                 ++Np;
//             }
//         }

//         ++part;
//     }
// }

#endif // CALC_FORCES_HPP
