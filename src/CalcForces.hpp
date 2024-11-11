#ifndef CALC_FORCES_HPP
#define CALC_FORCES_HPP

#include "Definitions.hpp"
#include "VectorUtilities.hpp"
#include "Kernel.hpp"
#include "Physics.hpp"
#include <cstdio>
#include <stdio.h>
#include <unistd.h>

template <typename vd_type>
inline __device__ void interact_fluid_fluid(vd_type &vd, unsigned int b, unsigned int a,
                                            Point<DIM, real_number> va,
                                            Point<DIM, real_number> dr,
                                            real_number r2,
                                            real_number massa,
                                            real_number massb,
                                            real_number rhoa,
                                            real_number Pa)
{
    const real_number rhob = vd.template getProp<vd1_rho>(b);
    const real_number Pb = vd.template getProp<vd2_pressure>(b);
    const Point<DIM, real_number> vb = vd.template getProp<vd4_velocity>(b);

    const real_number r = sqrt(r2);

    const Point<DIM, real_number> v_rel = va - vb;

    const Point<DIM, real_number> DW = DWab(dr, r, _params_gpu_.H, _params_gpu_.Kquintic);

    const real_number Va2 = (massa / rhoa) * (massa / rhoa);
    const real_number Vb2 = (massb / rhob) * (massb / rhob);

    const Point<DIM, real_number> ViscosityTerm = Pi_physical(dr, r, v_rel, DW, _params_gpu_.eta);
    const real_number PressureTerm = PressureForce(rhoa, rhob, Pa, Pb);

    const Point<DIM, real_number> vta = vd.template getProp<vd5_velocity_t>(a);
    const Point<DIM, real_number> vtb = vd.template getProp<vd5_velocity_t>(b);
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
        vd.template getProp<vd6_force>(a)[xyz] += (Va2 + Vb2) * (PressureTerm * DW.get(xyz) + ViscosityTerm.get(xyz) + DivATerm.get(xyz)) / massa;
        vd.template getProp<vd7_force_t>(a)[xyz] += -1.0 * (Va2 + Vb2) * _params_gpu_.Pbackground * DW.get(xyz) / massa;
    }
    // if (_params_gpu_.DENSITY_TYPE == DENSITY_DIFFERENTIAL)
    // {
    vd.template getProp<vd3_drho>(a) += rhoa * massb * dotProduct(v_rel, DW) / rhob;
    // }
}

template <typename vd_type, typename NN_type>
__global__ void calc_forcesGPU_new2D(vd_type vd,
                                     NN_type NN)
{
    // auto a = GET_PARTICLE(vd); // old unsorted code
    unsigned int a;
    GET_PARTICLE_SORT(a, NN);

    const real_number dp = _params_gpu_.dp;

    // if the particle is FLUID
    if (vd.template getProp<vd0_type>(a) == FLUID)
    {

        const Point<DIM, real_number> xa = vd.getPos(a);
        const real_number massa = _params_gpu_.MassFluid;
        const real_number massb = _params_gpu_.MassBound;

        const real_number rhoa = vd.template getProp<vd1_rho>(a);
        const real_number Pa = vd.template getProp<vd2_pressure>(a);
        Point<DIM, real_number> va = vd.template getProp<vd4_velocity>(a);

        // Reset the force counter (0 + gravity)
        for (int xyz = 0; xyz < DIM; xyz++)
        {
            vd.template getProp<vd6_force>(a)[xyz] = _params_gpu_.gravity_vector[xyz];
            vd.template getProp<vd7_force_t>(a)[xyz] = 0.0;
        }

        // reset accumulator for drag and lift
        vd.template getProp<vd13_force_red_x>(a) = 0.0;
        vd.template getProp<vd14_force_red_y>(a) = 0.0;

        // Reset the density difference
        vd.template getProp<vd3_drho>(a) = 0.0;

        // Get an iterator over the neighborhood particles of p
        auto Np = NN.getNNIteratorBox(NN.getCell(xa));

        // For each neighborhood particle
        while (Np.isNext() == true)
        {

            auto b = Np.get_sort();

            // Get the position xb of the particle
            Point<DIM, real_number> xb = vd.getPos(b);

            // Get the distance between a and b
            // in fluid - boundary its xf-xb i.e. vector pointing at fluid from boundary
            Point<DIM, real_number> dr = xa - xb;

            // take the norm (squared) of this vector
            real_number r2 = norm2(dr);

            if (r2 < _params_gpu_.r_cut2 && r2 > 1e-16)
            {
                auto typeb = vd.template getProp<vd0_type>(b);
                if (typeb == BOUNDARY || typeb == OBSTACLE)
                {
                    // Points from fluid to wall
                    Point<DIM, real_number> r_fluid_to_wall = -1.0 * dr;
                    // real_number dist2marker = getVectorNorm(r_fluid_to_wall);
                    Point<DIM, real_number> vw = vd.template getProp<vd4_velocity>(b);

                    real_number ang_vel = vd.template getProp<vd10_omega>(b);

                    // marker particles store centre of solid body in force_transport since it is unused
                    // vector pointing from centre of rotation to marker particle
                    const Point<DIM, real_number> radial_vec = {xb.get(0) - vd.template getProp<vd7_force_t>(b)[0],
                                                                xb.get(1) - vd.template getProp<vd7_force_t>(b)[1]};
                    const real_number radius = getVectorNorm(radial_vec);
                    // get vector tangential to the radial vector, rotation velocity is in this direction
                    const Point<DIM, real_number> tangential_rotation = getPerpendicularUnit2D(radial_vec);

                    // Wall velocity is linear velocity + w*R*tangential
                    vw.get(0) += radius * ang_vel * tangential_rotation.get(0);
                    vw.get(1) += radius * ang_vel * tangential_rotation.get(1);

                    // Get normal and tangential vectors for velocity mirroring
                    const Point<DIM, real_number> normal = vd.template getProp<vd8_normal>(b);
                    const Point<DIM, real_number> tangential = getPerpendicularUnit2D(normal);

                    // wall acceleration
                    const Point<DIM, real_number> aw = vd.template getProp<vd6_force>(b);

                    // Difference between fluid transport and momentum velocity
                    const Point<DIM, real_number> vtf = vd.template getProp<vd5_velocity_t>(a);
                    const Point<DIM, real_number> vdiff_f = vtf - va;

                    // Project va and vw on tangential and normal directions
                    real_number vt = dotProduct(va, tangential);
                    real_number vn = dotProduct(va, normal);
                    real_number vwt = dotProduct(vw, tangential);
                    real_number vwn = dotProduct(vw, normal);

                    // vertical distance from fluid particle to wall
                    real_number lf = dotProduct(r_fluid_to_wall, normal);
                    lf = (lf < 0.0 ? -1.0 * lf : lf); // absolute value

                    // Get array of vectors from fluid to 3 boundary particles and its norms
                    std::array<Point<DIM, real_number>, 3> r_boundary = getBoundaryPositions(r_fluid_to_wall, normal, dp);
                    std::array<real_number, 3> r_boundary_norm = {getVectorNorm(r_boundary[0]), getVectorNorm(r_boundary[1]), getVectorNorm(r_boundary[2])};

                    // const real_number dist2third = r_boundary_norm[2];
                    // distance from 3 boundary particles to marker
                    std::array<real_number, 3> lwall = {0.5 * dp, 1.5 * dp, 2.5 * dp};

                    // project forces on normal direction
                    real_number g_normal = dotProduct(_params_gpu_.gravity_vector, normal);
                    real_number a_normal = dotProduct(aw, normal);

                    // For gradient of Af tensor
                    const std::array<Point<DIM, real_number>, DIM> Af = dyadicProduct(rhoa * va, vdiff_f);
                    // to avoid division by zero
                    lf = std::max(lf, 0.25 * dp);

                    const Point<3, real_number> vol = vd.template getProp<vd9_volume>(b);

                    for (int i = 0; i < 3; i++) // for the 3 boundary particles
                    {
                        if (r_boundary_norm[i] < _params_gpu_.r_cut)
                        {

                            // const real_number Mass_boundary = vol[i] * _params_gpu_.rho0;

                            const Point<DIM, real_number> v_boundary_visc = ((vwt - vt) * (lwall[i] / lf) + vwt) * tangential + vn * normal;
                            const Point<DIM, real_number> v_boundary_cont = ((vwn - vn) * (lwall[i] / lf) + vwn) * normal + vt * tangential;
                            const real_number p_boundary = Pa - rhoa * (g_normal - a_normal) * (lf + lwall[i]); //  dot(r_boundary,normal) = -(lf+lw)
                            const real_number rho_boundary = InvEqState_particle(p_boundary, _params_gpu_.rho0, _params_gpu_.B, _params_gpu_.gamma, _params_gpu_.xi);

                            // flip sign of r_boundary to get vector pointing from boundary to fluid (Force routines use the vector pointing from b to a)
                            r_boundary[i] = -1.0 * r_boundary[i];

                            // Evaluate kernel gradient
                            const Point<DIM, real_number> DW = DWab(r_boundary[i], r_boundary_norm[i], _params_gpu_.H, _params_gpu_.Kquintic);

                            // Compute forces
                            const Point<DIM, real_number> v_rel_visc = va - v_boundary_visc;
                            const Point<DIM, real_number> v_rel_cont = va - v_boundary_cont;

                            // const real_number Vb = (Mass_boundary / rho_boundary); // vb is mass/rho instead of directly vol[i] because it allows to variate with density
                            const real_number Vb = vol[i];

                            const Point<DIM, real_number> ViscosityTerm = Pi_physical(r_boundary[i], r_boundary_norm[i], v_rel_visc, DW, _params_gpu_.eta);
                            const real_number PressureTerm = PressureForce(rhoa, rho_boundary, Pa, p_boundary); //-p_boundary - Pf;
                            const Point<DIM, real_number> DivATerm = 0.5 * matVec(Af, DW);

                            Point<DIM, real_number> force_tmp;
                            for (int xyz = 0; xyz < DIM; ++xyz)
                            {
                                force_tmp[xyz] = 2.0 * (Vb / rhoa) * (PressureTerm * DW.get(xyz) + ViscosityTerm.get(xyz) + DivATerm.get(xyz));
                                vd.template getProp<vd6_force>(a)[xyz] += force_tmp[xyz];
                                vd.template getProp<vd7_force_t>(a)[xyz] += -2.0 * (Vb / rhoa) * (_params_gpu_.Pbackground) * DW.get(xyz);
                            }

                            real_number scal = (typeb == OBSTACLE) ? 1.0 : 0.0;
                            vd.template getProp<vd13_force_red_x>(a) += -scal * force_tmp[0];
                            vd.template getProp<vd14_force_red_y>(a) += -scal * force_tmp[1];

                            // if (_params_gpu_.DENSITY_TYPE == DENSITY_DIFFERENTIAL)
                            // {
                            vd.template getProp<vd3_drho>(a) += rhoa * Vb * dotProduct(v_rel_cont, DW);
                            // }
                        }
                    }
                }
                else // INTERACT FLUID FLUID
                {
                    interact_fluid_fluid(vd, b, a, va, dr, r2, massa, massb, rhoa, Pa);
                }
            }
            ++Np;
        }
    }
}

template <typename vd_type, typename NN_type>
__global__ void calc_forcesGPU_new3D(vd_type vd,
                                     NN_type NN)
{
    // auto a = GET_PARTICLE(vd); // old unsorted code
    unsigned int a;
    GET_PARTICLE_SORT(a, NN);

    const real_number dp = _params_gpu_.dp;

    // if the particle is FLUID
    if (vd.template getProp<vd0_type>(a) == FLUID)
    {

        const Point<DIM, real_number> xa = vd.getPos(a);
        const real_number massa = _params_gpu_.MassFluid;
        const real_number massb = _params_gpu_.MassBound;

        const real_number rhoa = vd.template getProp<vd1_rho>(a);
        const real_number Pa = vd.template getProp<vd2_pressure>(a);
        Point<DIM, real_number> va = vd.template getProp<vd4_velocity>(a);

        // Reset the force counter (0 + gravity)
        for (int xyz = 0; xyz < DIM; xyz++)
        {
            vd.template getProp<vd6_force>(a)[xyz] = _params_gpu_.gravity_vector[xyz];
            vd.template getProp<vd7_force_t>(a)[xyz] = 0.0;
        }

        // reset accumulator for drag and lift
        vd.template getProp<vd13_force_red_x>(a) = 0.0;
        vd.template getProp<vd14_force_red_y>(a) = 0.0;

        // Reset the density difference
        vd.template getProp<vd3_drho>(a) = 0.0;

        // Get an iterator over the neighborhood particles of p
        auto Np = NN.getNNIteratorBox(NN.getCell(xa));

        // For each neighborhood particle
        while (Np.isNext() == true)
        {

            auto b = Np.get_sort();

            // Get the position xb of the particle
            Point<DIM, real_number> xb = vd.getPos(b);

            // Get the distance between a and b
            // in fluid - boundary its xf-xb i.e. vector pointing at fluid from boundary
            Point<DIM, real_number> dr = xa - xb;

            // take the norm (squared) of this vector
            real_number r2 = norm2(dr);

            if (r2 < _params_gpu_.r_cut2 && r2 > 1e-16)
            {
                auto typeb = vd.template getProp<vd0_type>(b);
                if (typeb == BOUNDARY || typeb == OBSTACLE)
                {
                    // Points from fluid to wall
                    Point<DIM, real_number> r_fluid_to_wall = -1.0 * dr;
                    // real_number dist2marker = getVectorNorm(r_fluid_to_wall);
                    Point<DIM, real_number> vw = vd.template getProp<vd4_velocity>(b);

                    // Get normal and tangential vectors for velocity mirroring
                    const Point<DIM, real_number> normal = vd.template getProp<vd8_normal>(b);
                    Point<DIM, real_number> tangential1, tangential2;

                    getPerpendicularUnit3D(normal, tangential1, tangential2);

                    // wall acceleration
                    const Point<DIM, real_number> aw = vd.template getProp<vd6_force>(b);

                    // Difference between fluid transport and momentum velocity
                    const Point<DIM, real_number> vtf = vd.template getProp<vd5_velocity_t>(a);
                    const Point<DIM, real_number> vdiff_f = vtf - va;

                    // Project va and vw on tangential and normal directions
                    real_number vt1 = dotProduct(va, tangential1);
                    real_number vt2 = dotProduct(va, tangential2);
                    real_number vn = dotProduct(va, normal);
                    real_number vwt1 = dotProduct(vw, tangential1);
                    real_number vwt2 = dotProduct(vw, tangential2);
                    real_number vwn = dotProduct(vw, normal);

                    // vertical distance from fluid particle to wall
                    real_number lf = dotProduct(r_fluid_to_wall, normal);
                    lf = (lf < 0.0 ? -1.0 * lf : lf); // absolute value

                    // Get array of vectors from fluid to 3 boundary particles and its norms
                    std::array<Point<DIM, real_number>, 3> r_boundary = getBoundaryPositions(r_fluid_to_wall, normal, dp);
                    std::array<real_number, 3> r_boundary_norm = {getVectorNorm(r_boundary[0]), getVectorNorm(r_boundary[1]), getVectorNorm(r_boundary[2])};

                    // const real_number dist2third = r_boundary_norm[2];
                    // distance from 3 boundary particles to marker
                    std::array<real_number, 3> lwall = {0.5 * dp, 1.5 * dp, 2.5 * dp};

                    // project forces on normal direction
                    real_number g_normal = dotProduct(_params_gpu_.gravity_vector, normal);
                    real_number a_normal = dotProduct(aw, normal);

                    // For gradient of Af tensor
                    const std::array<Point<DIM, real_number>, DIM> Af = dyadicProduct(rhoa * va, vdiff_f);
                    // to avoid division by zero
                    lf = std::max(lf, 0.25 * dp);

                    const Point<3, real_number> vol = vd.template getProp<vd9_volume>(b);

                    for (int i = 0; i < 3; i++) // for the 3 boundary particles
                    {
                        if (r_boundary_norm[i] < _params_gpu_.r_cut)
                        {

                            // const real_number Mass_boundary = vol[i] * _params_gpu_.rho0;

                            const Point<DIM, real_number> v_boundary_visc = ((vwt1 - vt1) * (lwall[i] / lf) + vwt1) * tangential1 + ((vwt2 - vt2) * (lwall[i] / lf) + vwt2) * tangential2 + vn * normal;
                            const Point<DIM, real_number> v_boundary_cont = ((vwn - vn) * (lwall[i] / lf) + vwn) * normal + vt1 * tangential1 + vt2 * tangential2;
                            const real_number p_boundary = Pa - rhoa * (g_normal - a_normal) * (lf + lwall[i]); //  dot(r_boundary,normal) = -(lf+lw)
                            const real_number rho_boundary = InvEqState_particle(p_boundary, _params_gpu_.rho0, _params_gpu_.B, _params_gpu_.gamma, _params_gpu_.xi);

                            // flip sign of r_boundary to get vector pointing from boundary to fluid (Force routines use the vector pointing from b to a)
                            r_boundary[i] = -1.0 * r_boundary[i];

                            // Evaluate kernel gradient
                            const Point<DIM, real_number> DW = DWab(r_boundary[i], r_boundary_norm[i], _params_gpu_.H, _params_gpu_.Kquintic);

                            // Compute forces
                            const Point<DIM, real_number> v_rel_visc = va - v_boundary_visc;
                            const Point<DIM, real_number> v_rel_cont = va - v_boundary_cont;

                            // const real_number Vb = (Mass_boundary / rho_boundary); // vb is mass/rho instead of directly vol[i] because it allows to variate with density

                            const real_number Vb = vol[i];

                            const Point<DIM, real_number> ViscosityTerm = Pi_physical(r_boundary[i], r_boundary_norm[i], v_rel_visc, DW, _params_gpu_.eta);
                            const real_number PressureTerm = PressureForce(rhoa, rho_boundary, Pa, p_boundary); //-p_boundary - Pf;
                            const Point<DIM, real_number> DivATerm = 0.5 * matVec(Af, DW);

                            Point<DIM, real_number> force_tmp;
                            for (int xyz = 0; xyz < DIM; ++xyz)
                            {
                                force_tmp[xyz] = 2.0 * (Vb / rhoa) * (PressureTerm * DW.get(xyz) + ViscosityTerm.get(xyz) + DivATerm.get(xyz));
                                vd.template getProp<vd6_force>(a)[xyz] += force_tmp[xyz];
                                vd.template getProp<vd7_force_t>(a)[xyz] += -2.0 * (Vb / rhoa) * (_params_gpu_.Pbackground) * DW.get(xyz);
                            }

                            real_number scal = (typeb == OBSTACLE) ? 1.0 : 0.0;
                            vd.template getProp<vd13_force_red_x>(a) += -scal * force_tmp[0];
                            vd.template getProp<vd14_force_red_y>(a) += -scal * force_tmp[1];

                            // if (_params_gpu_.DENSITY_TYPE == DENSITY_DIFFERENTIAL)
                            // {
                            vd.template getProp<vd3_drho>(a) += rhoa * Vb * dotProduct(v_rel_cont, DW);
                            // }
                        }
                    }
                }
                else // INTERACT FLUID FLUID
                {
                    interact_fluid_fluid(vd, b, a, va, dr, r2, massa, massb, rhoa, Pa);
                }
            }
            ++Np;
        }
    }
}

template <typename vd_type, typename NN_type>
__global__ void calc_forcesGPU_old(vd_type vd,
                                   NN_type NN)
{

    unsigned int a;
    GET_PARTICLE_SORT(a, NN);

    // if the particle is FLUID
    if (vd.template getProp<vd0_type>(a) == FLUID)
    {

        const Point<DIM, real_number> xa = vd.getPos(a);
        const real_number massa = _params_gpu_.MassFluid;
        const real_number massb = _params_gpu_.MassBound;

        const real_number rhoa = vd.template getProp<vd1_rho>(a);
        const real_number Pa = vd.template getProp<vd2_pressure>(a);
        Point<DIM, real_number> va = vd.template getProp<vd4_velocity>(a);

        // Reset the force counter (0 + gravity)
        for (int xyz = 0; xyz < DIM; xyz++)
        {
            vd.template getProp<vd6_force>(a)[xyz] = _params_gpu_.gravity_vector[xyz];
            vd.template getProp<vd7_force_t>(a)[xyz] = 0.0;
        }

        // Reset the density difference
        vd.template getProp<vd3_drho>(a) = 0.0;

        vd.template getProp<vd13_force_red_x>(a) = 0.0;
        vd.template getProp<vd14_force_red_y>(a) = 0.0;

        // Get an iterator over the neighborhood particles of p
        auto Np = NN.getNNIteratorBox(NN.getCell(xa));

        // For each neighborhood particle
        while (Np.isNext() == true)
        {
            auto b = Np.get_sort();

            // Get the position xb of the particle
            Point<DIM, real_number> xb = vd.getPos(b);

            // Get the distance between a and b
            // in fluid - boundary its xf-xb i.e. vector pointing at fluid from boundary
            Point<DIM, real_number> dr = xa - xb;

            // take the norm (squared) of this vector
            real_number r2 = norm2(dr);

            if (r2 < _params_gpu_.r_cut2 && r2 > 1e-16)
            {
                auto typeb = vd.template getProp<vd0_type>(b);
                if (typeb == BOUNDARY || typeb == OBSTACLE) // FLUID - BOUNDARY INTERACTION
                {
                    const real_number rhob = vd.template getProp<vd1_rho>(b);
                    const real_number Pb = vd.template getProp<vd2_pressure>(b);
                    const Point<DIM, real_number> vb = vd.template getProp<vd4_velocity>(b);
                    const Point<DIM, real_number> vb_noslip = vd.template getProp<vd5_velocity_t>(b); // here we store the extrapolated velocity for no slip BC

                    const real_number r = sqrt(r2);

                    const Point<DIM, real_number> v_rel = va - vb;
                    const Point<DIM, real_number> v_rel_aux = va - vb_noslip;

                    const Point<DIM, real_number> DW = DWab(dr, r, _params_gpu_.H, _params_gpu_.Kquintic);

                    const real_number Va2 = (massa / rhoa) * (massa / rhoa);
                    const real_number Vb2 = (massb / rhob) * (massb / rhob);

                    const Point<DIM, real_number> ViscosityTerm = Pi_physical(dr, r, v_rel_aux, DW, _params_gpu_.eta);
                    const real_number PressureTerm = PressureForce(rhoa, rhob, Pa, Pb);

                    const Point<DIM, real_number> vtf = vd.template getProp<vd5_velocity_t>(a);
                    const Point<DIM, real_number> vdiff_f = vtf - va;
                    // Boundary particles have no transport velocity difference
                    std::array<Point<DIM, real_number>, DIM> Af = dyadicProduct(rhoa * va, vdiff_f);
                    Point<DIM, real_number> DivATerm = 0.5 * matVec(Af, DW);

                    Point<DIM, real_number> force_tmp;

                    for (int xyz = 0; xyz < DIM; ++xyz)
                    {
                        force_tmp[xyz] = (Va2 + Vb2) * (PressureTerm * DW.get(xyz) + ViscosityTerm.get(xyz) + DivATerm.get(xyz)) / massa;
                        vd.template getProp<vd6_force>(a)[xyz] += force_tmp[xyz];
                        vd.template getProp<vd7_force_t>(a)[xyz] += -1.0 * (Va2 + Vb2) * _params_gpu_.Pbackground * DW.get(xyz) / massa;
                    }
                    real_number scal = (typeb == OBSTACLE) ? 1.0 : 0.0;
                    vd.template getProp<vd13_force_red_x>(a) += -scal * force_tmp[0];
                    vd.template getProp<vd14_force_red_y>(a) += -scal * force_tmp[1];

                    // if (_params_gpu_.DENSITY_TYPE == DENSITY_DIFFERENTIAL)
                    // {
                    vd.template getProp<vd3_drho>(a) += rhoa * massb * dotProduct(v_rel_aux, DW) / rhob;
                    // }
                }
                else // INTERACT FLUID FLUID
                {
                    interact_fluid_fluid(vd, b, a, va, dr, r2, massa, massb, rhoa, Pa);
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
    auto it = vd.getDomainIteratorGPU(32);

    if (params.BC_TYPE == NO_SLIP)
    {
        vd.template updateCellListGPU<vd0_type, vd1_rho, vd2_pressure, vd3_drho, vd4_velocity, vd5_velocity_t, vd6_force, vd7_force_t, vd8_normal, vd9_volume, vd10_omega, vd13_force_red_x, vd14_force_red_y>(NN);
        CUDA_LAUNCH(calc_forcesGPU_old, it, vd.toKernel(), NN.toKernel());
    }
    else if (params.BC_TYPE == NEW_NO_SLIP)
    {
        vd.template updateCellListGPU<vd0_type, vd1_rho, vd2_pressure, vd3_drho, vd4_velocity, vd5_velocity_t, vd6_force, vd7_force_t, vd8_normal, vd9_volume, vd10_omega, vd13_force_red_x, vd14_force_red_y>(NN);
        if constexpr (DIM == 2)
            CUDA_LAUNCH(calc_forcesGPU_new2D, it, vd.toKernel(), NN.toKernel());
        else if constexpr (DIM == 3)
            CUDA_LAUNCH(calc_forcesGPU_new3D, it, vd.toKernel(), NN.toKernel());
    }

    vd.template restoreOrder<vd0_type, vd1_rho, vd2_pressure, vd3_drho, vd4_velocity, vd5_velocity_t, vd6_force, vd7_force_t, vd8_normal, vd9_volume, vd10_omega, vd13_force_red_x, vd14_force_red_y>(NN);
}

#endif // CALC_FORCES_HPP
