#ifndef CALCULATIONS_H
#define CALCULATIONS_H

#include "Definitions.hpp"
#include "VectorUtilities.hpp"

template <typename CellList>
void CalcFluidVec(particles &vd, CellList &NN, const Parameters &params)
{
    // This function computes the vector pointing to the average position of the fluid particles for each boundary particle
    // Knowing this vector is useful to later compute the normal vector in the correct direction
    auto part = vd.getDomainIterator();

    // Update the cell-list
    vd.updateCellList(NN);

    // For each particle ...
    while (part.isNext())
    {
        // Key of the particle a
        vect_dist_key_dx a = part.get();

        // if particle BOUNDARY
        if (vd.getProp<type>(a) == BOUNDARY)
        {
            // Get the position xa of the particle a
            Point<DIM, real_number> xa = vd.getPos(a);

            // initialize sum
            Point<DIM, real_number> r_fluid_sum = {0.0, 0.0};

            // neighborhood particles
            auto Np = NN.getNNIteratorBox(NN.getCell(vd.getPos(a)));

            // iterate the neighborhood particles
            while (Np.isNext() == true)
            {
                // Key of b particle
                const unsigned long b = Np.get();

                if (vd.getProp<type>(b) == FLUID)
                {
                    // Get the position xb of the particle b
                    const Point<DIM, real_number> xb = vd.getPos(b);

                    // Get the vector pointing at fluid from boundary
                    Point<DIM, real_number> dr = xb - xa;

                    // take the norm squared of this vector
                    const real_number r2 = norm2(dr);

                    // If the particles interact ...
                    if (r2 < params.r_threshold * params.r_threshold)
                    {
                        normalizeVector(dr);
                        r_fluid_sum += dr;
                    }
                }

                ++Np;
            }
            normalizeVector(r_fluid_sum);

            // store in normal temporally, this vector is only needed for computing normals oriented "outside"
            // later is overwritten by the normal vector
            for (int xyz = 0; xyz < DIM; ++xyz)
            {
                vd.template getProp<normal_vector>(a)[xyz] = r_fluid_sum.get(xyz);
            }
        }

        ++part;
    }
}

template <typename CellList>
void CalcVorticity(particles &vd, CellList &NN, const Parameters &params)
{
    auto part = vd.getDomainIterator();

    // Update the cell-list
    vd.updateCellList(NN);

    // For each particle ...
    while (part.isNext())
    {
        // Key of the particle a
        vect_dist_key_dx a = part.get();

        // if particle FLUID
        if (vd.getProp<type>(a) == FLUID)
        {
            // Get the position xa of the particle a
            Point<DIM, real_number> xa = vd.getPos(a);

            // Get the velocity of the particle a
            Point<DIM, real_number> va = vd.getProp<velocity>(a);

            // initialize vorticity sum
            real_number vorticity = 0.0;

            // neighborhood particles
            auto Np = NN.getNNIteratorBox(NN.getCell(vd.getPos(a)));

            // iterate the neighborhood particles
            while (Np.isNext() == true)
            {
                // Key of b particle
                const unsigned long b = Np.get();

                // Only consider other fluid particles
                if (vd.getProp<type>(b) == FLUID)
                {
                    // Get the position ( and vel ) of particle b
                    const Point<DIM, real_number> xb = vd.getPos(b);

                    // Get the vector pointing at A from B
                    Point<DIM, real_number> dr = xa - xb;
                    // get the norm squared of this vector
                    const real_number r2 = norm2(dr);

                    if (r2 < params.r_threshold * params.r_threshold)
                    {
                        // Get the velocity of particle b
                        const Point<DIM, real_number> vb = vd.getProp<velocity>(b);
                        // evaluate the kernel gradient
                        Point<DIM, real_number> Dwab = DWab(dr, sqrt(r2), params.H, params.Kquintic);
                        vorticity += -params.MassFluid * crossProduct2D(vb - va, Dwab);
                    }
                }

                ++Np;
            }

            // Store in vorticity property
            vd.template getProp<vd_vorticity>(a) = vorticity / vd.getProp<rho>(a);
        }

        ++part;
    }
}

template <typename CellList>
void CalcNormalVec(particles &vd, CellList &NN, const Parameters &params)
{
    // This function computes the normal vector for a boundary particle based on the other boundary particles inside its kernel.
    // it computes the average (kernel weighted) of the vectors perpendicular to the vectors pointing to the other boundary particles ( oriented to the fluid )

    const real_number refinement = params.rf;
    const real_number global_dp = params.dp;
    // const real_number local_dp = params.dp / refinement;
    // const real_number local_H = params.H / refinement;
    const real_number local_dp = global_dp;
    const real_number local_H = params.H;
    const real_number local_Kquintic = (DIM == 3) ? 1.0 / 120.0 / M_PI / local_H / local_H / local_H : 7.0 / 478.0 / M_PI / local_H / local_H;

    auto part = vd.getDomainIterator();

    // Update the cell-list
    vd.updateCellList(NN);

    // For each particle ...
    while (part.isNext())
    {
        // Key of the particle a
        vect_dist_key_dx a = part.get();

        // if particle BOUNDARY
        if (vd.getProp<type>(a) == BOUNDARY)
        {
            // Get the position xa of the particle a
            Point<DIM, real_number> xa = vd.getPos(a);

            // get vector that points at fluid
            Point<DIM, real_number> Rfluid = vd.getProp<normal_vector>(a);

            // initialize sum
            Point<DIM, real_number> n_sum = {0.0};

            // neighborhood particles
            auto Np = NN.getNNIteratorBox(NN.getCell(vd.getPos(a)));

            // iterate the neighborhood particles
            while (Np.isNext() == true)
            {
                // Key of b particle
                const unsigned long b = Np.get();
                if (a.getKey() == b)
                {
                    ++Np;
                    continue;
                }

                if (vd.getProp<type>(b) == BOUNDARY)
                {
                    // Get the position xb of the particle b
                    const Point<DIM, real_number> xb = vd.getPos(b);

                    // Get the vector pointing at a from b
                    const Point<DIM, real_number> dr = xa - xb;

                    // take the norm squared of this vector
                    const real_number r2 = norm2(dr);

                    // If the particles interact ...
                    if (r2 < 9.0 * local_H * local_H)
                    {
                        // get perpendicular vector to dr
                        Point<DIM, real_number> perp = getPerpendicularUnit2D(dr); // this is normalized

                        // get scalar product of perp and Rfluid
                        const real_number perp_dot_fluid = dotProduct(perp, Rfluid);

                        // we want perp to point towards the fluid
                        if (perp_dot_fluid < 0.0)
                            perp = -1.0 * perp;

                        // evaluate kernel
                        real_number W = Wab(sqrt(r2), local_H, local_Kquintic);
                        n_sum += perp * W;
                    }
                }

                ++Np;
            }
            // normalize the summed vector
            normalizeVector(n_sum);

            // store in normal vector
            for (int xyz = 0; xyz < DIM; ++xyz)
            {
                vd.template getProp<normal_vector>(a)[xyz] = n_sum.get(xyz);
            }
        }

        ++part;
    }
}

template <typename CellList>
void CalcCurvature(particles &vd, CellList &NN, const Parameters &params)
{
    // This function computes the curvature of the boundary particles from the divergence of the normal vector

    const real_number refinement = params.rf;
    const real_number global_dp = params.dp;
    // const real_number local_dp = params.dp / refinement;
    // const real_number local_H = params.H / refinement;
    const real_number local_dp = global_dp;
    const real_number local_H = params.H;
    const real_number local_Kquintic = (DIM == 3) ? 1.0 / 120.0 / M_PI / local_H / local_H / local_H : 7.0 / 478.0 / M_PI / local_H / local_H;

    auto part = vd.getDomainIterator();

    // Update the cell-list
    vd.updateCellList(NN);

    // max curvature is determined form the derivation of the volume formula, for curvature higher than this the volume of the third particle
    // is no longer between the two circles. Actually for 1/(2.5*dp) it becomes 0 and later negative
    // const real_number max_curvature = 1.0 / (1.0 * global_dp);
    const real_number max_curvature = 1.0 / (3.0 * global_dp);

    // For each particle ...
    while (part.isNext())
    {
        // Key of the particle a
        vect_dist_key_dx a = part.get();

        // if particle BOUNDARY
        if (vd.getProp<type>(a) == BOUNDARY)
        {
            // Get the position xa of the particle a
            Point<DIM, real_number> xa = vd.getPos(a);

            // get normal of a
            Point<DIM, real_number> normal_a = vd.getProp<normal_vector>(a);

            // initialize sums
            real_number K_sum = 0.0;
            real_number w_sum = 0.0;

            auto Np = NN.getNNIteratorBox(NN.getCell(vd.getPos(a)));

            // iterate the neighborhood particles
            while (Np.isNext() == true)
            {
                // Key of b particle
                const unsigned long b = Np.get();

                if (vd.getProp<type>(b) == BOUNDARY)
                {
                    // Get the position xb of the particle b
                    const Point<DIM, real_number> xb = vd.getPos(b);

                    // Get the vector pointing at a from b
                    const Point<DIM, real_number> dr = xa - xb;

                    // take the norm squared of this vector
                    const real_number r2 = norm2(dr);

                    // If the particles interact ...
                    if (r2 < 9.0 * local_H * local_H)
                    {
                        if (a.getKey() != b)
                        {
                            // OLD CALCULATION
                            // Point<DIM, real_number> normal_b = vd.getProp<normal_vector>(b);
                            // real_number r = sqrt(r2);
                            // real_number W = Wab(r, local_H, local_Kquintic);
                            // Point<DIM, real_number> dW = DWab(dr, r, local_H, local_Kquintic);

                            // real_number local_k = dotProduct(normal_b - normal_a, dW); // divergence of the normal vector
                            // if (local_k > max_curvature)
                            // {
                            //     local_k = max_curvature;
                            // }

                            // K_sum += local_k;
                            // w_sum += W;

                            // NEW CALCULATION
                            Point<DIM, real_number> normal_b = vd.getProp<normal_vector>(b);
                            real_number r = sqrt(r2);
                            real_number W = Wab(r, local_H, local_Kquintic);
                            Point<DIM, real_number> eab = -1.0 * dr;
                            eab = eab / r;
                            real_number local_k = dotProduct(normal_b - normal_a, eab) / r;

                            if (local_k > max_curvature)
                            {
                                local_k = max_curvature;
                            }

                            K_sum += local_k * W;
                            w_sum += W;
                        }
                    }
                }

                ++Np;
            }

            K_sum = K_sum / w_sum;
            // store in curvature
            vd.template getProp<curvature_boundary>(a) = K_sum;
        }
        ++part;
    }
}

void CalcVolume(particles &vd, real_number dp)
{
    // This function computes the volume of the virtual particles for the new no-slip BC

    auto part = vd.getDomainIterator();

    // For each particle ...
    while (part.isNext())
    {
        // Key of the particle a
        vect_dist_key_dx a = part.get();

        // if particle BOUNDARY
        if (vd.getProp<type>(a) == BOUNDARY)
        {
            real_number kappa = vd.getProp<curvature_boundary>(a);
            real_number dxwall = vd.getProp<arc_length>(a);

            for (int i = 0; i < 3; i++)
            {
                vd.template getProp<vd_volume>(a)[i] = std::max(0.0, 0.5 * (2.0 * dp + dp * dp * kappa - 2.0 * (i + 1.0) * dp * dp * kappa) * dxwall);
            }
        }
        ++part;
    }
}

template <typename CellList>
void ExtrapolateVelocity(particles &vd, CellList &NN, const Parameters &params)
{
    // This function fills the value of v_transport for the boundary particles, with the extrapolated velocity for the old no slip BC

    auto part = vd.getDomainIterator();

    // Update the cell-list
    vd.updateCellList(NN);

    // For each particle ...
    while (part.isNext())
    {
        // Key of the b particle
        vect_dist_key_dx b = part.get();

        // if particle boundary
        if (vd.getProp<type>(b) != FLUID)
        {

            // Get the position xb of the boundary particle
            Point<DIM, real_number> xb = vd.getPos(b);

            // initialize sums
            Point<DIM, real_number> sum_vW = (DIM == 2) ? Point<DIM, real_number>{0.0, 0.0} : Point<DIM, real_number>{0.0, 0.0, 0.0};

            real_number sum_pW = 0.0;
            real_number sum_W = 0.0;

            auto Np = NN.getNNIteratorBox(NN.getCell(vd.getPos(b)));

            // iterate the neighborhood fluid particles
            while (Np.isNext() == true)
            {
                // Key of fluid particle
                unsigned long f = Np.get();

                // if (b == f) skip this particle
                if (b.getKey() == f)
                {
                    ++Np;
                    continue;
                };
                // Skip other boundary particles
                if (vd.getProp<type>(f) != FLUID)
                {
                    ++Np;
                    continue;
                }

                // Get the position xf of the fluid particle
                Point<DIM, real_number> xf = vd.getPos(f);

                // Get the velocity of the fluid particle
                Point<DIM, real_number> vf = vd.getProp<velocity>(f);

                // Get the density of the fluid particle
                real_number rhof = vd.getProp<rho>(f);

                // Get the pressure of the fluid particle
                real_number Pf = vd.getProp<pressure>(f);

                // Get the vector pointing at xb from xf rwf
                Point<DIM, real_number> dr = xb - xf;

                // take the norm squared of this vector
                real_number r2 = norm2(dr);

                // If the particles interact ...
                if (r2 < params.r_threshold * params.r_threshold)
                {
                    // calculate distance
                    real_number r = sqrt(r2);

                    // evaluate kernel
                    real_number w = Wab(r, params.H, params.Kquintic);
                    // compute v*W
                    sum_vW += w * vf;

                    // compute Pf +rhof*(g-a) dot rwf
                    // at the moment a = 0
                    const real_number dot = dotProduct(dr, params.gravity_vector);
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
                const Point<DIM, real_number> radial_vec = {xb.get(0) - vd.getProp<force_transport>(b)[0], xb.get(1) - vd.getProp<force_transport>(b)[1]};
                const real_number radial = getVectorNorm(radial_vec);
                // vector tangential to the radial vector, rotation velocity is in this direction
                const Point<DIM, real_number> rotation_tangential = getPerpendicularUnit2D(radial_vec);
                Point<DIM, real_number> vw = vd.template getProp<velocity>(b);
                vw.get(0) += radial * vd.getProp<vd_omega>(b) * rotation_tangential.get(0);
                vw.get(1) += radial * vd.getProp<vd_omega>(b) * rotation_tangential.get(1);

                for (int xyz = 0; xyz < DIM; ++xyz)
                {
                    vd.template getProp<v_transport>(b)[xyz] = 2.0 * vw[xyz] - sum_vW.get(xyz) / sum_W;
                }

                // Set the pressure of the boundary particle b
                vd.template getProp<pressure>(b) = sum_pW / sum_W;
                // Compute density from inverted Eq of state
                vd.template getProp<rho>(b) = InvEqState_particle(vd.template getProp<pressure>(b), params.rho_zero, params.B, params.gamma_, params.xi);
            }
            else
            {

                for (int xyz = 0; xyz < DIM; ++xyz)
                {
                    vd.template getProp<v_transport>(b)[xyz] = 2.0 * vd.template getProp<velocity>(b)[xyz];
                }

                vd.template getProp<pressure>(b) = 0.0;
                vd.template getProp<rho>(b) = params.rho_zero;
            }
        }

        ++part;
    }
}

#endif // CALCULATIONS_H
