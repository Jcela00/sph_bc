#ifndef CALCULATIONS_H
#define CALCULATIONS_H

#include "Definitions.hpp"
#include "VectorUtilities.hpp"

template <typename CellList>
void CalcFluidVec(particles &vd, CellList &NN, const Parameters &params)
{
    // This function computes the vector pointing to the average position of the fluid particles with respect to the boundary particle

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
            Point<DIM, double> xa = vd.getPos(a);

            // initialize sum
            Point<DIM, double> r_fluid_sum = {0.0, 0.0};

            // neighborhood particles
            auto Np = NN.getNNIterator(NN.getCell(vd.getPos(a)));

            // iterate the neighborhood particles
            while (Np.isNext() == true)
            {
                // Key of b particle
                const unsigned long b = Np.get();

                if (vd.getProp<type>(b) == FLUID)
                {
                    // Get the position xb of the particle b
                    const Point<DIM, double> xb = vd.getPos(b);

                    // Get the vector pointing at fluid from boundary
                    Point<DIM, double> dr = xb - xa;

                    // take the norm squared of this vector
                    const double r2 = norm2(dr);

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

            for (int xyz = 0; xyz < DIM; ++xyz)
            {
                vd.template getProp<normal_vector>(a)[xyz] = r_fluid_sum.get(xyz);
            }
            // store in normal temporally, this vector is only needed for computing normals oriented "outside"
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
            Point<DIM, double> xa = vd.getPos(a);
            Point<DIM, double> va = vd.getProp<velocity>(a);

            // initialize vorticity sum
            double vorticity = 0.0;

            // neighborhood particles
            auto Np = NN.getNNIterator(NN.getCell(vd.getPos(a)));

            // iterate the neighborhood particles
            while (Np.isNext() == true)
            {
                // Key of b particle
                const unsigned long b = Np.get();

                // Only consider other fluid particles
                if (vd.getProp<type>(b) == FLUID)
                {
                    // Get the position ( and vel ) of particle b
                    const Point<DIM, double> xb = vd.getPos(b);

                    // Get the vector pointing at A from B
                    Point<DIM, double> dr = xa - xb;
                    // get the norm squared of this vector
                    const double r2 = norm2(dr);

                    if (r2 < params.r_threshold * params.r_threshold)
                    {
                        // Get the velocity of particle b
                        const Point<DIM, double> vb = vd.getProp<velocity>(b);
                        // evaluate the kernel gradient
                        Point<DIM, double> Dwab = DWab(dr, sqrt(r2), params.H, params.Kquintic);
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
    // it computes the average of the perpendicular vectors to the vectors pointing to the other boundary particles ( oriented to the fluid )

    const double refinement = params.rf;
    const double global_dp = params.dp;
    // const double local_dp = params.dp / refinement;
    // const double local_H = params.H / refinement;
    const double local_dp = global_dp;
    const double local_H = params.H;
    const double local_Kquintic = (DIM == 3) ? 1.0 / 120.0 / M_PI / local_H / local_H / local_H : 7.0 / 478.0 / M_PI / local_H / local_H;

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
            Point<DIM, double> xa = vd.getPos(a);

            // get vector that points at fluid
            Point<DIM, double> Rfluid = vd.getProp<normal_vector>(a);

            // initialize sum
            Point<DIM, double> n_sum = {0.0};

            // neighborhood particles
            auto Np = NN.getNNIterator(NN.getCell(vd.getPos(a)));

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
                    const Point<DIM, double> xb = vd.getPos(b);

                    // Get the vector pointing at a from b
                    const Point<DIM, double> dr = xa - xb;

                    // take the norm squared of this vector
                    const double r2 = norm2(dr);

                    // If the particles interact ...
                    if (r2 < 9.0 * local_H * local_H)
                    {
                        // get perpendicular vector to dr
                        Point<DIM, double> perp = getPerpendicularUnit2D(dr); // this is normalized

                        // get scalar product of perp and Rfluid
                        const double perp_dot_fluid = dotProduct(perp, Rfluid);

                        // we want perp to point towards the fluid
                        if (perp_dot_fluid < 0.0)
                            perp = -1.0 * perp;

                        // evaluate kernel
                        double W = Wab(sqrt(r2), local_H, local_Kquintic);
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
void CalcCurvature(particles &vd, CellList &NN, Vcluster<> &v_cl, const Parameters &params)
{
    // This function computes the curvature of the boundary particles from the divergence of the normal vector

    const double refinement = params.rf;
    const double global_dp = params.dp;
    // const double local_dp = params.dp / refinement;
    // const double local_H = params.H / refinement;
    const double local_dp = global_dp;
    const double local_H = params.H;
    const double local_Kquintic = (DIM == 3) ? 1.0 / 120.0 / M_PI / local_H / local_H / local_H : 7.0 / 478.0 / M_PI / local_H / local_H;

    auto part = vd.getDomainIterator();

    // Update the cell-list
    vd.updateCellList(NN);

    // max curvature is determined form the derivation of the volume formula, for curvature higher than this the volume of the third particle
    // is no longer between the two circles. Actually for 1/(2.5*dp) it becomes 0 and later negative
    // const double max_curvature = 1.0 / (1.0 * global_dp);
    const double max_curvature = 1.0 / (3.0 * global_dp);

    // For each particle ...
    while (part.isNext())
    {
        // Key of the particle a
        vect_dist_key_dx a = part.get();

        // if particle BOUNDARY
        if (vd.getProp<type>(a) == BOUNDARY)
        {
            // Get the position xa of the particle a
            Point<DIM, double> xa = vd.getPos(a);

            // get normal of a
            Point<DIM, double> normal_a = vd.getProp<normal_vector>(a);

            // initialize sums
            double K_sum = 0.0;
            double w_sum = 0.0;

            auto Np = NN.getNNIterator(NN.getCell(vd.getPos(a)));

            // iterate the neighborhood particles
            while (Np.isNext() == true)
            {
                // Key of b particle
                const unsigned long b = Np.get();

                if (vd.getProp<type>(b) == BOUNDARY)
                {
                    // Get the position xb of the particle b
                    const Point<DIM, double> xb = vd.getPos(b);

                    // Get the vector pointing at a from b
                    const Point<DIM, double> dr = xa - xb;

                    // take the norm squared of this vector
                    const double r2 = norm2(dr);

                    // If the particles interact ...
                    if (r2 < 9.0 * local_H * local_H)
                    {
                        if (a.getKey() != b)
                        {
                            // OLD CALCULATION
                            // Point<DIM, double> normal_b = vd.getProp<normal_vector>(b);
                            // double r = sqrt(r2);
                            // double W = Wab(r, local_H, local_Kquintic);
                            // Point<DIM, double> dW = DWab(dr, r, local_H, local_Kquintic);

                            // double local_k = dotProduct(normal_b - normal_a, dW); // divergence of the normal vector
                            // if (local_k > max_curvature)
                            // {
                            //     local_k = max_curvature;
                            // }

                            // K_sum += local_k;
                            // w_sum += W;

                            // NEW CALCULATION
                            Point<DIM, double> normal_b = vd.getProp<normal_vector>(b);
                            double r = sqrt(r2);
                            double W = Wab(r, local_H, local_Kquintic);
                            Point<DIM, double> eab = -1.0 * dr / r;
                            double local_k = dotProduct(normal_b - normal_a, eab) / r;

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

void CalcVolume(particles &vd, double dp)
{
    // This function computes the curvature of the boundary particles from the divergence of the normal vector

    auto part = vd.getDomainIterator();

    // For each particle ...
    while (part.isNext())
    {
        // Key of the particle a
        vect_dist_key_dx a = part.get();

        // if particle BOUNDARY
        if (vd.getProp<type>(a) == BOUNDARY)
        {
            double kappa = vd.getProp<curvature_boundary>(a);
            double dxwall = vd.getProp<arc_length>(a);

            for (int i = 0; i < 3; i++)
            {
                vd.template getProp<vd_volume>(a)[i] = std::max(0.0, 0.5 * (2.0 * dp + dp * dp * kappa - 2.0 * (i + 1.0) * dp * dp * kappa) * dxwall);
            }
        }
        ++part;
    }
}

template <typename CellList>
void CalcDensity(particles &vd, CellList &NN, const Parameters &params)
{
    // This function computes the density of particles from the summation of the kernel

    double dp = params.dp;
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

            // Get the position xb of the particle a
            Point<DIM, double> xa = vd.getPos(a);

            // initialize density sum
            double rho_sum = 0.0;
            // double w_sum = 0.0;
            auto Np = NN.getNNIterator(NN.getCell(vd.getPos(a)));

            // iterate the neighborhood particles
            while (Np.isNext() == true)
            {
                // Key of b particle
                unsigned long b = Np.get();

                // Get the position xb of the particle b
                const Point<DIM, double> xb = vd.getPos(b);

                // Get the vector pointing at xa from xb
                const Point<DIM, double> dr = xa - xb;

                // take the norm squared of this vector
                const double r2 = norm2(dr);
                const double r = sqrt(r2);

                // If the particles interact ...
                if (r < params.r_threshold)
                {

                    if (vd.getProp<type>(b) == FLUID)
                    {

                        // evaluate kernel
                        const double w = Wab(r, params.H, params.Kquintic);
                        // w_sum += w * dp * dp;

                        rho_sum += w * params.MassFluid;
                    }
                    else
                    {
                        if (params.BC_TYPE == NO_SLIP)
                        {
                            const double r = sqrt(r2);

                            // evaluate kernel
                            const double w = Wab(r, params.H, params.Kquintic);
                            // w_sum += w * dp * dp;
                            rho_sum += w * params.MassBound;
                        }
                        else if (params.BC_TYPE == NEW_NO_SLIP) // need to evaluate kernel at dummy particles
                        {
                            // get normal vector of b
                            const Point<DIM, double> normal = vd.getProp<normal_vector>(b);
                            // get volumes, and curvature of dummy particles
                            const Point<3, double> vol = vd.getProp<vd_volume>(b);
                            const double kappa = vd.getProp<curvature_boundary>(b);
                            // Apply offsets to dr to get 3 vectrors pointing to dummy particles
                            const std::array<Point<DIM, double>, 3> R_dummy = getBoundaryPositions(-1.0 * dr, normal, params.dp);

                            // distance to the marker particle
                            const double dist2marker = sqrt(r2);

                            // if (dotProduct(dr, normal) > 0.0)
                            // {

                            for (int i = 0; i < 3; i++)
                            {
                                double rmax = sqrt(3.0 * 3.0 - (0.5 + (double)i) * (0.5 + (double)i)) * dp;
                                double rmin = (3.0 - (0.5 + (double)i)) * dp;
                                double kappa_max = 1.0 / (3.0 * dp);

                                // kappa = 0 gets rmax, kappa = kappa_max gets rmin
                                double r_interp = (rmin - rmax) / kappa_max * kappa + rmax;

                                if (dist2marker < r_interp)
                                {
                                    const double W = Wab(getVectorNorm(R_dummy[i]), params.H, params.Kquintic);
                                    rho_sum += W * vol[i] * params.rho_zero; // W*mass
                                }
                            }
                            // }
                        }
                    }
                }

                ++Np;
            }
            if (rho_sum != 0.0)
            {
                vd.template getProp<rho>(a) = rho_sum;
            }
            else
            {
                std::cout << "WARNING: NO particles around, density summation zero" << std::endl;
                vd.template getProp<rho>(a) = params.rho_zero;
            }
        }

        ++part;
    }
}

template <typename CellList>
void ExtrapolateVelocity(particles &vd, CellList &NN, const Parameters &params)
{
    // This function fills the value of v_transport for the boundary particles, with the velocity for no slip BC

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
            Point<DIM, double> xb = vd.getPos(b);

            // initialize sums
            Point<DIM, double> sum_vW = (DIM == 2) ? Point<DIM, double>{0.0, 0.0} : Point<DIM, double>{0.0, 0.0, 0.0};

            double sum_pW = 0.0;
            double sum_W = 0.0;

            auto Np = NN.getNNIterator(NN.getCell(vd.getPos(b)));

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
                Point<DIM, double> xf = vd.getPos(f);

                // Get the velocity of the fluid particle
                Point<DIM, double> vf = vd.getProp<velocity>(f);

                // Get the density of the fluid particle
                double rhof = vd.getProp<rho>(f);

                // Get the pressure of the fluid particle
                double Pf = vd.getProp<pressure>(f);

                // Get the vector pointing at xb from xf rwf
                Point<DIM, double> dr = xb - xf;

                // take the norm squared of this vector
                double r2 = norm2(dr);

                // If the particles interact ...
                if (r2 < params.r_threshold * params.r_threshold)
                {
                    // calculate distance
                    double r = sqrt(r2);

                    // evaluate kernel
                    double w = Wab(r, params.H, params.Kquintic);
                    // compute v*W
                    sum_vW += w * vf;

                    // compute Pf +rhof*(g-a) dot rwf
                    // at the moment a = 0
                    const double dot = dotProduct(dr, params.gravity_vector);
                    sum_pW += w * (Pf + rhof * dot);
                    sum_W += w;
                }

                ++Np;
            }
            if (sum_W != 0.0)
            {
                // Set the velocity of the boundary particle b ( to be used in the momentum equation to impose BC)
                // We use the v_transport field because boundary particles dont use this array

                // marker particles store centre of solid body in force_transport since it is unused
                // vector pointing from centre of rotation to marker particle
                const Point<DIM, double> radial_vec = {xb.get(0) - vd.getProp<force_transport>(b)[0], xb.get(1) - vd.getProp<force_transport>(b)[1]};
                const double radial = getVectorNorm(radial_vec);
                // vector tangential to the radial vector, rotation velocity is in this direction
                const Point<DIM, double> rotation_tangential = getPerpendicularUnit2D(radial_vec);
                Point<DIM, double> vw = vd.template getProp<velocity>(b);
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
