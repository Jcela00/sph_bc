#ifndef CALCULATIONS_H
#define CALCULATIONS_H

#include "Definitions.hpp"
#include "VectorUtilities.hpp"

template <typename CellList>
void CalcFluidVec2D(particles &vd, CellList &NN, const Parameters &params)
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
        if (vd.getProp<vd0_type>(a) != FLUID)
        {

            // check if the normal vector is uninitialized
            Point<DIM, real_number> normal = vd.template getProp<vd8_normal>(a);
            if (normal.get(0) == 0.0 && normal.get(1) == 0.0)
            {
                // Get the position xa of the particle a
                Point<DIM, real_number> xa = vd.getPos(a);

                // initialize sum
                Point<DIM, real_number> r_fluid_sum = {0.0};

                // neighborhood particles
                auto Np = NN.getNNIteratorBox(NN.getCell(vd.getPos(a)));

                // iterate the neighborhood particles
                while (Np.isNext() == true)
                {
                    // Key of b particle
                    const unsigned long b = Np.get();

                    if (vd.getProp<vd0_type>(b) == FLUID)
                    {
                        // Get the position xb of the particle b
                        const Point<DIM, real_number> xb = vd.getPos(b);

                        // Get the vector pointing at fluid from boundary
                        Point<DIM, real_number> dr = xb - xa;

                        // take the norm squared of this vector
                        const real_number r2 = norm2(dr);

                        // If the particles interact ...
                        if (r2 < params.r_cut * params.r_cut)
                        {
                            normalizeVector(dr);
                            r_fluid_sum += dr;
                        }
                    }

                    ++Np;
                }
                normalizeVector(r_fluid_sum);

                // store in normal temporally, this vector is only needed to have a reference when computing normals oriented "outside"
                // later is overwritten by the normal vector
                for (int xyz = 0; xyz < DIM; ++xyz)
                {
                    vd.template getProp<vd8_normal>(a)[xyz] = r_fluid_sum.get(xyz);
                }
            }
        }

        ++part;
    }
}

template <typename CellList>
void CalcFluidVec3D(particles &vd, CellList &NN, const Parameters &params)
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
        if (vd.getProp<vd0_type>(a) != FLUID)
        {

            // check if the normal vector is uninitialized
            Point<DIM, real_number> normal = vd.template getProp<vd8_normal>(a);
            if (normal.get(0) == 0.0 && normal.get(1) == 0.0 && normal.get(2) == 0.0)
            {
                // Get the position xa of the particle a
                Point<DIM, real_number> xa = vd.getPos(a);

                // initialize sum
                Point<DIM, real_number> r_fluid_sum = {0.0};

                // neighborhood particles
                auto Np = NN.getNNIteratorBox(NN.getCell(vd.getPos(a)));

                // iterate the neighborhood particles
                while (Np.isNext() == true)
                {
                    // Key of b particle
                    const unsigned long b = Np.get();

                    if (vd.getProp<vd0_type>(b) == FLUID)
                    {
                        // Get the position xb of the particle b
                        const Point<DIM, real_number> xb = vd.getPos(b);

                        // Get the vector pointing at fluid from boundary
                        Point<DIM, real_number> dr = xb - xa;

                        // take the norm squared of this vector
                        const real_number r2 = norm2(dr);

                        // If the particles interact ...
                        if (r2 < params.r_cut * params.r_cut)
                        {
                            normalizeVector(dr);
                            r_fluid_sum += dr;
                        }
                    }

                    ++Np;
                }
                normalizeVector(r_fluid_sum);

                // store in normal temporally, this vector is only needed to have a reference when computing normals oriented "outside"
                // later is overwritten by the normal vector
                for (int xyz = 0; xyz < DIM; ++xyz)
                {
                    vd.template getProp<vd8_normal>(a)[xyz] = r_fluid_sum.get(xyz);
                }
            }
        }

        ++part;
    }
}

template <typename CellList>
void CalcFluidVec(particles &vd, CellList &NN, const Parameters &params)
{
    if constexpr (DIM == 2)
    {
        CalcFluidVec2D(vd, NN, params);
    }
    else if constexpr (DIM == 3)
    {
        CalcFluidVec3D(vd, NN, params);
    }
}

template <typename CellList>
void CalcNormalVec2D(particles &vd, CellList &NN, const Parameters &params)
{
    // This function computes the normal vector for a boundary particle based on the other boundary particles inside its kernel.
    // it computes the average (kernel weighted) of the vectors perpendicular to the vectors pointing to the other boundary particles ( oriented to the fluid )

    const real_number local_H = params.H;
    const real_number local_Kquintic = params.Kquintic;

    auto part = vd.getDomainIterator();

    // Update the cell-list
    vd.updateCellList(NN);

    // For each particle ...
    while (part.isNext())
    {
        // Key of the particle a
        vect_dist_key_dx a = part.get();

        // if particle BOUNDARY
        if (vd.getProp<vd0_type>(a) != FLUID)
        {

            // get vector that points at fluid
            Point<DIM, real_number> Rfluid = vd.getProp<vd8_normal>(a);
            real_number normRfluid = getVectorNorm(Rfluid);

            // When we initialize the normals manually we give them a norm > 2.0 to distinguish the manually initialized from the ones we want to automatically initialize
            if (normRfluid < 2.0) // AUTOMATIC INITIALIZATION
            {
                // Get the position xa of the particle a
                Point<DIM, real_number> xa = vd.getPos(a);

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

                    if (vd.getProp<vd0_type>(b) != FLUID)
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
                    vd.template getProp<vd8_normal>(a)[xyz] = n_sum.get(xyz);
                }
            }
            else // Manually initialized, we just need to normalize it
            {
                normalizeVector(Rfluid);

                // store in normal vector
                for (int xyz = 0; xyz < DIM; ++xyz)
                {
                    vd.template getProp<vd8_normal>(a)[xyz] = Rfluid.get(xyz);
                }
            }
        }

        ++part;
    }
}

template <typename CellList>
void CalcNormalVec3D(particles &vd, CellList &NN, const Parameters &params)
{
    // This function computes the normal vector for a boundary particle based on the other boundary particles inside its kernel.
    // it computes the average (kernel weighted) of the vectors perpendicular to the vectors pointing to the other boundary particles ( oriented to the fluid )

    const real_number local_H = params.H;
    const real_number local_Kquintic = params.Kquintic;

    auto part = vd.getDomainIterator();

    // Update the cell-list
    vd.updateCellList(NN);

    // For each particle ...
    while (part.isNext())
    {
        // Key of the particle a
        vect_dist_key_dx a = part.get();

        // if particle BOUNDARY
        if (vd.getProp<vd0_type>(a) != FLUID)
        {

            // get vector that points at fluid
            Point<DIM, real_number> Rfluid = vd.getProp<vd8_normal>(a);
            real_number normRfluid = getVectorNorm(Rfluid);

            // When we initialize the normals manually we give them a norm > 2.0 to distinguish the manually initialized from the ones we want to automatically initialize
            if (normRfluid < 2.0) // AUTOMATIC INITIALIZATION
            {
                // Get the position xa of the particle a
                Point<DIM, real_number> xa = vd.getPos(a);

                // initialize sum
                Point<DIM, real_number> n_sum = {0.0};

                // neighborhood particles
                auto Np = NN.getNNIteratorBox(NN.getCell(xa));

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

                    if (vd.getProp<vd0_type>(b) != FLUID)
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
                            // // get perpendicular vector to dr
                            // Point<DIM, real_number> perp = getPerpendicularUnit2D(dr); // this is normalized

                            // // get scalar product of perp and Rfluid
                            // const real_number perp_dot_fluid = dotProduct(perp, Rfluid);

                            // // we want perp to point towards the fluid
                            // if (perp_dot_fluid < 0.0)
                            //     perp = -1.0 * perp;

                            // project rfluid onto the plane perpendicular to dr
                            Point<DIM, real_number> proj = Rfluid - (dotProduct(Rfluid, dr) / r2) * dr;
                            // normalize this projection
                            normalizeVector(proj);

                            // evaluate kernel
                            real_number W = Wab(sqrt(r2), local_H, local_Kquintic);
                            n_sum += proj * W;
                        }
                    }

                    ++Np;
                }
                // normalize the summed vector
                normalizeVector(n_sum);

                // store in normal vector
                for (int xyz = 0; xyz < DIM; ++xyz)
                {
                    vd.template getProp<vd8_normal>(a)[xyz] = n_sum.get(xyz);
                }
            }
            else // Manually initialized, we just need to normalize it
            {
                normalizeVector(Rfluid);

                // store in normal vector
                for (int xyz = 0; xyz < DIM; ++xyz)
                {
                    vd.template getProp<vd8_normal>(a)[xyz] = Rfluid.get(xyz);
                }
            }
        }

        ++part;
    }
}

template <typename CellList>
void CalcNormalVec(particles &vd, CellList &NN, const Parameters &params)
{
    if constexpr (DIM == 2)
    {
        CalcNormalVec2D(vd, NN, params);
    }
    else if constexpr (DIM == 3)
    {
        CalcNormalVec3D(vd, NN, params);
    }
}

template <typename CellList>
void CalcCurvature(particles &vd, CellList &NN, const Parameters &params)
{
    // This function computes the curvature of the boundary particles from the divergence of the normal vector

    const real_number dp = params.dp;
    const real_number local_H = params.H;
    const real_number local_Kquintic = params.Kquintic;

    auto part = vd.getDomainIterator();

    // Update the cell-list
    vd.updateCellList(NN);

    // max curvature is determined form the derivation of the volume formula, for curvature higher than this the volume of the third particle
    // is no longer between the two circles. Actually for 1/(2.5*dp) it becomes 0 and later negative

    const real_number max_curvature = 1.0 / (1.5 * dp);
    std::cout << "Max curvature: " << max_curvature << std::endl;

    // For each particle ...
    while (part.isNext())
    {
        // Key of the particle a
        vect_dist_key_dx a = part.get();

        // if particle BOUNDARY
        if (vd.getProp<vd0_type>(a) != FLUID)
        {

            Point<3, real_number> vol = vd.template getProp<vd9_volume>(a);
            if (vol[1] == 0.0) // if curvature is uninitialized
            {
                // Get the position xa of the particle a
                Point<DIM, real_number> xa = vd.getPos(a);

                // get normal of a
                Point<DIM, real_number> normal_a = vd.getProp<vd8_normal>(a);

                // initialize sums
                real_number K_sum = 0.0;
                real_number w_sum = 0.0;

                auto Np = NN.getNNIteratorBox(NN.getCell(vd.getPos(a)));

                // iterate the neighborhood particles
                while (Np.isNext() == true)
                {
                    // Key of b particle
                    const unsigned long b = Np.get();

                    if (vd.getProp<vd0_type>(b) != FLUID)
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
                                // Point<DIM, real_number> normal_b = vd.getProp<vd8_normal>(b);
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
                                Point<DIM, real_number> normal_b = vd.getProp<vd8_normal>(b);
                                real_number r = sqrt(r2);
                                real_number W = Wab(r, local_H, local_Kquintic);
                                Point<DIM, real_number> eab = dr / r;
                                real_number local_k = dotProduct(normal_a - normal_b, eab) / r;

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
                vd.template getProp<vd9_volume>(a)[1] = K_sum;
            }
            else // curvature was manually initialized but we want to set it to 0
            {
                vd.template getProp<vd9_volume>(a)[1] = 0.0;
            }
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
        if (vd.getProp<vd0_type>(a) != FLUID)
        {
            real_number dxwall = vd.getProp<vd9_volume>(a)[0];
            real_number kappa = vd.getProp<vd9_volume>(a)[1];

            for (int i = 0; i < 3; i++)
            {

                // if constexpr (DIM == 2)
                //     vd.template getProp<vd9_volume>(a)[i] = dp * dp;
                // else if constexpr (DIM == 3)
                //     vd.template getProp<vd9_volume>(a)[i] = dp * dp * dp;

                // n=i+1
                real_number n = static_cast<real_number>(i + 1);
                if constexpr (DIM == 2)
                    vd.template getProp<vd9_volume>(a)[i] = std::max(0.0, 0.5 * (2.0 * dp + dp * dp * kappa - 2.0 * n * dp * dp * kappa) * dxwall);
                else if constexpr (DIM == 3)
                    vd.template getProp<vd9_volume>(a)[i] = std::max(0.0, (1.0 / 3.0) * (3.0 * dp + kappa * dp * dp * (3.0 - 6.0 * n) + kappa * kappa * dp * dp * dp * (3 * n * n - 3 * n + 1)) * dxwall);
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
        if (vd.getProp<vd0_type>(a) == FLUID)
        {
            // Get the position xa of the particle a
            Point<DIM, real_number> xa = vd.getPos(a);

            // Get the velocity of the particle a
            Point<DIM, real_number> va = vd.getProp<vd4_velocity>(a);

            // initialize vorticity sum
            real_number vorticity = 0.0;

            // neighborhood particles
            auto Np = NN.getNNIteratorBox(NN.getCell(vd.getPos(a)));

            real_number rhoa = vd.getProp<vd1_rho>(a);

            // iterate the neighborhood particles
            while (Np.isNext() == true)
            {
                // Key of b particle
                const unsigned long b = Np.get();

                // Only consider other fluid particles
                if (vd.getProp<vd0_type>(b) == FLUID)
                {
                    // Get the position ( and vel ) of particle b
                    const Point<DIM, real_number> xb = vd.getPos(b);

                    const real_number rhob = vd.getProp<vd1_rho>(b);
                    // Get the vector pointing at A from B
                    Point<DIM, real_number> dr = xa - xb;
                    // get the norm squared of this vector
                    const real_number r2 = norm2(dr);

                    if (r2 < params.r_cut * params.r_cut)
                    {
                        // Get the velocity of particle b
                        const Point<DIM, real_number> vb = vd.getProp<vd4_velocity>(b);
                        // evaluate the kernel gradient
                        Point<DIM, real_number> Dwab = DWab(dr, sqrt(r2), params.H, params.Kquintic);
                        Point<DIM, real_number> dv = vb / rhob - va / rhoa;
                        vorticity += -params.MassFluid * crossProduct2D(dv, Dwab);
                    }
                }

                ++Np;
            }

            // Store in vorticity property
            // // vd.template getProp<vd11_vorticity>(a) = vorticity / rhoa;
            vd.template getProp<vd11_vorticity>(a) = vorticity;
        }

        ++part;
    }
}

#endif // CALCULATIONS_H
