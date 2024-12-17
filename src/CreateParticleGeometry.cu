#include "CreateParticleGeometry.hpp"

void CreateParticleGeometry(particles &vd, std::vector<std::pair<probe_particles, int>> &vp_vec, Obstacle *&obstacle_ptr, Parameters &params, AuxiliarParameters &auxParams)
{

    Vcluster<> &v_cl = create_vcluster();

    // Non periodic situation grid of 5 fluid particles and 3 boundary particles
    // We need a virtual grid of 5 + 2*(3+1) particles,
    // therefore the domain is discretized with 13 grid points,
    // when we use DrawParticles::DrawBox we will draw only the particles at the grid positons strictly inside the box,
    // the () repesent the recipient box, and the || represent the fluid box, we can see how this distribution places exactly 5 fluid particles inside and 3 boundary particles
    //           D-(-o--o--o-|-x--x--x--x--x--|-o-o-o-)-D
    // D: domain, o: boundary, x: fluid, --: dp distance
    // in a periodic situation we have the following
    // .....--x--x--D-|-x--x--x--x--x--|-D--x--x--......
    // therefore we need a grid of 5 + 2 particles, and the domain is discretized with 7 grid points

    // Size of the virtual cartesian grid that defines where to place the particles
    size_t sz[DIM];

    obstacle_ptr = new EmptyObstacle(params);
    real_number dp = params.dp;
    // Initialize obstacle in scenarios where needed
    if (params.SCENARIO == CYLINDER_ARRAY)
        obstacle_ptr = new CylinderObstacle(params.LengthScale, params.ObstacleCenter, params, params.ObstacleVelocity, params.ObstacleOmega, params.rf);
    else if (params.SCENARIO == CYLINDER_LATTICE)
        obstacle_ptr = new CylinderObstacle(params.LengthScale, params.ObstacleCenter, params, params.ObstacleVelocity, params.ObstacleOmega, params.rf);
    else if (params.SCENARIO == SQUARE)
        obstacle_ptr = new RectangleObstacle(params.ObstacleCenter, params, params.ObstacleBase, params.ObstacleHeight, params.ObstacleVelocity, params.ObstacleOmega, params.rf);
    else if (params.SCENARIO == TRIANGLE)
        obstacle_ptr = new TriangleObstacle(params.ObstacleCenter, params, params.ObstacleBase, params.ObstacleHeight, params.ObstacleVelocity, params.ObstacleOmega, params.rf);
    else if (params.SCENARIO == TRIANGLE_TEST)
        obstacle_ptr = new TriangleTestObstacle(params.ObstacleCenter, params, params.ObstacleBase, params.ObstacleHeight, params.ObstacleVelocity, params.ObstacleOmega, params.rf);
    else if (params.SCENARIO == TRIANGLE_EQUILATERAL)
        obstacle_ptr = new TriangleEqui(params.ObstacleCenter, params, params.ObstacleBase, params.ObstacleVelocity, params.ObstacleOmega, params.rf);
    else if (params.SCENARIO == ELLIPSE)
        obstacle_ptr = new EllipticObstacle(params.ObstacleBase, params.ObstacleHeight, params.ObstacleTilt, params.ObstacleCenter, params, params.ObstacleVelocity, params.ObstacleOmega, params.rf);
    else if (params.SCENARIO == FLOWER)
    {
        // obstacle_ptr = new CurveObstacle(params.flowerA, params.flowerB, params.flowerK, params.flowerM, params.ObstacleCenter, params, params.ObstacleVelocity, params.ObstacleOmega, params.rf);
        obstacle_ptr = new EpiCycloid(params.flowerA, params.flowerK, params.ObstacleCenter, params, params.ObstacleVelocity, params.ObstacleOmega, params.rf);
        // obstacle_ptr = new HipoCycloid(params.flowerA, params.flowerK, params.ObstacleCenter, params, params.ObstacleVelocity, params.ObstacleOmega, params.rf);
    }
    else if (params.SCENARIO == CUSTOM)
    {
        if (params.CustomObstacle == CYLINDER_LATTICE)
            obstacle_ptr = new CylinderObstacle(params.ObstacleBase, params.ObstacleCenter, params, params.ObstacleVelocity, params.ObstacleOmega, params.rf);
        else if (params.CustomObstacle == SQUARE)
            obstacle_ptr = new RectangleObstacle(params.ObstacleCenter, params, params.ObstacleBase, params.ObstacleHeight, params.ObstacleVelocity, params.ObstacleOmega, params.rf);
        else if (params.CustomObstacle == TRIANGLE)
            obstacle_ptr = new TriangleObstacle(params.ObstacleCenter, params, params.ObstacleBase, params.ObstacleHeight, params.ObstacleVelocity, params.ObstacleOmega, params.rf);
        else if (params.CustomObstacle == TRIANGLE_EQUILATERAL)
            obstacle_ptr = new TriangleEqui(params.ObstacleCenter, params, params.ObstacleBase, params.ObstacleVelocity, params.ObstacleOmega, params.rf);
        else if (params.CustomObstacle == ELLIPSE)
            obstacle_ptr = new EllipticObstacle(params.ObstacleBase, params.ObstacleHeight, params.ObstacleTilt, params.ObstacleCenter, params, params.ObstacleVelocity, params.ObstacleOmega, params.rf);
    }

    real_number refine_factor = params.rf;

    // Now define the iterator boxes
    // We define the boxes in terms of offstes with respect to the fluid box that goes from 0 to length
    real_number offset_domain_left[DIM] = {0.0};
    real_number offset_domain_right[DIM] = {0.0};
    real_number offset_recipient[DIM] = {0.0};
    real_number offset_periodic_fluid[DIM] = {0.0};
    real_number offset_periodic_recipient[DIM] = {0.0};

    for (int xyz = 0; xyz < DIM; xyz++)
    {
        if (params.bc[xyz] == NON_PERIODIC) // non periodic, fluid covered by boundary
        {
            sz[xyz] = params.Nfluid[xyz] + 2 * (params.Nboundary[xyz] + 1);
            offset_domain_left[xyz] = (0.5 + params.Nboundary[xyz]) * dp;
            offset_domain_right[xyz] = (0.5 + params.Nboundary[xyz]) * dp;

            if (params.BC_TYPE == NEW_NO_SLIP) // Nboundary should only be 0 or 1 if we are using the new bc
                offset_recipient[xyz] = 0.1 * params.Nboundary[xyz] * dp;
            else if (params.BC_TYPE == NO_SLIP)
                offset_recipient[xyz] = params.Nboundary[xyz] * dp;
        }
        else // periodic, open ended
        {
            sz[xyz] = params.Nfluid[xyz] + 1;

            offset_domain_left[xyz] = 0.0;
            offset_domain_right[xyz] = dp;
            offset_periodic_fluid[xyz] = 0.75 * dp;
            offset_periodic_recipient[xyz] = 0.85 * dp;
        }
    }

    // Define the boxes
    Box<DIM, real_number> domain({-offset_domain_left[0],
                                  -offset_domain_left[1]},
                                 {params.length[0] + offset_domain_right[0],
                                  params.length[1] + offset_domain_right[1]});

    Box<DIM, real_number> fluid_box({0.0,
                                     0.0},
                                    {params.length[0] + offset_periodic_fluid[0],
                                     params.length[1] + offset_periodic_fluid[1]});

    Box<DIM, real_number> recipient({-offset_recipient[0],
                                     -offset_recipient[1]},
                                    {params.length[0] + offset_recipient[0] + offset_periodic_recipient[0],
                                     params.length[1] + offset_recipient[1] + offset_periodic_recipient[1]});

    // Will only be used in the new bc
    // Box<DIM, real_number> recipient_hole({offset_recipient[0],
    //                                       offset_recipient[1]},
    //                                      {params.length[0] - offset_recipient[0] + offset_periodic_fluid[0],
    //                                       params.length[1] - offset_recipient[1] + offset_periodic_fluid[1]});

    for (int xyz = 0; xyz < DIM; xyz++) // correct length in periodic case
    {
        if (params.bc[xyz] == PERIODIC)
            params.length[xyz] += dp;
    }

    // extended boundary around the domain, and the processor domain
    Ghost<DIM, real_number> g(params.r_cut);

    // create particle object
    particles vd_loc(0, domain, params.bc, g, DEC_GRAN(128));
    // vd is argument passed as reference we want to fill with particles
    vd = vd_loc;

    // place probes
    if (params.PROBES_ENABLED)
    {
        // In practice this is only used in the Cylinder Lattice scenario
        // we want to place probes  in a vertical line at this locations

        Point<DIM, real_number> EndChannel = {params.length[0], 0.0};
        Point<DIM, real_number> HalfChannel = {params.length[0] / 2.0, 0.0};
        Point<DIM, real_number> VerticalOffset = {0.0, dp};
        Point<DIM, real_number> HorizontalOffset = {dp, 0.0};
        int k0 = 0;
        int kendHeight = params.Nfluid[1] + 1;

        std::vector<Point<DIM, real_number>> ProbePoints; // start points for the PlaceProbes function
        std::vector<int> ProbeComponents;                 // velocity component to measure 0 for x 1 for y
        std::vector<Point<DIM, real_number>> Offsets;
        std::vector<int> maxIters;

        if (params.SCENARIO == CYLINDER_LATTICE)
        {
            ProbePoints.push_back(HalfChannel);
            ProbePoints.push_back(EndChannel);

            ProbeComponents.push_back(0); // measure x velocity
            ProbeComponents.push_back(0); // measure x velocity

            Offsets.push_back(VerticalOffset);
            Offsets.push_back(VerticalOffset);

            maxIters.push_back(kendHeight);
            maxIters.push_back(kendHeight);
        }

        for (unsigned int k = 0; k < ProbePoints.size(); k++) // for each probe in the list
        {
            // create probe object
            Ghost<DIM, real_number> gp(0);
            size_t bc_p[DIM] = {NON_PERIODIC, NON_PERIODIC};
            probe_particles vp_loc(0, domain, bc_p, gp, DEC_GRAN(512));
            if (ProbeComponents[k] == 0)
            {
                openfpm::vector<std::string> names_p = {"vx"};
                vp_loc.setPropNames(names_p);
            }
            else if (ProbeComponents[k] == 1)
            {
                openfpm::vector<std::string> names_p = {"vy"};
                vp_loc.setPropNames(names_p);
            }

            if (v_cl.getProcessUnitID() == 0)
            {
                PlaceProbes(vp_loc, k0, maxIters[k], ProbePoints[k], Offsets[k], obstacle_ptr);
            }
            std::pair<probe_particles, int> tmp = std::make_pair(vp_loc, ProbeComponents[k]);
            vp_vec.push_back(tmp);
            auxParams.probe_filenames.push_back("probes_" + std::to_string(k) + "_" + auxParams.filename);
        }
    }

    // For the new bc we place the wall/obstacle particles manually
    // Add the obstacle/walls as marker particles only on processor 0
    if (params.BC_TYPE == NEW_NO_SLIP && v_cl.getProcessUnitID() == 0)
    {
        // Add obstacle
        obstacle_ptr->AddObstacle(vd);

        // Add walls
        if (params.bc[0] == PERIODIC && params.bc[1] == NON_PERIODIC) // Channel like scenario
        {
            real_number dx_wall = dp / refine_factor;
            int Nwall = ceil(params.length[0] / dx_wall);
            dx_wall = params.length[0] / Nwall;
            Point<DIM, real_number> X_Offset = {dx_wall, 0.0};

            Point<DIM, real_number> LL_corner = {0.0, 0.0};
            Point<DIM, real_number> UL_corner = {0.0, params.length[1]};
            // Top And Bottom Walls

            if (params.SCENARIO == CUSTOM)
            {
                AddFlatWallNewBC(vd, 0, Nwall, LL_corner, X_Offset, dx_wall, {0.0, 0.0}, params.vw_bottom, params, FREE_SLIP_BOUNDARY, 0.0);
                AddFlatWallNewBC(vd, 0, Nwall, UL_corner, X_Offset, dx_wall, {0.0, 0.0}, params.vw_top, params, FREE_SLIP_BOUNDARY, 0.0);
            }
            else
            {
                AddFlatWallNewBC(vd, 0, Nwall, LL_corner, X_Offset, dx_wall, {0.0, 0.0}, params.vw_bottom, params, BOUNDARY, 0.0);
                AddFlatWallNewBC(vd, 0, Nwall, UL_corner, X_Offset, dx_wall, {0.0, 0.0}, params.vw_top, params, BOUNDARY, 0.0);
            }
        }
        else if (params.bc[0] == NON_PERIODIC && params.bc[1] == NON_PERIODIC) // Box like scenario
        {
            real_number dx_wall_x = dp / refine_factor;
            int Nwall_x = ceil(params.length[0] / dx_wall_x);
            dx_wall_x = params.length[0] / Nwall_x;
            Point<DIM, real_number> X_Offset = {dx_wall_x, 0.0};

            real_number dx_wall_y = dp / refine_factor;
            int Nwall_y = ceil(params.length[1] / dx_wall_y);
            dx_wall_y = params.length[1] / Nwall_y;
            Point<DIM, real_number> Y_Offset = {0.0, dx_wall_y};

            Point<DIM, real_number> LL_corner = {0.0, 0.0};
            Point<DIM, real_number> LR_corner = {params.length[0], 0.0};

            Point<DIM, real_number> UL_corner = {0.0, params.length[1]};

            // Top And Bottom Walls
            AddFlatWallNewBC(vd, 0, Nwall_x + 1, LL_corner, X_Offset, dx_wall_x, {0.0, 0.0}, params.vw_bottom, params, BOUNDARY, 0.0);
            AddFlatWallNewBC(vd, 0, Nwall_x + 1, UL_corner, X_Offset, dx_wall_x, {0.0, 0.0}, params.vw_top, params, BOUNDARY, 0.0);

            // Left And Right Walls
            AddFlatWallNewBC(vd, 1, Nwall_y, LL_corner, Y_Offset, dx_wall_y, {0.0, 0.0}, {0.0, 0.0}, params, BOUNDARY, 0.0);
            AddFlatWallNewBC(vd, 1, Nwall_y, LR_corner, Y_Offset, dx_wall_y, {0.0, 0.0}, {0.0, 0.0}, params, BOUNDARY, 0.0);
        }
    }

    // Now place fluid particles
    // For the old bc, some of the fluid positions will be inside the obstacle
    // if that is true we need to set them as boundary particles
    // For the new bc, we will not add particles that are inside the obstacle, since we already added them

    // return an iterator to the fluid particles to add to vd
    auto fluid_it = DrawParticles::DrawBox(vd, sz, domain, fluid_box);

    // for each particle inside the fluid box ...
    while (fluid_it.isNext())
    {

        Point<DIM, real_number> iterator_position = fluid_it.get();

        if ((*obstacle_ptr).isInside(iterator_position)) // if inside the obstacle region
        {
            if (params.BC_TYPE == NO_SLIP) // add particle but set it as boundary
            {
                // ... add a particle ...
                vd.add();
                vd.template getLastProp<vd0_type>() = OBSTACLE;
                vd.template getLastProp<vd10_omega>() = (*obstacle_ptr).AngularVelocity_;
                for (int xyz = 0; xyz < DIM; xyz++)
                {
                    vd.template getLastProp<vd4_velocity>()[xyz] = ((*obstacle_ptr).LinearVelocity_).get(xyz);
                    vd.template getLastProp<vd7_force_t>()[xyz] = ((*obstacle_ptr).Centre_).get(xyz);
                }
            }
            else if (params.BC_TYPE == NEW_NO_SLIP) // not add particle because already added
            {
                ++fluid_it;
                continue;
            }
        }
        else // if not inside obstacle at all just add fluid particles
        {
            // ... add a particle ...
            vd.add();
            vd.template getLastProp<vd0_type>() = FLUID;
            vd.template getLastProp<vd10_omega>() = 0.0;
            for (int xyz = 0; xyz < DIM; xyz++)
            {

                // if (params.SCENARIO == CUSTOM)
                // {
                //     vd.template getLastProp<vd4_velocity>()[xyz] = params.Vinflow[xyz];
                // }
                // else
                // {
                vd.template getLastProp<vd4_velocity>()[xyz] = 0.0;
                // }
                // else
                // {
                // vd.template getLastProp<vd4_velocity>()[xyz] = 0.0;
                // }

                vd.template getLastProp<vd7_force_t>()[xyz] = 0.0;
            }
        }

        // Set properties common to fluid and obstacle particles
        vd.template getLastProp<vd1_rho>() = params.rho0;
        vd.template getLastProp<vd2_pressure>() = 0.0;
        vd.template getLastProp<vd3_drho>() = 0.0;

        for (int xyz = 0; xyz < DIM; xyz++)
        {
            //  random displacement ( useful for debugging )
            // double random_number = static_cast<double>(std::rand()) / RAND_MAX * 2.0 - 1.0;

            vd.getLastPos()[xyz] = iterator_position.get(xyz); //+ 0.1 * dp * random_number;
            vd.template getLastProp<vd5_velocity_t>()[xyz] = 0.0;
            vd.template getLastProp<vd8_normal>()[xyz] = 0.0;
            vd.template getLastProp<vd6_force>()[xyz] = 0.0;
        }

        // fluid particles carry dp in volume[0], it is useful to read it in postprocessing
        // this array is unused for fluid particles anyway
        vd.template getLastProp<vd9_volume>()[0] = dp;
        vd.template getLastProp<vd9_volume>()[1] = 0.0;
        vd.template getLastProp<vd9_volume>()[2] = 0.0;

        vd.template getLastProp<vd11_vorticity>() = 0.0;
        vd.template getLastProp<vd12_vel_red>() = 0.0;
        vd.template getLastProp<vd13_force_red_x>() = 0.0;
        vd.template getLastProp<vd14_force_red_y>() = 0.0;

        // next fluid particle
        ++fluid_it;
    }

    // Now place solid walls using iterators (only for OLD BC)

    if (params.BC_TYPE == NO_SLIP)
    {

        openfpm::vector<Box<DIM, real_number>> holes;
        holes.add(fluid_box);
        Box<DIM, real_number> hole_box = holes.get(0);
        auto bound_box = DrawParticles::DrawSkin(vd, sz, domain, holes, recipient);

        if (params.bc[0] != PERIODIC || params.bc[1] != PERIODIC) // no walls in all periodic scenario
        {
            while (bound_box.isNext())
            {
                Point<DIM, real_number> position = bound_box.get();

                // periodic bc, with no boundary particles in y direction has a bug, it puts 3 extra particles outside in the y direction
                // When running on multiple cores, with this we check if particle is outside the recipient box
                // Another bug places boundary particles in the correct plane, but inside the fluid box;
                // if (bc[0] == PERIODIC && position.get(0) > dp / 2.0 && position.get(0) < length[0] - dp / 2.0)
                // {
                // 	++bound_box;
                // 	continue;
                // }

                if (!recipient.isInside((position)))
                {
                    ++bound_box;
                    continue;
                }
                if (hole_box.isInside(position))
                {
                    ++bound_box;
                    continue;
                }

                vd.add();

                vd.template getLastProp<vd0_type>() = BOUNDARY;
                vd.template getLastProp<vd1_rho>() = params.rho0;
                vd.template getLastProp<vd2_pressure>() = 0.0;
                vd.template getLastProp<vd3_drho>() = 0.0;

                for (int xyz = 0; xyz < DIM; xyz++)
                {
                    vd.getLastPos()[xyz] = bound_box.get().get(xyz);
                    vd.template getLastProp<vd6_force>()[xyz] = 0.0;
                    vd.template getLastProp<vd7_force_t>()[xyz] = 0.0;
                    vd.template getLastProp<vd5_velocity_t>()[xyz] = 0.0;
                    vd.template getLastProp<vd8_normal>()[xyz] = 0.0;
                    if (position.get(1) < dp / 4.0) // bottom wall
                    {
                        vd.template getLastProp<vd4_velocity>()[xyz] = params.vw_bottom[xyz];
                    }
                    else if (position.get(1) > params.length[1] - dp / 4.0) // top wall
                    {
                        vd.template getLastProp<vd4_velocity>()[xyz] = params.vw_top[xyz];
                    }
                }

                vd.template getLastProp<vd9_volume>()[0] = dp;
                vd.template getLastProp<vd9_volume>()[1] = 0.0;
                vd.template getLastProp<vd9_volume>()[2] = 0.0;
                vd.template getLastProp<vd11_vorticity>() = 0.0;
                vd.template getLastProp<vd12_vel_red>() = 0.0;
                vd.template getLastProp<vd13_force_red_x>() = 0.0;
                vd.template getLastProp<vd14_force_red_y>() = 0.0;

                ++bound_box;
            }
        }
    }
}

void CreateParticleGeometryPoiseuilleTank(particles &vd, Parameters &params, AuxiliarParameters &auxParams)
{

    Vcluster<> &v_cl = create_vcluster();

    // Non periodic situation grid of 5 fluid particles and 3 boundary particles
    // We need a virtual grid of 5 + 2*(3+1) particles,
    // therefore the domain is discretized with 13 grid points,
    // when we use DrawParticles::DrawBox we will draw only the particles at the grid positons strictly inside the box,
    // the () repesent the recipient box, and the || represent the fluid box, we can see how this distribution places exactly 5 fluid particles inside and 3 boundary particles
    //           D-(-o--o--o-|-x--x--x--x--x--|-o-o-o-)-D
    // D: domain, o: boundary, x: fluid, --: dp distance
    // in a periodic situation we have the following
    // .....--x--x--D-|-x--x--x--x--x--|-D--x--x--......
    // therefore we need a grid of 5 + 2 particles, and the domain is discretized with 7 grid points

    // Size of the virtual cartesian grid that defines where to place the particles
    size_t sz[DIM];

    real_number dp = params.dp;

    real_number refine_factor = params.rf;

    // Now define the iterator boxes
    // We define the boxes in terms of offstes with respect to the fluid box that goes from 0 to length
    real_number offset_domain_left[DIM] = {0.0};
    real_number offset_domain_right[DIM] = {0.0};
    real_number offset_recipient[DIM] = {0.0};
    real_number offset_periodic_fluid[DIM] = {0.0};
    real_number offset_periodic_recipient[DIM] = {0.0};

    for (int xyz = 0; xyz < DIM; xyz++)
    {
        if (params.bc[xyz] == NON_PERIODIC) // non periodic, fluid covered by boundary
        {
            sz[xyz] = params.Nfluid[xyz] + 2 * (params.Nboundary[xyz] + 1);
            offset_domain_left[xyz] = (0.5 + params.Nboundary[xyz]) * dp;
            offset_domain_right[xyz] = (0.5 + params.Nboundary[xyz]) * dp;

            if (params.BC_TYPE == NEW_NO_SLIP) // Nboundary should only be 0 or 1 if we are using the new bc
                offset_recipient[xyz] = 0.1 * params.Nboundary[xyz] * dp;
            else if (params.BC_TYPE == NO_SLIP)
                offset_recipient[xyz] = params.Nboundary[xyz] * dp;
        }
        else // periodic, open ended
        {
            sz[xyz] = params.Nfluid[xyz] + 1;

            offset_domain_left[xyz] = 0.0;
            offset_domain_right[xyz] = dp;
            offset_periodic_fluid[xyz] = 0.75 * dp;
            offset_periodic_recipient[xyz] = 0.85 * dp;
        }
    }

    // auxiliary tank length and heigth
    real_number tank_length = params.length[0] - 10.0 * dp;
    real_number tank_height = 10.0 * dp;

    sz[1] += 20;

    Box<DIM, real_number> second_fluid_box({5.0 * dp,
                                            params.length[1] + tank_height},
                                           {params.length[0] - 5.0 * dp,
                                            params.length[1] + 2.0 * tank_height});

    // Define the boxes
    Box<DIM, real_number> domain({-offset_domain_left[0],
                                  -offset_domain_left[1]},
                                 {params.length[0] + offset_domain_right[0],
                                  params.length[1] + offset_domain_right[1] + 2.0 * tank_height});

    Box<DIM, real_number> fluid_box({0.0,
                                     0.0},
                                    {params.length[0] + offset_periodic_fluid[0],
                                     params.length[1] + offset_periodic_fluid[1]});

    Box<DIM, real_number> recipient({-offset_recipient[0],
                                     -offset_recipient[1]},
                                    {params.length[0] + offset_recipient[0] + offset_periodic_recipient[0],
                                     params.length[1] + offset_recipient[1] + offset_periodic_recipient[1]});

    for (int xyz = 0; xyz < DIM; xyz++) // correct length in periodic case
    {
        if (params.bc[xyz] == PERIODIC)
            params.length[xyz] += dp;
    }

    // extended boundary around the domain, and the processor domain
    Ghost<DIM, real_number> g(params.r_cut);

    // create particle object
    particles vd_loc(0, domain, params.bc, g, DEC_GRAN(128));
    // vd is argument passed as reference we want to fill with particles
    vd = vd_loc;

    // For the new bc we place the wall/obstacle particles manually
    // Add the obstacle/walls as marker particles only on processor 0
    if (params.BC_TYPE == NEW_NO_SLIP && v_cl.getProcessUnitID() == 0)
    {
        // Add walls
        if (params.bc[0] == PERIODIC && params.bc[1] == NON_PERIODIC) // Channel like scenario
        {
            real_number dx_wall = dp / refine_factor;
            int Nwall = ceil(params.length[0] / dx_wall);
            dx_wall = params.length[0] / Nwall;
            Point<DIM, real_number> X_Offset = {dx_wall, 0.0};
            Point<DIM, real_number> Y_Offset = {0.0, dx_wall};

            Point<DIM, real_number> LL_corner = {0.0, 0.0};
            Point<DIM, real_number> UL_corner = {0.0, params.length[1]};

            // Top And Bottom Walls
            AddFlatWallNewBC(vd, 0, Nwall, LL_corner, X_Offset, dx_wall, {0.0, 0.0}, params.vw_bottom, params, BOUNDARY, 0.0);
            AddFlatWallNewBC(vd, 0, Nwall, UL_corner, X_Offset, dx_wall, {0.0, 0.0}, params.vw_top, params, BOUNDARY, 0.0);

            // Add auxiliary tank on top of the channel
            Point<DIM, real_number> tank_LL_corner = {0.0 + 5.0 * dp - 0.5 * dp, params.length[1] + tank_height};
            Point<DIM, real_number> tank_UR_corner = {params.length[0] - 5.0 * dp - 0.5 * dp, params.length[1] + 2.0 * tank_height};

            int Nwall_tank = ceil(tank_length / dx_wall);
            int Nwall_tank_y = ceil(tank_height / dx_wall);
            // Top And Bottom Walls
            AddFlatWallNewBC(vd, 0, Nwall_tank + 1, tank_LL_corner, X_Offset, dx_wall, {0.0, 0.0}, {0.0, 0.0}, params, BOUNDARY, 0.0);
            AddFlatWallNewBC(vd, 1, Nwall_tank_y + 1, tank_LL_corner, Y_Offset, dx_wall, {0.0, 0.0}, {0.0, 0.0}, params, BOUNDARY, 0.0);

            AddFlatWallNewBC(vd, 0, Nwall_tank + 1, tank_UR_corner, -X_Offset, dx_wall, {0.0, 0.0}, {0.0, 0.0}, params, BOUNDARY, 0.0);
            AddFlatWallNewBC(vd, 1, Nwall_tank_y + 1, tank_UR_corner, -Y_Offset, dx_wall, {0.0, 0.0}, {0.0, 0.0}, params, BOUNDARY, 0.0);
        }
    }

    // Now place fluid particles
    // For the old bc, some of the fluid positions will be inside the obstacle
    // if that is true we need to set them as boundary particles
    // For the new bc, we will not add particles that are inside the obstacle, since we already added them

    // return an iterator to the fluid particles to add to vd
    auto fluid_it = DrawParticles::DrawBox(vd, sz, domain, fluid_box);

    // for each particle inside the fluid box ...
    while (fluid_it.isNext())
    {

        Point<DIM, real_number> iterator_position = fluid_it.get();

        // ... add a particle ...
        vd.add();
        vd.template getLastProp<vd0_type>() = FLUID;
        vd.template getLastProp<vd10_omega>() = 0.0;
        for (int xyz = 0; xyz < DIM; xyz++)
        {
            vd.template getLastProp<vd4_velocity>()[xyz] = 0.0;
            vd.template getLastProp<vd7_force_t>()[xyz] = 0.0;
        }

        // Set properties common to fluid and obstacle particles
        vd.template getLastProp<vd1_rho>() = params.rho0;
        vd.template getLastProp<vd2_pressure>() = 0.0;
        vd.template getLastProp<vd3_drho>() = 0.0;

        for (int xyz = 0; xyz < DIM; xyz++)
        {
            vd.getLastPos()[xyz] = iterator_position.get(xyz);
            vd.template getLastProp<vd5_velocity_t>()[xyz] = 0.0;
            vd.template getLastProp<vd8_normal>()[xyz] = 0.0;
            vd.template getLastProp<vd6_force>()[xyz] = 0.0;
        }

        // fluid particles carry dp in volume[0], it is useful to read it in postprocessing
        // this array is unused for fluid particles anyway
        vd.template getLastProp<vd9_volume>()[0] = dp;
        vd.template getLastProp<vd9_volume>()[1] = 0.0;
        vd.template getLastProp<vd9_volume>()[2] = 0.0;

        vd.template getLastProp<vd11_vorticity>() = 0.0;
        vd.template getLastProp<vd12_vel_red>() = 0.0;
        vd.template getLastProp<vd13_force_red_x>() = 0.0;
        vd.template getLastProp<vd14_force_red_y>() = 0.0;

        // next fluid particle
        ++fluid_it;
    }

    auto fluid_it2 = DrawParticles::DrawBox(vd, sz, domain, second_fluid_box);
    // for each particle inside the fluid box ...
    while (fluid_it2.isNext())
    {

        Point<DIM, real_number> iterator_position = fluid_it2.get();

        // ... add a particle ...
        vd.add();
        vd.template getLastProp<vd0_type>() = FLUID;
        vd.template getLastProp<vd10_omega>() = 0.0;
        for (int xyz = 0; xyz < DIM; xyz++)
        {
            vd.template getLastProp<vd4_velocity>()[xyz] = 0.0;
            vd.template getLastProp<vd7_force_t>()[xyz] = 0.0;
        }

        // Set properties
        vd.template getLastProp<vd1_rho>() = params.rho0;
        vd.template getLastProp<vd2_pressure>() = 0.0;
        vd.template getLastProp<vd3_drho>() = 0.0;

        for (int xyz = 0; xyz < DIM; xyz++)
        {
            vd.getLastPos()[xyz] = iterator_position.get(xyz);
            vd.template getLastProp<vd5_velocity_t>()[xyz] = 0.0;
            vd.template getLastProp<vd8_normal>()[xyz] = 0.0;
        }

        vd.template getLastProp<vd9_volume>()[0] = dp;
        vd.template getLastProp<vd9_volume>()[1] = 0.0;
        vd.template getLastProp<vd9_volume>()[2] = 0.0;

        // next fluid particle
        ++fluid_it2;
    }

    // Create a Cell list object for CPU calculations
    auto NN = vd.getCellList(params.r_cut);
    vd.updateCellList(NN);

    CalcFluidVec(vd, NN, params);
    vd.ghost_get<vd0_type, vd1_rho, vd2_pressure, vd4_velocity, vd5_velocity_t, vd6_force, vd7_force_t, vd8_normal, vd9_volume, vd10_omega, vd11_vorticity>();

    CalcNormalVec(vd, NN, params);
    vd.ghost_get<vd0_type, vd1_rho, vd2_pressure, vd4_velocity, vd5_velocity_t, vd6_force, vd7_force_t, vd8_normal, vd9_volume, vd10_omega, vd11_vorticity>();

    CalcCurvature(vd, NN, params);
    vd.ghost_get<vd0_type, vd1_rho, vd2_pressure, vd4_velocity, vd5_velocity_t, vd6_force, vd7_force_t, vd8_normal, vd9_volume, vd10_omega, vd11_vorticity>();

    CalcVolume(vd, params.dp);
    vd.ghost_get<vd0_type, vd1_rho, vd2_pressure, vd4_velocity, vd5_velocity_t, vd6_force, vd7_force_t, vd8_normal, vd9_volume, vd10_omega, vd11_vorticity>();

    vd.write_frame("frametest", 0, 0.0, VTK_WRITER);
    // now move aux tank beside the dam break
    real_number ythreshold = params.length[1] + 0.5 * tank_height;
    auto it = vd.getDomainIterator();
    while (it.isNext())
    {
        auto akey = it.get();

        Point<DIM, real_number> position = vd.getPos(akey);
        if (position.get(1) > ythreshold)
        {
            vd.getPos(akey)[1] -= (tank_height - 0.5 * dp);
        }

        ++it;
    }
}

void CreateParticleGeometryTaylorCouette(particles &vd, std::vector<std::pair<probe_particles, int>> &vp_vec, Obstacle *&obstacle_ptr, Parameters params, AuxiliarParameters &auxParams)
{
    Vcluster<> &v_cl = create_vcluster();

    // Size of the virtual cartesian grid that defines where to place the particles
    size_t sz[DIM];

    real_number dp = params.dp;
    // We define the boxes in terms of offstes with respect to the fluid box that goes from 0 to length
    real_number offset_domain_left[DIM] = {0.0};
    real_number offset_domain_right[DIM] = {0.0};

    // non periodic situation grid of 5 fluid particles and 3 boundary particles
    // We need a virtual grid of 5 + 2*(3+1) particles,
    // therefore the domain is discretized with 13 grid points,
    // when we use DrawParticles::DrawBox we will draw only the particles at the grid positons strictly inside the box,
    // the () repesent the recipient box, and the || represent the fluid box, we can see how this distribution places exactly 5 fluid particles inside and 3 boundary particles
    //           D-(-o--o--o-|-x--x--x--x--x--|-o-o-o-)-D
    // D: domain, o: boundary, x: fluid, --: dp distance
    // in a periodic situation we have the following
    // .....--x--x--D-|-x--x--x--x--x--|-D--x--x--......
    // therefore we need a grid of 5 + 2 particles, and the domain is discretized with 7 grid points

    real_number Rin = params.Rin;
    real_number Rout = params.Rout;
    real_number Win = params.Win;
    real_number Wout = params.Wout;

    real_number a_tc = -((Rout * Rout * Rin * Rin) / (Rout * Rout - Rin * Rin)) * (Wout - Win);
    real_number b_tc = (Wout * Rout * Rout - Win * Rin * Rin) / (Rout * Rout - Rin * Rin);

    size_t Nbound = (params.BC_TYPE == NEW_NO_SLIP) ? 1 : 3;

    for (int dim = 0; dim < DIM; dim++)
    {
        params.length[dim] = dp * params.Nfluid[dim];
        sz[dim] = params.Nfluid[dim] + 2 * (Nbound + 1);
        offset_domain_left[dim] = (0.5 + Nbound) * dp;
        offset_domain_right[dim] = (0.5 + Nbound) * dp;
    }

    // Define the boxes
    Box<DIM, real_number> domain({-params.length[0] / 2.0 - offset_domain_left[0],
                                  -params.length[1] / 2.0 - offset_domain_left[1]},
                                 {params.length[0] / 2.0 + offset_domain_right[0],
                                  params.length[1] / 2.0 + offset_domain_right[1]});
    Box<DIM, real_number> fluid_box({-params.length[0] / 2.0,
                                     -params.length[1] / 2.0},
                                    {params.length[0] / 2.0,
                                     params.length[1] / 2.0});

    // extended boundary around the domain, and the processor domain
    Ghost<DIM, real_number> g(params.r_cut);

    // create particle object
    particles vd_loc(0, domain, params.bc, g, DEC_GRAN(128));
    vd = vd_loc;

    // Write constants on file
    real_number rf = params.rf;

    // Set cylindrical object parameters
    Point<DIM, real_number> CylinderCentre = {0.0, 0.0};

    obstacle_ptr = new EmptyObstacle(params);

    const Point<DIM, real_number> vel = {0.0, 0.0};

    CylinderObstacle *obstacle_ptr_out = new CylinderObstacle(Rout, CylinderCentre, params, vel, Wout, rf);
    CylinderObstacle *obstacle_ptr_in = new CylinderObstacle(Rin, CylinderCentre, params, vel, Win, rf);

    CylinderObstacle *obstacle_ptr_out_aux = new CylinderObstacle(Rout + 3.0 * dp, CylinderCentre, params, vel, Wout, rf);
    CylinderObstacle *obstacle_ptr_in_aux = new CylinderObstacle(Rin - 3.0 * dp, CylinderCentre, params, vel, Win, rf);

    // Add the obstacle as marker particles only on processor 0
    if (params.BC_TYPE == NEW_NO_SLIP && v_cl.getProcessUnitID() == 0)
    {
        obstacle_ptr_out->AddObstacle(vd);
        obstacle_ptr_in->AddObstacle(vd);
    }

    Box<DIM, real_number> fluid_box_aux({-Rout - 3.5 * dp,
                                         -Rout - 3.5 * dp},
                                        {Rout + 3.5 * dp,
                                         Rout + 3.5 * dp});

    // Outer Cylinder boundary particles
    auto out_it = DrawParticles::DrawBox(vd, sz, domain, fluid_box_aux);
    if (params.BC_TYPE == NO_SLIP)
    {
        while (out_it.isNext())
        {

            Point<DIM, real_number> iterator_position = out_it.get();
            if (!(*obstacle_ptr_out).isInside_minEps(iterator_position) && (*obstacle_ptr_out_aux).isInside(iterator_position)) // if outside the outer cylinder and inside outer cylinder aux
            {
                if (params.BC_TYPE == NO_SLIP)
                {
                    vd.add();
                    // Set properties
                    vd.template getLastProp<vd0_type>() = OBSTACLE;
                    vd.template getLastProp<vd10_omega>() = (*obstacle_ptr_out).AngularVelocity_;
                    for (int xyz = 0; xyz < DIM; xyz++)
                    {
                        vd.template getLastProp<vd4_velocity>()[xyz] = ((*obstacle_ptr_out).LinearVelocity_).get(xyz);
                        vd.template getLastProp<vd7_force_t>()[xyz] = ((*obstacle_ptr_out).Centre_).get(xyz);
                    }
                    vd.template getLastProp<vd1_rho>() = params.rho0;
                    vd.template getLastProp<vd2_pressure>() = 0.0;
                    vd.template getLastProp<vd3_drho>() = 0.0;

                    for (int xyz = 0; xyz < DIM; xyz++)
                    {
                        vd.getLastPos()[xyz] = iterator_position.get(xyz);
                        vd.template getLastProp<vd5_velocity_t>()[xyz] = 0.0;
                        vd.template getLastProp<vd8_normal>()[xyz] = 0.0;
                    }

                    vd.template getLastProp<vd9_volume>()[0] = dp;

                    // next fluid particle
                    ++out_it;
                    continue;
                }
                else
                {
                    ++out_it;
                    continue;
                }
            }
            else // skip fluid particle
            {
                ++out_it;
                continue;
            }
        }
    }

    // Inner Cylinder boundary particles and fluid particles
    auto fluid_it = DrawParticles::DrawBox(vd, sz, domain, fluid_box);
    while (fluid_it.isNext())
    {

        Point<DIM, real_number> iterator_position = fluid_it.get();
        if ((*obstacle_ptr_out).isInside_minEps(iterator_position)) // if inside the outer cylinder
        {
            if ((*obstacle_ptr_in).isInside(iterator_position)) // if inside the inner cylinder region
            {
                if (!(*obstacle_ptr_in_aux).isInside_minEps(iterator_position))
                {
                    if (params.BC_TYPE == NO_SLIP) // add particle but set it as boundary
                    {
                        // ... add a particle ...
                        vd.add();
                        vd.template getLastProp<vd0_type>() = OBSTACLE;
                        vd.template getLastProp<vd10_omega>() = (*obstacle_ptr_in).AngularVelocity_;
                        for (int xyz = 0; xyz < DIM; xyz++)
                        {
                            vd.template getLastProp<vd4_velocity>()[xyz] = ((*obstacle_ptr_in).LinearVelocity_).get(xyz);
                            vd.template getLastProp<vd7_force_t>()[xyz] = ((*obstacle_ptr_in).Centre_).get(xyz);
                        }
                    }
                    else if (params.BC_TYPE == NEW_NO_SLIP) // not add particle because already added
                    {
                        ++fluid_it;
                        continue;
                    }
                }
                else
                {
                    ++fluid_it;
                    continue;
                }
            }
            else // if no cylinder at all just add fluid particles
            {
                // ... add a particle ...
                vd.add();
                vd.template getLastProp<vd0_type>() = FLUID;
                vd.template getLastProp<vd10_omega>() = 0.0;

                // real_number r = iterator_position.get(0) * iterator_position.get(0) + iterator_position.get(1) * iterator_position.get(1);
                // r = sqrt(r);
                // real_number uth = a_tc / r + b_tc * r;

                // real_number ux = uth * (-iterator_position.get(1) / r);
                // real_number uy = uth * (iterator_position.get(0) / r);

                vd.template getLastProp<vd4_velocity>()[0] = 0.0; // ux;
                vd.template getLastProp<vd4_velocity>()[1] = 0.0; // uy;

                for (int xyz = 0; xyz < DIM; xyz++)
                {
                    vd.template getLastProp<vd7_force_t>()[xyz] = 0.0;
                }
            }
        }
        else // skip fluid particle
        {
            ++fluid_it;
            continue;
        }
        // Set properties
        vd.template getLastProp<vd1_rho>() = params.rho0;
        vd.template getLastProp<vd2_pressure>() = 0.0;
        vd.template getLastProp<vd3_drho>() = 0.0;

        for (int xyz = 0; xyz < DIM; xyz++)
        {
            vd.getLastPos()[xyz] = iterator_position.get(xyz);
            vd.template getLastProp<vd5_velocity_t>()[xyz] = 0.0;
            vd.template getLastProp<vd8_normal>()[xyz] = 0.0;
        }

        vd.template getLastProp<vd9_volume>()[0] = dp;

        // next fluid particle
        ++fluid_it;
    }
}
void CreateParticleGeometryStep(particles &vd, std::vector<std::pair<probe_particles, int>> &vp_vec, Parameters params, AuxiliarParameters &auxParams)
{

    Vcluster<> &v_cl = create_vcluster();
    Obstacle *obstacle_ptr = new EmptyObstacle(params);

    // Size of the virtual cartesian grid that defines where to place the particles
    size_t sz[DIM];

    real_number length_small[DIM];
    real_number length_big[DIM];
    // In the case of the new bc we need particles at the wall, for this we need sz_aux
    // We want to put one virtual grid point between each pair of the old ones,
    // so that the new spacing is dp/2, and we can put a fluid particle exactly at the wall
    size_t sz_aux[DIM];

    // Boundary conditions
    size_t bc[DIM];

    // Number of boundary particles in each direction
    size_t Nboundary_big[DIM];

    // Number of fluid particles in each direction
    size_t Nfluid_big[DIM];
    size_t Nfluid_small[DIM];

    // We define the boxes in terms of offstes with respect to the fluid box that goes from 0 to length
    real_number offset_domain[DIM] = {0.0};
    real_number offset_recipient[DIM] = {0.0};
    real_number offset_periodic_fluid[DIM] = {0.0};
    real_number offset_periodic_recipient[DIM] = {0.0};

    // non periodic situation grid of 5 fluid particles and 3 boundary particles
    // We need a virtual grid of 5 + 2*(3+1) particles,
    // therefore the domain is discretized with 13 grid points,
    // when we use DrawParticles::DrawBox we will draw only the particles at the grid positons strictly inside the box,
    // the () repesent the recipient box, and the || represent the fluid box, we can see how this distribution places exactly 5 fluid particles inside and 3 boundary particles
    //           D-(-o--o--o-|-x--x--x--x--x--|-o-o-o-)-D
    // D: domain, o: boundary, x: fluid, --: dp distance
    // in a periodic situation we have the following
    // .....--x--x--D-|-x--x--x--x--x--|-D--x--x--......
    // therefore we need a grid of 5 + 2 particles, and the domain is discretized with 7 grid points

    real_number StepHeight = 4.9;
    params.LengthScale = StepHeight;

    Nfluid_big[0] = 275;
    Nfluid_big[1] = 31;

    Nfluid_small[0] = 100;
    Nfluid_small[1] = 16;

    size_t Nbound = (params.BC_TYPE == NEW_NO_SLIP) ? 1 : 3;
    Nboundary_big[0] = 0;
    Nboundary_big[1] = Nbound;
    // real_number Nboundary_small_up = Nbound;
    // real_number Nboundary_small_down = Nbound + Nfluid_big[1] - Nfluid_small[1];

    bc[0] = PERIODIC;
    bc[1] = NON_PERIODIC;
    params.dp = params.LengthScale / ((real_number)Nfluid_big[1] - Nfluid_small[1]);
    real_number dp = params.dp;
    params.umax = 1.4 * 1e-1;

    params.H = params.Hconst * dp;
    // r_cut = (KERNEL == CUBIC ? 2.0 * H : 3.0 * H);
    params.r_cut = 3.0 * params.H;
    params.Kquintic = (DIM == 3) ? 1.0 / 120.0 / M_PI / params.H / params.H / params.H : 7.0 / 478.0 / M_PI / params.H / params.H;
    params.MassFluid = params.rho0 * (DIM == 3 ? dp * dp * dp : dp * dp);
    params.MassBound = params.rho0 * (DIM == 3 ? dp * dp * dp : dp * dp);
    params.cbar = params.coeff_sound * params.umax;
    params.B = params.rho0 * params.cbar * params.cbar / params.gamma;
    params.Pbackground = params.Bfactor * params.B;
    params.eta = params.nu * params.rho0;
    params.Re = params.umax * 2.0 * 5.2 / params.nu;

    params.gravity = getVectorNorm(params.gravity_vector);

    for (int dim = 0; dim < DIM; dim++)
    {
        if (bc[dim] == NON_PERIODIC) // non periodic, fluid covered by boundary
        {
            params.length[dim] = dp * (Nfluid_big[dim]);
            length_small[dim] = dp * (Nfluid_small[dim]);
            length_big[dim] = dp * (Nfluid_big[dim]);
            sz[dim] = Nfluid_big[dim] + 2 * (Nboundary_big[dim] + 1);
            offset_domain[dim] = (0.5 + Nboundary_big[dim]) * dp;

            if (Nboundary_big[dim] != 0)
                sz_aux[dim] = 2 * Nfluid_big[dim] - 1 + 2 * (2 * Nboundary_big[dim] + 1 + 1);
            else // for a direction with no boundary particles we dont need to add anything
                sz_aux[dim] = sz[dim];

            if (params.BC_TYPE == NEW_NO_SLIP) // Nboundary should only be 0 or 1 if we are using the new bc
                offset_recipient[dim] = 0.25 * Nboundary_big[dim] * dp;
            else if (params.BC_TYPE == NO_SLIP)
                offset_recipient[dim] = Nboundary_big[dim] * dp;
        }
        else // periodic, open ended
        {
            Nfluid_big[dim] -= 1;
            params.length[dim] = dp * (Nfluid_big[dim] + Nfluid_small[dim]);
            length_small[dim] = dp * (Nfluid_small[dim]);
            length_big[dim] = dp * (Nfluid_big[dim]);

            sz[dim] = (Nfluid_big[dim] + Nfluid_small[dim]) + 2;
            sz_aux[dim] = sz[dim];

            offset_domain[dim] = 0.5 * dp;
            offset_periodic_fluid[dim] = 0.75 * dp;
            offset_periodic_recipient[dim] = 0.85 * dp;
        }
    }

    // Define the boxes
    Box<DIM, real_number> domain({-offset_domain[0],
                                  -offset_domain[1]},
                                 {params.length[0] + offset_domain[0],
                                  params.length[1] + offset_domain[1]});

    Box<DIM, real_number> fluid_box_small({0.0,
                                           (Nfluid_big[1] - Nfluid_small[1]) * dp},
                                          {length_small[0],
                                           (Nfluid_big[1] - Nfluid_small[1]) * dp + length_small[1]});

    Box<DIM, real_number> fluid_box_big({length_small[0],
                                         0.0},
                                        {length_small[0] + length_big[0] + offset_periodic_fluid[0],
                                         length_big[1]});

    Box<DIM, real_number> recipient({-offset_recipient[0],
                                     -offset_recipient[1]},
                                    {params.length[0] + offset_recipient[0] + offset_periodic_recipient[0],
                                     params.length[1] + offset_recipient[1] + offset_periodic_recipient[1]});

    // Will only be used in the new bc
    Box<DIM, real_number> recipient_hole_small({offset_recipient[0],
                                                (Nfluid_big[1] - Nfluid_small[1]) * dp + offset_recipient[1]},
                                               {length_small[0] - offset_recipient[0],
                                                (Nfluid_big[1] - Nfluid_small[1]) * dp + length_small[1] - offset_recipient[1]});

    Box<DIM, real_number> recipient_hole_big({length_small[0] + offset_recipient[0],
                                              offset_recipient[1]},
                                             {length_small[0] + length_big[0] - offset_recipient[0] + offset_periodic_fluid[0],
                                              length_big[1] - offset_recipient[1]});

    Box<DIM, real_number> CornerHole{{3 * dp, -3 * dp}, {(3 + Nfluid_small[0] - 6) * dp, (Nfluid_big[1] - Nfluid_small[1] - 3) * dp}};
    Box<DIM, real_number> CornerHole_New{{dp, -1 * dp}, {(1 + Nfluid_small[0] - 2) * dp, (Nfluid_big[1] - Nfluid_small[1]) * dp - 0.5 * dp}};

    // extended boundary around the domain, and the processor domain
    Ghost<DIM, real_number> g(params.r_cut);

    // create particle object
    particles vd_loc(0, domain, bc, g, DEC_GRAN(128));
    vd = vd_loc;

    // correct the number of particles in case of periodicity, we substracted 1 before to accomodate the periodic boundary
    for (int dim = 0; dim < DIM; dim++)
    {
        if (bc[dim] == PERIODIC)
        {
            Nfluid_big[dim] += 1;
            params.length[dim] += dp;
        }
    }

    // Write constants on file
    // WriteParameters(v_cl, params);

    // we want to place probes  in a vertical line at this locations
    Point<DIM, real_number> P1 = {0.75 * Nfluid_small[0] * dp, 0.0};
    Point<DIM, real_number> P2 = {1.2 * Nfluid_small[0] * dp, 0.0};
    Point<DIM, real_number> P3 = {1.4 * Nfluid_small[0] * dp, 0.0};
    Point<DIM, real_number> P4 = {Nfluid_small[0] * dp + Nfluid_big[0] * dp * 0.7, 0.0};

    std::vector<Point<DIM, real_number>> ProbePoints = {P1, P2, P3, P4};

    Point<DIM, real_number> VerticalOffset = {0.0, dp};
    int k0 = 0;
    int kendHeight = Nfluid_big[1];

    // place probes
    if (params.PROBES_ENABLED)
    {
        for (int k = 0; k < 4; k++)
        {
            // create probe object
            Ghost<DIM, real_number> gp(0);
            size_t bc_p[DIM] = {NON_PERIODIC, NON_PERIODIC};
            probe_particles vp_loc(0, domain, bc_p, gp, DEC_GRAN(128));
            openfpm::vector<std::string> names_p({"vx"});
            vp_loc.setPropNames(names_p);

            PlaceProbes(vp_loc, k0, kendHeight, ProbePoints[k], VerticalOffset, obstacle_ptr);
            std::pair<probe_particles, int> tmp = std::make_pair(vp_loc, 0);
            vp_vec.push_back(tmp);
            auxParams.probe_filenames.push_back("probes_" + std::to_string(k) + "_" + auxParams.filename);
        }
    }

    // return an iterator to the fluid particles to add to vd
    auto fluid_it1 = DrawParticles::DrawBox(vd, sz, domain, fluid_box_big);
    auto fluid_it2 = DrawParticles::DrawBox(vd, sz, domain, fluid_box_small);

    // for each particle inside the fluid box ...
    while (fluid_it1.isNext())
    {

        Point<DIM, real_number> iterator_position = fluid_it1.get();

        // ... add a particle ...
        vd.add();
        vd.template getLastProp<vd0_type>() = FLUID;

        // Set properties
        vd.template getLastProp<vd1_rho>() = params.rho0;
        vd.template getLastProp<vd2_pressure>() = 0.0;
        vd.template getLastProp<vd3_drho>() = 0.0;

        for (int xyz = 0; xyz < DIM; xyz++)
        {
            vd.getLastPos()[xyz] = iterator_position.get(xyz);
            vd.template getLastProp<vd7_force_t>()[xyz] = 0.0;
            vd.template getLastProp<vd5_velocity_t>()[xyz] = 0.0;
            vd.template getLastProp<vd8_normal>()[xyz] = 0.0;
            vd.template getLastProp<vd4_velocity>()[xyz] = 0.0;
        }

        vd.template getLastProp<vd9_volume>()[0] = dp;

        // next fluid particle
        ++fluid_it1;
    }
    // for each particle inside the fluid box ...
    while (fluid_it2.isNext())
    {

        Point<DIM, real_number> iterator_position = fluid_it2.get();

        // ... add a particle ...
        vd.add();
        vd.template getLastProp<vd0_type>() = FLUID;

        // Set properties
        vd.template getLastProp<vd1_rho>() = params.rho0;
        vd.template getLastProp<vd2_pressure>() = 0.0;
        vd.template getLastProp<vd3_drho>() = 0.0;

        for (int xyz = 0; xyz < DIM; xyz++)
        {
            vd.getLastPos()[xyz] = iterator_position.get(xyz);
            vd.template getLastProp<vd7_force_t>()[xyz] = 0.0;
            vd.template getLastProp<vd5_velocity_t>()[xyz] = 0.0;
            vd.template getLastProp<vd8_normal>()[xyz] = 0.0;
            vd.template getLastProp<vd4_velocity>()[xyz] = 0.0;
        }

        vd.template getLastProp<vd9_volume>()[0] = dp;

        // next fluid particle
        ++fluid_it2;
    }

    // Now place solid walls
    openfpm::vector<Box<DIM, real_number>> holes;

    if (params.BC_TYPE == NEW_NO_SLIP)
    {
        holes.add(recipient_hole_small);
        holes.add(recipient_hole_big);
        holes.add(CornerHole_New);
        sz[0] = sz_aux[0];
        sz[1] = sz_aux[1];
    }
    else if (params.BC_TYPE == NO_SLIP)
    {
        holes.add(fluid_box_big);
        holes.add(fluid_box_small);
        holes.add(CornerHole);
    }
    Box<DIM, real_number> hole_get0 = holes.get(0);
    Box<DIM, real_number> hole_get1 = holes.get(1);
    auto bound_box = DrawParticles::DrawSkin(vd, sz, domain, holes, recipient);

    if (bc[0] != PERIODIC || bc[1] != PERIODIC) // no walls in all periodic scenario
    {
        while (bound_box.isNext())
        {
            Point<DIM, real_number> position = bound_box.get();

            // periodic bc, with no boundary particles in y direction has a bug, it puts 3 extra particles outside in the y direction
            // When running on multiple cores, with this we check if particle is outside the recipient box
            // Another bug places boundary particles in the correct plane, but inside the fluid box;
            // the first bug seems to be fixed
            // if (!recipient.isInside((position)))
            // {
            // 	++bound_box;
            // 	continue;
            // }
            // if (bc[0] == PERIODIC && position.get(0) > dp / 2.0 && position.get(0) < length[0] - dp / 2.0)
            // {
            // 	++bound_box;
            // 	continue;
            // }

            if (hole_get0.isInside(position) || hole_get1.isInside(position))
            {
                ++bound_box;
                continue;
            }

            if (params.BC_TYPE == NEW_NO_SLIP)
            {
                // Check if x and y coordinates are multiples of dp, keep multiples, discard the rest
                // real_number remx = fmod(position.get(0), dp);
                real_number remy = fmod(position.get(1), dp);
                real_number tol = 0.5 * dp * 10e-2;

                // if (remx > tol && remx < dp - tol)
                // {
                // 	++bound_box;
                // 	continue;
                // }
                if (remy > tol && remy < dp - tol)
                {
                    ++bound_box;
                    continue;
                }
            }
            vd.add();

            vd.template getLastProp<vd0_type>() = BOUNDARY;
            vd.template getLastProp<vd1_rho>() = params.rho0;
            vd.template getLastProp<vd2_pressure>() = 0.0;
            vd.template getLastProp<vd3_drho>() = 0.0;

            for (int xyz = 0; xyz < DIM; xyz++)
            {
                vd.getLastPos()[xyz] = bound_box.get().get(xyz);
                vd.template getLastProp<vd6_force>()[xyz] = 0.0;
                vd.template getLastProp<vd7_force_t>()[xyz] = 0.0;
                vd.template getLastProp<vd5_velocity_t>()[xyz] = 0.0;
                vd.template getLastProp<vd8_normal>()[xyz] = 0.0;
                if (position.get(1) < dp / 4.0) // bottom wall
                {
                    vd.template getLastProp<vd4_velocity>()[xyz] = params.vw_bottom[xyz];
                }
                else if (position.get(1) > params.length[1] - dp / 4.0) // top wall
                {
                    vd.template getLastProp<vd4_velocity>()[xyz] = params.vw_top[xyz];
                }
            }

            vd.template getLastProp<vd9_volume>()[0] = dp;

            ++bound_box;
        }
    }
}

void CreateParticleGeometryDamBreak(particles &vd, std::vector<std::pair<probe_particles, int>> &vp_vec, Parameters &params, AuxiliarParameters &auxParams)
{

    Vcluster<> &v_cl = create_vcluster();

    // Non periodic situation grid of 5 fluid particles and 3 boundary particles
    // We need a virtual grid of 5 + 2*(3+1) particles,
    // therefore the domain is discretized with 13 grid points,
    // when we use DrawParticles::DrawBox we will draw only the particles at the grid positons strictly inside the box,
    // the () repesent the recipient box, and the || represent the fluid box, we can see how this distribution places exactly 5 fluid particles inside and 3 boundary particles
    //           D-(-o--o--o-|-x--x--x--x--x--|-o-o-o-)-D
    // D: domain, o: boundary, x: fluid, --: dp distance
    // in a periodic situation we have the following
    // .....--x--x--D-|-x--x--x--x--x--|-D--x--x--......
    // therefore we need a grid of 5 + 2 particles, and the domain is discretized with 7 grid points

    // Size of the virtual cartesian grid that defines where to place the particles
    size_t sz[DIM];

    real_number dp = params.dp;
    real_number refine_factor = params.rf;

    // Now define the iterator boxes
    // We define the boxes in terms of offstes with respect to the fluid box that goes from 0 to length
    real_number offset_domain_left[DIM] = {0.0};
    real_number offset_domain_right[DIM] = {0.0};
    real_number offset_recipient[DIM] = {0.0};
    real_number offset_periodic_fluid[DIM] = {0.0};
    real_number offset_periodic_recipient[DIM] = {0.0};

    for (int xyz = 0; xyz < DIM; xyz++)
    {
        if (params.bc[xyz] == NON_PERIODIC) // non periodic, fluid covered by boundary
        {
            sz[xyz] = params.Nfluid[xyz] + 2 * (params.Nboundary[xyz] + 1);
            offset_domain_left[xyz] = (0.5 + params.Nboundary[xyz]) * dp;
            offset_domain_right[xyz] = (0.5 + params.Nboundary[xyz]) * dp;

            if (params.BC_TYPE == NEW_NO_SLIP) // Nboundary should only be 0 or 1 if we are using the new bc
                offset_recipient[xyz] = 0.1 * params.Nboundary[xyz] * dp;
            else if (params.BC_TYPE == NO_SLIP)
                offset_recipient[xyz] = params.Nboundary[xyz] * dp;
        }
        else // periodic, open ended
        {
            sz[xyz] = params.Nfluid[xyz] + 1;

            offset_domain_left[xyz] = 0.0;
            offset_domain_right[xyz] = dp;
            offset_periodic_fluid[xyz] = 0.75 * dp;
            offset_periodic_recipient[xyz] = 0.85 * dp;
        }
    }

    // Define the boxes
    Box<DIM, real_number> domain({-offset_domain_left[0],
                                  -offset_domain_left[1]},
                                 {params.length[0] + offset_domain_right[0],
                                  params.length[1] + offset_domain_right[1]});

    Box<DIM, real_number> fluid_hole({0.0,
                                      0.0},
                                     {params.length[0] + offset_periodic_fluid[0],
                                      params.length[1] + offset_periodic_fluid[1]});

    real_number wlx = params.waterB;
    real_number wly = params.waterH;

    Box<DIM, real_number> fluid_box({0.0,
                                     0.0},
                                    {wlx,
                                     wly});

    Box<DIM, real_number> recipient({-offset_recipient[0],
                                     -offset_recipient[1]},
                                    {params.length[0] + offset_recipient[0] + offset_periodic_recipient[0],
                                     params.length[1] + offset_recipient[1] + offset_periodic_recipient[1]});

    for (int xyz = 0; xyz < DIM; xyz++) // correct length in periodic case
    {
        if (params.bc[xyz] == PERIODIC)
            params.length[xyz] += dp;
    }

    // extended boundary around the domain, and the processor domain
    Ghost<DIM, real_number> g(params.r_cut);

    // create particle object
    particles vd_loc(0, domain, params.bc, g, DEC_GRAN(128));
    // vd is argument passed as reference we want to fill with particles
    vd = vd_loc;

    // Write constants on file
    // WriteParameters(v_cl, params);

    // place probes
    // if (params.PROBES_ENABLED)
    // {
    //     // we want to place probes  in a vertical line at this locations
    //     Point<DIM, real_number> EndChannel = {params.length[0], 0.0};
    //     Point<DIM, real_number> HalfChannel = {params.length[0] / 2.0, 0.0};
    //     Point<DIM, real_number> HalfHeight = {0.0, params.length[1] / 2.0};
    //     Point<DIM, real_number> VerticalOffset = {0.0, dp};
    //     Point<DIM, real_number> HorizontalOffset = {dp, 0.0};
    //     int k0 = 0;
    //     int kendHeight = params.Nfluid[1] + 1;
    //     int kendWidth = params.Nfluid[0] + 1;

    //     std::vector<Point<DIM, real_number>> ProbePoints; // start points for the PlaceProbes function
    //     std::vector<int> ProbeComponents;                 // component to measure 0 for x 1 for y
    //     std::vector<Point<DIM, real_number>> Offsets;
    //     std::vector<int> maxIters;

    //     if (params.SCENARIO == CAVITY)
    //     {
    //         ProbePoints.push_back(HalfChannel);
    //         ProbePoints.push_back(HalfHeight);

    //         ProbeComponents.push_back(0); // measure x velocity
    //         ProbeComponents.push_back(1); // measure y velocity

    //         Offsets.push_back(VerticalOffset);
    //         Offsets.push_back(HorizontalOffset);

    //         maxIters.push_back(kendHeight);
    //         maxIters.push_back(kendWidth);
    //     }
    //     else
    //     {
    //         ProbePoints.push_back(HalfChannel);
    //         ProbePoints.push_back(EndChannel);

    //         ProbeComponents.push_back(0); // measure x velocity
    //         ProbeComponents.push_back(0); // measure x velocity

    //         Offsets.push_back(VerticalOffset);
    //         Offsets.push_back(VerticalOffset);

    //         maxIters.push_back(kendHeight);
    //         maxIters.push_back(kendHeight);
    //     }

    //     for (unsigned int k = 0; k < ProbePoints.size(); k++)
    //     {
    //         // create probe object
    //         Ghost<DIM, real_number> gp(0);
    //         size_t bc_p[DIM] = {NON_PERIODIC, NON_PERIODIC};
    //         probe_particles vp_loc(0, domain, bc_p, gp, DEC_GRAN(512));
    //         if (ProbeComponents[k] == 0)
    //         {
    //             openfpm::vector<std::string> names_p = {"vx"};
    //             vp_loc.setPropNames(names_p);
    //         }
    //         else if (ProbeComponents[k] == 1)
    //         {
    //             openfpm::vector<std::string> names_p = {"vy"};
    //             vp_loc.setPropNames(names_p);
    //         }

    //         if (v_cl.getProcessUnitID() == 0)
    //         {
    //             PlaceProbes(vp_loc, k0, maxIters[k], ProbePoints[k], Offsets[k]);
    //         }
    //         std::pair<probe_particles, int> tmp = std::make_pair(vp_loc, ProbeComponents[k]);
    //         vp_vec.push_back(tmp);
    //         params.probe_filenames.push_back("probes_" + std::to_string(k) + "_" + params.filename);
    //     }
    // }

    // Add the obstacle/walls as marker particles only on processor 0
    if (params.BC_TYPE == NEW_NO_SLIP && v_cl.getProcessUnitID() == 0)
    {

        real_number dx_wall_x = dp / refine_factor;
        int Nwall_x = ceil(params.length[0] / dx_wall_x);
        dx_wall_x = params.length[0] / Nwall_x;
        Point<DIM, real_number> X_Offset = {dx_wall_x, 0.0};

        real_number dx_wall_y = dp / refine_factor;
        int Nwall_y = ceil(params.length[1] / dx_wall_y);
        dx_wall_y = params.length[1] / Nwall_y;
        Point<DIM, real_number> Y_Offset = {0.0, dx_wall_y};

        Point<DIM, real_number> LL_corner = {0.0, 0.0};
        Point<DIM, real_number> LR_corner = {params.length[0], 0.0};

        Point<DIM, real_number> UL_corner = {0.0, params.length[1]};

        // Top And Bottom Walls
        AddFlatWallModNewBC(vd, 0, 1, LL_corner, X_Offset, dx_wall_x, {0.0, 0.0}, params.vw_bottom, params, BOUNDARY, {10.0 * 1.0, 10.0 * 1.0}, 0.0);
        AddFlatWallModNewBC(vd, Nwall_x, Nwall_x + 1, LL_corner, X_Offset, dx_wall_x, {0.0, 0.0}, params.vw_bottom, params, BOUNDARY, {10.0 * -1.0, 10.0 * 1.0}, 0.0);
        AddFlatWallModNewBC(vd, 1, Nwall_x, LL_corner, X_Offset, dx_wall_x, {0.0, 0.0}, params.vw_bottom, params, BOUNDARY, {0.0, 10.0 * 1.0}, 0.0);

        AddFlatWallModNewBC(vd, 0, 1, UL_corner, X_Offset, dx_wall_x, {0.0, 0.0}, params.vw_top, params, BOUNDARY, {10.0 * 1.0, 10.0 * -1.0}, 0.0);
        AddFlatWallModNewBC(vd, Nwall_x, Nwall_x + 1, UL_corner, X_Offset, dx_wall_x, {0.0, 0.0}, params.vw_top, params, BOUNDARY, {-10.0 * 1.0, 10.0 * -1.0}, 0.0);
        AddFlatWallModNewBC(vd, 1, Nwall_x, UL_corner, X_Offset, dx_wall_x, {0.0, 0.0}, params.vw_top, params, BOUNDARY, {0.0, 10.0 * -1.0}, 0.0);

        // Left And Right Walls
        AddFlatWallModNewBC(vd, 1, Nwall_y, LL_corner, Y_Offset, dx_wall_y, {0.0, 0.0}, {0.0, 0.0}, params, BOUNDARY, {10.0 * 1.0, 0.0}, 0.0);
        AddFlatWallModNewBC(vd, 1, Nwall_y, LR_corner, Y_Offset, dx_wall_y, {0.0, 0.0}, {0.0, 0.0}, params, BOUNDARY, {10.0 * -1.0, 0.0}, 0.0);
    }

    // return an iterator to the fluid particles to add to vd
    auto fluid_it = DrawParticles::DrawBox(vd, sz, domain, fluid_box);

    // for each particle inside the fluid box ...
    while (fluid_it.isNext())
    {

        Point<DIM, real_number> iterator_position = fluid_it.get();

        // ... add a particle ...
        vd.add();
        vd.template getLastProp<vd0_type>() = FLUID;
        vd.template getLastProp<vd10_omega>() = 0.0;
        for (int xyz = 0; xyz < DIM; xyz++)
        {
            vd.template getLastProp<vd4_velocity>()[xyz] = 0.0;
            vd.template getLastProp<vd7_force_t>()[xyz] = 0.0;
        }

        // Set properties
        vd.template getLastProp<vd1_rho>() = params.rho0;
        vd.template getLastProp<vd2_pressure>() = 0.0;
        vd.template getLastProp<vd3_drho>() = 0.0;

        for (int xyz = 0; xyz < DIM; xyz++)
        {
            vd.getLastPos()[xyz] = iterator_position.get(xyz);
            vd.template getLastProp<vd5_velocity_t>()[xyz] = 0.0;
            vd.template getLastProp<vd8_normal>()[xyz] = 0.0;
        }

        vd.template getLastProp<vd9_volume>()[0] = dp;
        vd.template getLastProp<vd9_volume>()[1] = 0.0;
        vd.template getLastProp<vd9_volume>()[2] = 0.0;

        // next fluid particle
        ++fluid_it;
    }

    // Now place solid walls using iterators (only for OLD BC)

    if (params.BC_TYPE == NO_SLIP)
    {

        openfpm::vector<Box<DIM, real_number>> holes;
        holes.add(fluid_hole);
        Box<DIM, real_number> hole_box = holes.get(0);
        auto bound_box = DrawParticles::DrawSkin(vd, sz, domain, holes, recipient);

        if (params.bc[0] != PERIODIC || params.bc[1] != PERIODIC) // no walls in all periodic scenario
        {
            while (bound_box.isNext())
            {
                Point<DIM, real_number> position = bound_box.get();

                // periodic bc, with no boundary particles in y direction has a bug, it puts 3 extra particles outside in the y direction
                // When running on multiple cores, with this we check if particle is outside the recipient box
                // Another bug places boundary particles in the correct plane, but inside the fluid box;
                // if (bc[0] == PERIODIC && position.get(0) > dp / 2.0 && position.get(0) < length[0] - dp / 2.0)
                // {
                // 	++bound_box;
                // 	continue;
                // }

                if (!recipient.isInside((position)))
                {
                    ++bound_box;
                    continue;
                }
                if (hole_box.isInside(position))
                {
                    ++bound_box;
                    continue;
                }

                // if (params.BC_TYPE == NEW_NO_SLIP && (params.bc[0] == NON_PERIODIC && params.bc[1] == NON_PERIODIC))
                // {
                // 	// Check if x and z coordinates are multiples of dp, keep multiples, discard the rest
                // 	real_number remx = fmod(position.get(0), dp);
                // 	real_number remz = fmod(position.get(1), dp);
                // 	real_number tol = 0.5 * dp * 10e-2;

                // 	if (remx > tol && remx < dp - tol)
                // 	{
                // 		++bound_box;
                // 		continue;
                // 	}
                // 	if (remz > tol && remz < dp - tol)
                // 	{
                // 		++bound_box;
                // 		continue;
                // 	}
                // }
                vd.add();

                vd.template getLastProp<vd0_type>() = BOUNDARY;
                vd.template getLastProp<vd1_rho>() = params.rho0;
                vd.template getLastProp<vd2_pressure>() = 0.0;
                vd.template getLastProp<vd3_drho>() = 0.0;

                for (int xyz = 0; xyz < DIM; xyz++)
                {
                    vd.getLastPos()[xyz] = bound_box.get().get(xyz);
                    vd.template getLastProp<vd6_force>()[xyz] = 0.0;
                    vd.template getLastProp<vd7_force_t>()[xyz] = 0.0;
                    vd.template getLastProp<vd5_velocity_t>()[xyz] = 0.0;
                    vd.template getLastProp<vd8_normal>()[xyz] = 0.0;
                    if (position.get(1) < dp / 4.0) // bottom wall
                    {
                        vd.template getLastProp<vd4_velocity>()[xyz] = params.vw_bottom[xyz];
                    }
                    else if (position.get(1) > params.length[1] - dp / 4.0) // top wall
                    {
                        vd.template getLastProp<vd4_velocity>()[xyz] = params.vw_top[xyz];
                    }
                }

                vd.template getLastProp<vd9_volume>()[0] = dp;

                ++bound_box;
            }
        }
    }
}

void CreateParticleGeometryDamBreakAdj(particles &vd, std::vector<std::pair<probe_particles, int>> &vp_vec, Parameters &params, AuxiliarParameters &auxParams)
{

    Vcluster<> &v_cl = create_vcluster();

    // Non periodic situation grid of 5 fluid particles and 3 boundary particles
    // We need a virtual grid of 5 + 2*(3+1) particles,
    // therefore the domain is discretized with 13 grid points,
    // when we use DrawParticles::DrawBox we will draw only the particles at the grid positons strictly inside the box,
    // the () repesent the recipient box, and the || represent the fluid box, we can see how this distribution places exactly 5 fluid particles inside and 3 boundary particles
    //           D-(-o--o--o-|-x--x--x--x--x--|-o-o-o-)-D
    // D: domain, o: boundary, x: fluid, --: dp distance
    // in a periodic situation we have the following
    // .....--x--x--D-|-x--x--x--x--x--|-D--x--x--......
    // therefore we need a grid of 5 + 2 particles, and the domain is discretized with 7 grid points

    // Size of the virtual cartesian grid that defines where to place the particles
    size_t sz[DIM];

    real_number dp = params.dp;
    real_number refine_factor = params.rf;

    // Now define the iterator boxes
    // We define the boxes in terms of offstes with respect to the fluid box that goes from 0 to length
    real_number offset_domain_left[DIM] = {0.0};
    real_number offset_domain_right[DIM] = {0.0};
    real_number offset_recipient[DIM] = {0.0};
    real_number offset_periodic_fluid[DIM] = {0.0};
    real_number offset_periodic_recipient[DIM] = {0.0};

    for (int xyz = 0; xyz < DIM; xyz++)
    {
        if (params.bc[xyz] == NON_PERIODIC) // non periodic, fluid covered by boundary
        {
            sz[xyz] = params.Nfluid[xyz] + 2 * (params.Nboundary[xyz] + 1);
            offset_domain_left[xyz] = (0.5 + params.Nboundary[xyz]) * dp;
            offset_domain_right[xyz] = (0.5 + params.Nboundary[xyz]) * dp;

            if (params.BC_TYPE == NEW_NO_SLIP) // Nboundary should only be 0 or 1 if we are using the new bc
                offset_recipient[xyz] = 0.1 * params.Nboundary[xyz] * dp;
            else if (params.BC_TYPE == NO_SLIP)
                offset_recipient[xyz] = params.Nboundary[xyz] * dp;
        }
        else // periodic, open ended
        {
            sz[xyz] = params.Nfluid[xyz] + 1;

            offset_domain_left[xyz] = 0.0;
            offset_domain_right[xyz] = dp;
            offset_periodic_fluid[xyz] = 0.75 * dp;
            offset_periodic_recipient[xyz] = 0.85 * dp;
        }
    }

    // Define the boxes
    real_number auxiliary_tank_width = 22.0 * dp;
    real_number auxiliary_tank_height = 6.0 * dp;

    sz[0] += 44;

    Box<DIM, real_number> domain({-offset_domain_left[0],
                                  -offset_domain_left[1]},
                                 {params.length[0] + 2.0 * auxiliary_tank_width + offset_domain_right[0],
                                  params.length[1] + offset_domain_right[1]});

    Box<DIM, real_number> fluid_hole({0.0,
                                      0.0},
                                     {params.length[0] + offset_periodic_fluid[0],
                                      params.length[1] + offset_periodic_fluid[1]});

    real_number wlx = params.waterB;
    real_number wly = params.waterH;

    Box<DIM, real_number> fluid_box({0.0,
                                     0.0},
                                    {wlx,
                                     wly});

    Box<DIM, real_number> second_fluid_box = {{params.length[0] + dp + auxiliary_tank_width,
                                               auxiliary_tank_height},
                                              {params.length[0] + 2.0 * auxiliary_tank_width + offset_periodic_fluid[0],
                                               params.length[1] + offset_periodic_fluid[1] - 2.0 * auxiliary_tank_height}};

    Box<DIM, real_number> recipient({-offset_recipient[0],
                                     -offset_recipient[1]},
                                    {params.length[0] + offset_recipient[0] + offset_periodic_recipient[0],
                                     params.length[1] + offset_recipient[1] + offset_periodic_recipient[1]});

    for (int xyz = 0; xyz < DIM; xyz++) // correct length in periodic case
    {
        if (params.bc[xyz] == PERIODIC)
            params.length[xyz] += dp;
    }

    // extended boundary around the domain, and the processor domain
    Ghost<DIM, real_number> g(params.r_cut);

    // create particle object
    particles vd_loc(0, domain, params.bc, g, DEC_GRAN(128));
    // vd is argument passed as reference we want to fill with particles
    vd = vd_loc;

    // Add the obstacle/walls as marker particles only on processor 0
    if (params.BC_TYPE == NEW_NO_SLIP && v_cl.getProcessUnitID() == 0)
    {

        real_number dx_wall_x = dp / refine_factor;
        int Nwall_x = ceil(params.length[0] / dx_wall_x);
        dx_wall_x = params.length[0] / Nwall_x;
        Point<DIM, real_number> X_Offset = {dx_wall_x, 0.0};

        real_number dx_wall_y = dp / refine_factor;
        int Nwall_y = ceil(params.length[1] / dx_wall_y);
        dx_wall_y = params.length[1] / Nwall_y;
        Point<DIM, real_number> Y_Offset = {0.0, dx_wall_y};

        Point<DIM, real_number> LL_corner = {0.0, 0.0};
        Point<DIM, real_number> LR_corner = {params.length[0], 0.0};
        Point<DIM, real_number> UL_corner = {0.0, params.length[1]};

        Point<DIM, real_number> ColumnLL = {params.length[0] + dp + auxiliary_tank_width, auxiliary_tank_height};
        Point<DIM, real_number> ColumnLR = {params.length[0] + 2.0 * auxiliary_tank_width, auxiliary_tank_height};
        Point<DIM, real_number> ColumnUL = {params.length[0] + dp + auxiliary_tank_width, params.length[1] - auxiliary_tank_height};

        int Nwall_x_small = ceil(auxiliary_tank_width / dx_wall_x);
        int Nwall_y_small = ceil((params.length[1] - 2 * auxiliary_tank_height) / dx_wall_y);
        real_number fac = 0.5;
        // Top And Bottom Walls
        AddFlatWallModNewBC(vd, 0, 1, LL_corner, X_Offset, dx_wall_x, {0.0, 0.0}, params.vw_bottom, params, BOUNDARY, {fac * 1.0, fac * 1.0}, 0.0, 0.0);
        AddFlatWallModNewBC(vd, Nwall_x, Nwall_x + 1, LL_corner, X_Offset, dx_wall_x, {0.0, 0.0}, params.vw_bottom, params, BOUNDARY, {fac * -1.0, fac * 1.0}, 0.0, 0.0);
        AddFlatWallModNewBC(vd, 1, Nwall_x, LL_corner, X_Offset, dx_wall_x, {0.0, 0.0}, params.vw_bottom, params, BOUNDARY, {0.0, fac * 1.0}, 0.0, 0.0);

        AddFlatWallModNewBC(vd, 0, 1, UL_corner, X_Offset, dx_wall_x, {0.0, 0.0}, params.vw_top, params, BOUNDARY, {fac * 1.0, fac * -1.0}, 0.0), 0.0;
        AddFlatWallModNewBC(vd, Nwall_x, Nwall_x + 1, UL_corner, X_Offset, dx_wall_x, {0.0, 0.0}, params.vw_top, params, BOUNDARY, {-fac * 1.0, fac * -1.0}, 0.0, 0.0);
        AddFlatWallModNewBC(vd, 1, Nwall_x, UL_corner, X_Offset, dx_wall_x, {0.0, 0.0}, params.vw_top, params, BOUNDARY, {0.0, fac * -1.0}, 0.0, 0.0);

        // Left And Right Walls
        AddFlatWallModNewBC(vd, 1, Nwall_y, LL_corner, Y_Offset, dx_wall_y, {0.0, 0.0}, {0.0, 0.0}, params, BOUNDARY, {fac * 1.0, 0.0}, 0.0, 0.0);
        AddFlatWallModNewBC(vd, 1, Nwall_y, LR_corner, Y_Offset, dx_wall_y, {0.0, 0.0}, {0.0, 0.0}, params, BOUNDARY, {fac * -1.0, 0.0}, 0.0, 0.0);

        // Auxiliary Tank
        AddFlatWallModNewBC(vd, 0, 1, ColumnLL, X_Offset, dx_wall_x, {0.0, 0.0}, {0.0, 0.0}, params, BOUNDARY, {fac * 1.0, fac * 1.0}, 0.0, 0.0);
        AddFlatWallModNewBC(vd, 1, Nwall_x_small - 1, ColumnLL, X_Offset, dx_wall_x, {0.0, 0.0}, {0.0, 0.0}, params, BOUNDARY, {0.0, fac * 1.0}, 0.0, 0.0);
        AddFlatWallModNewBC(vd, Nwall_x_small - 1, Nwall_x_small, ColumnLL, X_Offset, dx_wall_x, {0.0, 0.0}, {0.0, 0.0}, params, BOUNDARY, {fac * -1.0, fac * 1.0}, 0.0, 0.0);

        AddFlatWallModNewBC(vd, 0, 1, ColumnUL, X_Offset, dx_wall_x, {0.0, 0.0}, {0.0, 0.0}, params, BOUNDARY, {fac * 1.0, -fac * 1.0}, 0.0, 0.0);
        AddFlatWallModNewBC(vd, 1, Nwall_x_small - 1, ColumnUL, X_Offset, dx_wall_x, {0.0, 0.0}, {0.0, 0.0}, params, BOUNDARY, {0.0, -fac * 1.0}, 0.0, 0.0);
        AddFlatWallModNewBC(vd, Nwall_x_small - 1, Nwall_x_small, ColumnUL, X_Offset, dx_wall_x, {0.0, 0.0}, {0.0, 0.0}, params, BOUNDARY, {fac * -1.0, -fac * 1.0}, 0.0, 0.0);

        AddFlatWallModNewBC(vd, 1, Nwall_y_small, ColumnLL, Y_Offset, dx_wall_y, {0.0, 0.0}, {0.0, 0.0}, params, BOUNDARY, {fac * 1.0, 0.0}, 0.0, 0.0);
        AddFlatWallModNewBC(vd, 1, Nwall_y_small, ColumnLR, Y_Offset, dx_wall_y, {0.0, 0.0}, {0.0, 0.0}, params, BOUNDARY, {fac * -1.0, 0.0}, 0.0, 0.0);
    }

    // return an iterator to the fluid particles to add to vd
    auto fluid_it = DrawParticles::DrawBox(vd, sz, domain, fluid_box);
    auto fluid_it2 = DrawParticles::DrawBox(vd, sz, domain, second_fluid_box);

    // for each particle inside the fluid box ...
    while (fluid_it.isNext())
    {

        Point<DIM, real_number> iterator_position = fluid_it.get();

        // ... add a particle ...
        vd.add();
        vd.template getLastProp<vd0_type>() = FLUID;
        vd.template getLastProp<vd10_omega>() = 0.0;
        for (int xyz = 0; xyz < DIM; xyz++)
        {
            vd.template getLastProp<vd4_velocity>()[xyz] = 0.0;
            vd.template getLastProp<vd7_force_t>()[xyz] = 0.0;
        }

        // Set properties
        vd.template getLastProp<vd1_rho>() = params.rho0;
        vd.template getLastProp<vd2_pressure>() = 0.0;
        vd.template getLastProp<vd3_drho>() = 0.0;

        for (int xyz = 0; xyz < DIM; xyz++)
        {
            vd.getLastPos()[xyz] = iterator_position.get(xyz);
            vd.template getLastProp<vd5_velocity_t>()[xyz] = 0.0;
            vd.template getLastProp<vd8_normal>()[xyz] = 0.0;
        }

        vd.template getLastProp<vd9_volume>()[0] = dp;
        vd.template getLastProp<vd9_volume>()[1] = 0.0;
        vd.template getLastProp<vd9_volume>()[2] = 0.0;

        // next fluid particle
        ++fluid_it;
    }
    // for each particle inside the fluid box ...
    while (fluid_it2.isNext())
    {

        Point<DIM, real_number> iterator_position = fluid_it2.get();

        // ... add a particle ...
        vd.add();
        vd.template getLastProp<vd0_type>() = FLUID;
        vd.template getLastProp<vd10_omega>() = 0.0;
        for (int xyz = 0; xyz < DIM; xyz++)
        {
            vd.template getLastProp<vd4_velocity>()[xyz] = 0.0;
            vd.template getLastProp<vd7_force_t>()[xyz] = 0.0;
        }

        // Set properties
        vd.template getLastProp<vd1_rho>() = params.rho0;
        vd.template getLastProp<vd2_pressure>() = 0.0;
        vd.template getLastProp<vd3_drho>() = 0.0;

        for (int xyz = 0; xyz < DIM; xyz++)
        {
            vd.getLastPos()[xyz] = iterator_position.get(xyz);
            vd.template getLastProp<vd5_velocity_t>()[xyz] = 0.0;
            vd.template getLastProp<vd8_normal>()[xyz] = 0.0;
        }

        vd.template getLastProp<vd9_volume>()[0] = dp;
        vd.template getLastProp<vd9_volume>()[1] = 0.0;
        vd.template getLastProp<vd9_volume>()[2] = 0.0;

        // next fluid particle
        ++fluid_it2;
    }

    // Create a Cell list object for CPU calculations
    auto NN = vd.getCellList(params.r_cut);
    vd.updateCellList(NN);

    CalcFluidVec(vd, NN, params);
    vd.ghost_get<vd0_type, vd1_rho, vd2_pressure, vd4_velocity, vd5_velocity_t, vd6_force, vd7_force_t, vd8_normal, vd9_volume, vd10_omega, vd11_vorticity>();

    CalcNormalVec(vd, NN, params);
    vd.ghost_get<vd0_type, vd1_rho, vd2_pressure, vd4_velocity, vd5_velocity_t, vd6_force, vd7_force_t, vd8_normal, vd9_volume, vd10_omega, vd11_vorticity>();

    CalcCurvature(vd, NN, params);
    vd.ghost_get<vd0_type, vd1_rho, vd2_pressure, vd4_velocity, vd5_velocity_t, vd6_force, vd7_force_t, vd8_normal, vd9_volume, vd10_omega, vd11_vorticity>();

    CalcVolume(vd, params.dp);
    vd.ghost_get<vd0_type, vd1_rho, vd2_pressure, vd4_velocity, vd5_velocity_t, vd6_force, vd7_force_t, vd8_normal, vd9_volume, vd10_omega, vd11_vorticity>();

    vd.write_frame("frametest", 0, 0.0, VTK_WRITER);
    // now move aux tank beside the dam break
    auto it = vd.getDomainIterator();
    while (it.isNext())
    {
        auto akey = it.get();

        Point<DIM, real_number> position = vd.getPos(akey);
        real_number xthreshold = params.length[0] + 0.5 * auxiliary_tank_width;
        if (position.get(0) > xthreshold)
        {
            vd.getPos(akey)[0] -= (auxiliary_tank_width + 0.99 * dp);
        }

        ++it;
    }
}

void CreateParticleGeometryCavity(particles &vd, std::vector<std::pair<probe_particles, int>> &vp_vec, Obstacle *&obstacle_ptr, Parameters &params, AuxiliarParameters &auxParams)
{

    Vcluster<> &v_cl = create_vcluster();
    obstacle_ptr = new EmptyObstacle(params);

    // Size of the virtual cartesian grid that defines where to place the particles
    size_t sz[DIM];
    real_number dp = params.dp;

    // Now define the iterator boxes
    // We define the boxes in terms of offstes with respect to the fluid box that goes from 0 to length
    real_number offset_domain_left[DIM] = {0.0};
    real_number offset_domain_right[DIM] = {0.0};
    real_number offset_recipient[DIM] = {0.0};
    // Nboundary[0] = 3;
    // Nboundary[1] = 3;

    sz[0] = params.Nfluid[0] + 2 * (params.Nboundary[0] + 1) + 6;
    sz[1] = params.Nfluid[1] + 2 * (params.Nboundary[1] + 1);

    offset_domain_left[0] = (3.5 + params.Nboundary[0]) * dp;
    offset_domain_left[1] = (0.5 + params.Nboundary[1]) * dp;

    offset_domain_right[0] = (3.5 + params.Nboundary[0]) * dp;
    offset_domain_right[1] = (0.5 + params.Nboundary[1]) * dp;

    if (params.BC_TYPE == NEW_NO_SLIP) // Nboundary should only be 0 or 1 if we are using the new bc
    {
        offset_recipient[0] = 0.1 * params.Nboundary[0] * dp;
        offset_recipient[1] = 0.1 * params.Nboundary[1] * dp;
    }
    else if (params.BC_TYPE == NO_SLIP)
    {
        offset_recipient[0] = params.Nboundary[0] * dp;
        offset_recipient[1] = params.Nboundary[1] * dp;
    }

    // Define the boxes
    Box<DIM, real_number> domain({-offset_domain_left[0],
                                  -offset_domain_left[1]},
                                 {params.length[0] + offset_domain_right[0],
                                  params.length[1] + offset_domain_right[1]});

    Box<DIM, real_number> fluid_box({0.0,
                                     0.0},
                                    {params.length[0],
                                     params.length[1]});

    Box<DIM, real_number> recipient({-offset_recipient[0],
                                     -offset_recipient[1]},
                                    {params.length[0] + offset_recipient[0],
                                     params.length[1] + offset_recipient[1]});

    // for (int xyz = 0; xyz < DIM; xyz++) // correct length in periodic case
    // {
    //     if (params.bc[xyz] == PERIODIC)
    //         params.length[xyz] += dp;
    // }

    // extended boundary around the domain, and the processor domain
    Ghost<DIM, real_number> g(params.r_cut);

    // create particle object
    particles vd_loc(0, domain, params.bc, g, DEC_GRAN(128));
    // vd is argument passed as reference we want to fill with particles
    vd = vd_loc;

    // place probes
    if (params.PROBES_ENABLED)
    {
        // we want to place probes in a vertical line at this locations
        Point<DIM, real_number> HalfChannel = {params.length[0] / 2.0, 0.0};
        Point<DIM, real_number> HalfHeight = {0.0, params.length[1] / 2.0};
        Point<DIM, real_number> VerticalOffset = {0.0, dp};
        Point<DIM, real_number> HorizontalOffset = {dp, 0.0};
        int k0 = 0;
        int kendHeight = params.Nfluid[1] + 1;
        int kendWidth = params.Nfluid[0] + 1;

        std::vector<Point<DIM, real_number>> ProbePoints; // start points for the PlaceProbes function
        std::vector<int> ProbeComponents;                 // velocity component to measure 0 for x 1 for y
        std::vector<Point<DIM, real_number>> Offsets;
        std::vector<int> maxIters;

        std::vector<int> FixedProbeIndices_horizontal = {k0, kendWidth - 1};
        std::vector<int> FixedProbeIndices_vertical = {k0, kendHeight - 1};
        std::vector<real_number> FixedProbeValues_horizontal = {0.0, 0.0};
        std::vector<real_number> FixedProbeValues_vertical = {0.0, params.vw_top[0]};

        std::vector<std::vector<int>> FixedProbeIndices;
        std::vector<std::vector<real_number>> FixedProbeValues;

        if (params.SCENARIO == CAVITY)
        {
            ProbePoints.push_back(HalfChannel);
            ProbePoints.push_back(HalfHeight);

            ProbeComponents.push_back(0); // measure x velocity
            ProbeComponents.push_back(1); // measure y velocity

            Offsets.push_back(VerticalOffset);
            Offsets.push_back(HorizontalOffset);

            maxIters.push_back(kendHeight);
            maxIters.push_back(kendWidth);

            FixedProbeIndices.push_back(FixedProbeIndices_vertical);
            FixedProbeIndices.push_back(FixedProbeIndices_horizontal);
            FixedProbeValues.push_back(FixedProbeValues_vertical);
            FixedProbeValues.push_back(FixedProbeValues_horizontal);
        }

        for (unsigned int k = 0; k < ProbePoints.size(); k++)
        {
            // create probe object
            Ghost<DIM, real_number> gp(0);
            size_t bc_p[DIM] = {NON_PERIODIC, NON_PERIODIC};
            probe_particles vp_loc(0, domain, bc_p, gp, DEC_GRAN(512));

            if (v_cl.getProcessUnitID() == 0)
            {
                PlaceProbes(vp_loc, k0, maxIters[k], ProbePoints[k], Offsets[k], obstacle_ptr, FixedProbeIndices[k], FixedProbeValues[k]);
            }
            if (ProbeComponents[k] == 0) // measuring x velocity
            {
                openfpm::vector<std::string> names_p({"type", "vx"});
                vp_loc.setPropNames(names_p);
            }
            else if (ProbeComponents[k] == 1) // measuring y velocity
            {
                openfpm::vector<std::string> names_p({"type", "vy"});
                vp_loc.setPropNames(names_p);
            }
            std::pair<probe_particles, int> tmp = std::make_pair(vp_loc, ProbeComponents[k]);
            vp_vec.push_back(tmp);
            auxParams.probe_filenames.push_back("probes_" + std::to_string(k) + "_" + auxParams.filename);
        }
    }

    // Add the obstacle/walls as marker particles only on processor 0
    if (params.BC_TYPE == NEW_NO_SLIP && v_cl.getProcessUnitID() == 0)
    {
        real_number refine_factor = params.rf;
        real_number dx_wall_x = dp / refine_factor;
        int Nwall_x_bot = ceil(params.length[0] / dx_wall_x);
        int Nwall_x_top = ceil((params.length[0] + 9.0 * dp) / dx_wall_x);

        dx_wall_x = params.length[0] / Nwall_x_bot;
        Point<DIM, real_number> X_Offset = {dx_wall_x, 0.0};

        real_number dx_wall_y = dp / refine_factor;
        int Nwall_y = ceil(params.length[1] / dx_wall_y);
        dx_wall_y = params.length[1] / Nwall_y;
        Point<DIM, real_number> Y_Offset = {0.0, dx_wall_y};

        Point<DIM, real_number> LL_corner = {0.0, 0.0};
        Point<DIM, real_number> LR_corner = {params.length[0], 0.0};

        Point<DIM, real_number> UL_corner = {-4.0 * dp, params.length[1]};

        // Top And Bottom Walls
        AddFlatWallNewBC(vd, 0, Nwall_x_bot + 1, LL_corner, X_Offset, dx_wall_x, {0.0, 0.0}, params.vw_bottom, params, BOUNDARY, 0.0);
        AddFlatWallModNewBC(vd, 0, Nwall_x_top + 1, UL_corner, X_Offset, dx_wall_x, {0.0, 0.0}, params.vw_top, params, BOUNDARY, {0.0, -10.0}, 1.0, 0.0);

        // Left And Right Walls
        real_number r_cut_int = ceil(params.r_cut / dx_wall_y); // r_cut expressed in number of particles ( ceil to get an integer )
        // from 1 to N_wall_y - r_cut_int add normal wall
        AddFlatWallNewBC(vd, 1, Nwall_y - r_cut_int, LL_corner, Y_Offset, dx_wall_y, {0.0, 0.0}, {0.0, 0.0}, params, BOUNDARY, 0.0);
        AddFlatWallNewBC(vd, 1, Nwall_y - r_cut_int, LR_corner, Y_Offset, dx_wall_y, {0.0, 0.0}, {0.0, 0.0}, params, BOUNDARY, 0.0);
        // From N_wall_y - r_cut_int to N_wall_y add wall with prescribed normal and curvature
        AddFlatWallModNewBC(vd, Nwall_y - r_cut_int, Nwall_y + 1, LL_corner, Y_Offset, dx_wall_y, {0.0, 0.0}, {0.0, 0.0}, params, BOUNDARY, {10.0, 0.0}, 1.0, 0.0);
        AddFlatWallModNewBC(vd, Nwall_y - r_cut_int, Nwall_y + 1, LR_corner, Y_Offset, dx_wall_y, {0.0, 0.0}, {0.0, 0.0}, params, BOUNDARY, {-10.0, 0.0}, 1.0, 0.0);
    }

    // return an iterator to the fluid particles to add to vd
    auto fluid_it = DrawParticles::DrawBox(vd, sz, domain, fluid_box);

    // for each particle inside the fluid box ...
    while (fluid_it.isNext())
    {

        Point<DIM, real_number> iterator_position = fluid_it.get();

        // ... add a particle ...
        vd.add();
        vd.template getLastProp<vd0_type>() = FLUID;
        vd.template getLastProp<vd10_omega>() = 0.0;
        for (int xyz = 0; xyz < DIM; xyz++)
        {
            vd.template getLastProp<vd4_velocity>()[xyz] = 0.0;
            vd.template getLastProp<vd7_force_t>()[xyz] = 0.0;
        }

        // Set properties
        vd.template getLastProp<vd1_rho>() = params.rho0;
        vd.template getLastProp<vd2_pressure>() = 0.0;
        vd.template getLastProp<vd3_drho>() = 0.0;

        for (int xyz = 0; xyz < DIM; xyz++)
        {
            vd.getLastPos()[xyz] = iterator_position.get(xyz);
            vd.template getLastProp<vd5_velocity_t>()[xyz] = 0.0;
            vd.template getLastProp<vd8_normal>()[xyz] = 0.0;
        }

        vd.template getLastProp<vd9_volume>()[0] = dp;
        vd.template getLastProp<vd9_volume>()[1] = 0.0;
        vd.template getLastProp<vd9_volume>()[2] = 0.0;

        // next fluid particle
        ++fluid_it;
    }

    // Now place solid walls using iterators (only for OLD BC)

    if (params.BC_TYPE == NO_SLIP)
    {

        openfpm::vector<Box<DIM, real_number>> holes;
        holes.add(fluid_box);
        Box<DIM, real_number> hole_box = holes.get(0);
        auto bound_box = DrawParticles::DrawSkin(vd, sz, domain, holes, recipient);

        while (bound_box.isNext())
        {
            Point<DIM, real_number> position = bound_box.get();

            // periodic bc, with no boundary particles in y direction has a bug, it puts 3 extra particles outside in the y direction
            // When running on multiple cores, with this we check if particle is outside the recipient box
            // Another bug places boundary particles in the correct plane, but inside the fluid box;
            // if (bc[0] == PERIODIC && position.get(0) > dp / 2.0 && position.get(0) < length[0] - dp / 2.0)
            // {
            // 	++bound_box;
            // 	continue;
            // }

            if (!recipient.isInside((position)))
            {
                ++bound_box;
                continue;
            }
            if (hole_box.isInside(position))
            {
                ++bound_box;
                continue;
            }

            vd.add();

            vd.template getLastProp<vd0_type>() = BOUNDARY;
            vd.template getLastProp<vd1_rho>() = params.rho0;
            vd.template getLastProp<vd2_pressure>() = 0.0;
            vd.template getLastProp<vd3_drho>() = 0.0;

            for (int xyz = 0; xyz < DIM; xyz++)
            {
                vd.getLastPos()[xyz] = bound_box.get().get(xyz);
                vd.template getLastProp<vd6_force>()[xyz] = 0.0;
                vd.template getLastProp<vd7_force_t>()[xyz] = 0.0;
                vd.template getLastProp<vd5_velocity_t>()[xyz] = 0.0;
                vd.template getLastProp<vd8_normal>()[xyz] = 0.0;
                if (position.get(1) < dp / 4.0) // bottom wall
                {
                    vd.template getLastProp<vd4_velocity>()[xyz] = params.vw_bottom[xyz];
                }
                else if (position.get(1) > params.length[1] - dp / 4.0) // top wall
                {
                    vd.template getLastProp<vd4_velocity>()[xyz] = params.vw_top[xyz];
                }
            }

            vd.template getLastProp<vd9_volume>()[0] = dp;
            vd.template getLastProp<vd9_volume>()[1] = 0.0;
            vd.template getLastProp<vd9_volume>()[2] = 0.0;
            vd.template getLastProp<vd11_vorticity>() = 0.0;
            vd.template getLastProp<vd12_vel_red>() = 0.0;
            vd.template getLastProp<vd13_force_red_x>() = 0.0;
            vd.template getLastProp<vd14_force_red_y>() = 0.0;

            ++bound_box;
        }

        if (v_cl.getProcessUnitID() == 0)
        {
            // RIGHT CORNER
            const Point<DIM, real_number> Corner = {params.length[0] + 2.5 * dp,
                                                    params.length[1] + 0.5 * dp};

            // Manually place chunck of particles to fill periodicity of top wall
            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                {

                    Point<DIM, real_number> offset = {(j + 1) * dp, (i)*dp};
                    Point<DIM, real_number> position = Corner + offset;

                    vd.add();

                    vd.template getLastProp<vd0_type>() = BOUNDARY;
                    vd.template getLastProp<vd1_rho>() = params.rho0;
                    vd.template getLastProp<vd2_pressure>() = 0.0;
                    vd.template getLastProp<vd3_drho>() = 0.0;

                    for (int xyz = 0; xyz < DIM; xyz++)
                    {
                        vd.getLastPos()[xyz] = position.get(xyz);
                        vd.template getLastProp<vd6_force>()[xyz] = 0.0;
                        vd.template getLastProp<vd7_force_t>()[xyz] = 0.0;
                        vd.template getLastProp<vd5_velocity_t>()[xyz] = 0.0;
                        vd.template getLastProp<vd8_normal>()[xyz] = 0.0;

                        vd.template getLastProp<vd4_velocity>()[xyz] = params.vw_top[xyz];
                    }

                    vd.template getLastProp<vd9_volume>()[0] = dp;
                }
            }

            // LEFT CORNER
            const Point<DIM, real_number> CornerUL = {0.0 - 2.5 * dp,
                                                      params.length[1] + 0.5 * dp};

            // Manually place chunck of particles to fill periodicity of top wall
            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 4; j++)
                {

                    Point<DIM, real_number> offset = {-(j + 1) * dp, (i)*dp};
                    Point<DIM, real_number> position = CornerUL + offset;

                    vd.add();

                    vd.template getLastProp<vd0_type>() = BOUNDARY;
                    vd.template getLastProp<vd1_rho>() = params.rho0;
                    vd.template getLastProp<vd2_pressure>() = 0.0;
                    vd.template getLastProp<vd3_drho>() = 0.0;

                    for (int xyz = 0; xyz < DIM; xyz++)
                    {
                        vd.getLastPos()[xyz] = position.get(xyz);
                        vd.template getLastProp<vd6_force>()[xyz] = 0.0;
                        vd.template getLastProp<vd7_force_t>()[xyz] = 0.0;
                        vd.template getLastProp<vd5_velocity_t>()[xyz] = 0.0;
                        vd.template getLastProp<vd8_normal>()[xyz] = 0.0;

                        vd.template getLastProp<vd4_velocity>()[xyz] = params.vw_top[xyz];
                    }

                    vd.template getLastProp<vd9_volume>()[0] = dp;
                }
            }
        }
    }
}

void CreateParticleGeometrySphere(particles &vd, std::vector<std::pair<probe_particles, int>> &vp_vec, Obstacle *&obstacle_ptr, Parameters &params, AuxiliarParameters &auxParams)
{

    Vcluster<> &v_cl = create_vcluster();

    // Non periodic situation grid of 5 fluid particles and 3 boundary particles
    // We need a virtual grid of 5 + 2*(3+1) particles,
    // therefore the domain is discretized with 13 grid points,
    // when we use DrawParticles::DrawBox we will draw only the particles at the grid positons strictly inside the box,
    // the () repesent the recipient box, and the || represent the fluid box, we can see how this distribution places exactly 5 fluid particles inside and 3 boundary particles
    //           D-(-o--o--o-|-x--x--x--x--x--|-o-o-o-)-D
    // D: domain, o: boundary, x: fluid, --: dp distance
    // in a periodic situation we have the following
    // .....--x--x--D-|-x--x--x--x--x--|-D--x--x--......
    // therefore we need a grid of 5 + 2 particles, and the domain is discretized with 7 grid points

    // Size of the virtual cartesian grid that defines where to place the particles
    size_t sz[DIM];

    obstacle_ptr = new EmptyObstacle(params);
    real_number dp = params.dp;
    // Initialize obstacle in scenarios where needed
    if (params.SCENARIO == SPHERE)
        obstacle_ptr = new SphereObstacle(params.ObstacleBase, params.ObstacleCenter, params, params.ObstacleVelocity, params.ObstacleOmega, params.rf, true);

    real_number refine_factor = params.rf;

    // Now define the iterator boxes
    // We define the boxes in terms of offstes with respect to the fluid box that goes from 0 to length
    real_number offset_domain_left[DIM] = {0.0};
    real_number offset_domain_right[DIM] = {0.0};
    real_number offset_recipient[DIM] = {0.0};
    real_number offset_periodic_fluid[DIM] = {0.0};
    real_number offset_periodic_recipient[DIM] = {0.0};

    for (int xyz = 0; xyz < DIM; xyz++)
    {
        if (params.bc[xyz] == NON_PERIODIC) // non periodic, fluid covered by boundary
        {
            sz[xyz] = params.Nfluid[xyz] + 2 * (params.Nboundary[xyz] + 1);
            offset_domain_left[xyz] = (0.5 + params.Nboundary[xyz]) * dp;
            offset_domain_right[xyz] = (0.5 + params.Nboundary[xyz]) * dp;

            if (params.BC_TYPE == NEW_NO_SLIP) // Nboundary should only be 0 or 1 if we are using the new bc
                offset_recipient[xyz] = 0.1 * params.Nboundary[xyz] * dp;
            else if (params.BC_TYPE == NO_SLIP)
                offset_recipient[xyz] = params.Nboundary[xyz] * dp;
        }
        else // periodic, open ended
        {
            sz[xyz] = params.Nfluid[xyz] + 1;

            offset_domain_left[xyz] = 0.0;
            offset_domain_right[xyz] = dp;
            offset_periodic_fluid[xyz] = 0.75 * dp;
            offset_periodic_recipient[xyz] = 0.85 * dp;
        }
    }

    // Define the boxes
    Box<DIM, real_number> domain({-offset_domain_left[0],
                                  -offset_domain_left[1],
                                  -offset_domain_left[2]},
                                 {params.length[0] + offset_domain_right[0],
                                  params.length[1] + offset_domain_right[1],
                                  params.length[2] + offset_domain_right[2]});

    Box<DIM, real_number> fluid_box({0.0,
                                     0.0,
                                     0.0},
                                    {params.length[0] + offset_periodic_fluid[0],
                                     params.length[1] + offset_periodic_fluid[1],
                                     params.length[2] + offset_periodic_fluid[2]});

    Box<DIM, real_number> recipient({-offset_recipient[0],
                                     -offset_recipient[1],
                                     -offset_recipient[2]},
                                    {params.length[0] + offset_recipient[0] + offset_periodic_recipient[0],
                                     params.length[1] + offset_recipient[1] + offset_periodic_recipient[1],
                                     params.length[2] + offset_recipient[2] + offset_periodic_recipient[2]});

    // Will only be used in the new bc
    Box<DIM, real_number> recipient_hole({offset_recipient[0],
                                          offset_recipient[1],
                                          offset_recipient[2]},
                                         {params.length[0] - offset_recipient[0] + offset_periodic_fluid[0],
                                          params.length[1] - offset_recipient[1] + offset_periodic_fluid[1],
                                          params.length[2] - offset_recipient[2] + offset_periodic_fluid[2]});

    for (int xyz = 0; xyz < DIM; xyz++) // correct length in periodic case
    {
        if (params.bc[xyz] == PERIODIC)
            params.length[xyz] += dp;
    }

    // extended boundary around the domain, and the processor domain
    Ghost<DIM, real_number> g(params.r_cut);

    // create particle object
    particles vd_loc(0, domain, params.bc, g, DEC_GRAN(128));
    // vd is argument passed as reference we want to fill with particles
    vd = vd_loc;

    // Add the obstacle/walls as marker particles only on processor 0
    if (params.BC_TYPE == NEW_NO_SLIP && v_cl.getProcessUnitID() == 0)
    {
        // Add obstacle
        obstacle_ptr->AddObstacle(vd);

        // Add walls
        if (params.bc[0] == PERIODIC && params.bc[1] == NON_PERIODIC && params.bc[2] == NON_PERIODIC) // Channel like scenario
        {
            real_number dx_wall = dp / refine_factor;
            int Nwall_x = ceil(params.length[0] / dx_wall);
            int Nwall_y = ceil(params.length[1] / dx_wall);
            int Nwall_z = ceil(params.length[2] / dx_wall);

            dx_wall = params.length[0] / Nwall_x;
            Point<DIM, real_number> X_Offset = {dx_wall, 0.0, 0.0};
            Point<DIM, real_number> Y_Offset = {0.0, dx_wall, 0.0};
            Point<DIM, real_number> Z_Offset = {0.0, 0.0, dx_wall};

            Point<DIM, real_number> LL_corner = {0.0, 0.0, 0.0};
            Point<DIM, real_number> UR_corner = {0.0, params.length[1], params.length[2]};

            AddFlatWallModNewBC3D(vd, 0, Nwall_x, 0, Nwall_z + 1, LL_corner, X_Offset, Z_Offset, dx_wall, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, params, FREE_SLIP_BOUNDARY, {0.0, 0.0, 0.0}, 0.0, 0.0);
            AddFlatWallModNewBC3D(vd, 0, Nwall_x, 1, Nwall_y + 1, LL_corner, X_Offset, Y_Offset, dx_wall, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, params, FREE_SLIP_BOUNDARY, {0.0, 0.0, 0.0}, 0.0, 0.0);

            AddFlatWallModNewBC3D(vd, 0, Nwall_x, 0, Nwall_z, UR_corner, X_Offset, -Z_Offset, dx_wall, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, params, FREE_SLIP_BOUNDARY, {0.0, 0.0, 0.0}, 0.0, 0.0);
            AddFlatWallModNewBC3D(vd, 0, Nwall_x, 1, Nwall_y, UR_corner, X_Offset, -Y_Offset, dx_wall, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, params, FREE_SLIP_BOUNDARY, {0.0, 0.0, 0.0}, 0.0, 0.0);
        }
    }

    // return an iterator to the fluid particles to add to vd
    auto fluid_it = DrawParticles::DrawBox(vd, sz, domain, fluid_box);

    // for each particle inside the fluid box ...
    while (fluid_it.isNext())
    {

        Point<DIM, real_number> iterator_position = fluid_it.get();

        if ((*obstacle_ptr).isInside(iterator_position)) // if inside the obstacle region
        {
            if (params.BC_TYPE == NO_SLIP) // add particle but set it as boundary
            {
                // ... add a particle ...
                vd.add();
                vd.template getLastProp<vd0_type>() = OBSTACLE;
                vd.template getLastProp<vd10_omega>() = (*obstacle_ptr).AngularVelocity_;
                for (int xyz = 0; xyz < DIM; xyz++)
                {
                    vd.template getLastProp<vd4_velocity>()[xyz] = ((*obstacle_ptr).LinearVelocity_).get(xyz);
                    vd.template getLastProp<vd7_force_t>()[xyz] = ((*obstacle_ptr).Centre_).get(xyz);
                }
            }
            else if (params.BC_TYPE == NEW_NO_SLIP) // not add particle because already added
            {
                ++fluid_it;
                continue;
            }
        }
        else // if not inside obstacle at all just add fluid particles
        {
            // ... add a particle ...
            vd.add();
            vd.template getLastProp<vd0_type>() = FLUID;
            vd.template getLastProp<vd10_omega>() = 0.0;
            for (int xyz = 0; xyz < DIM; xyz++)
            {
                vd.template getLastProp<vd7_force_t>()[xyz] = 0.0;
            }
        }

        // Set properties
        vd.template getLastProp<vd1_rho>() = params.rho0;
        vd.template getLastProp<vd2_pressure>() = 0.0;
        vd.template getLastProp<vd3_drho>() = 0.0;

        for (int xyz = 0; xyz < DIM; xyz++)
        {
            vd.getLastPos()[xyz] = iterator_position.get(xyz);
            vd.template getLastProp<vd5_velocity_t>()[xyz] = 0.0;
            vd.template getLastProp<vd8_normal>()[xyz] = 0.0;
        }

        vd.template getLastProp<vd4_velocity>()[0] = 0.0;
        vd.template getLastProp<vd4_velocity>()[1] = 0.0;
        vd.template getLastProp<vd4_velocity>()[2] = 0.0;

        vd.template getLastProp<vd9_volume>()[0] = dp;
        vd.template getLastProp<vd9_volume>()[1] = 0.0;
        vd.template getLastProp<vd9_volume>()[2] = 0.0;

        // next fluid particle
        ++fluid_it;
    }

    // Now place solid walls using iterators (only for OLD BC)

    if (params.BC_TYPE == NO_SLIP)
    {

        openfpm::vector<Box<DIM, real_number>> holes;
        holes.add(fluid_box);
        Box<DIM, real_number> hole_box = holes.get(0);
        auto bound_box = DrawParticles::DrawSkin(vd, sz, domain, holes, recipient);

        if (params.bc[0] != PERIODIC || params.bc[1] != PERIODIC) // no walls in all periodic scenario
        {
            while (bound_box.isNext())
            {
                Point<DIM, real_number> position = bound_box.get();

                // periodic bc, with no boundary particles in y direction has a bug, it puts 3 extra particles outside in the y direction
                // When running on multiple cores, with this we check if particle is outside the recipient box
                // Another bug places boundary particles in the correct plane, but inside the fluid box;
                // if (bc[0] == PERIODIC && position.get(0) > dp / 2.0 && position.get(0) < length[0] - dp / 2.0)
                // {
                // 	++bound_box;
                // 	continue;
                // }

                if (!recipient.isInside((position)))
                {
                    ++bound_box;
                    continue;
                }
                if (hole_box.isInside(position))
                {
                    ++bound_box;
                    continue;
                }

                // if (params.BC_TYPE == NEW_NO_SLIP && (params.bc[0] == NON_PERIODIC && params.bc[1] == NON_PERIODIC))
                // {
                // 	// Check if x and z coordinates are multiples of dp, keep multiples, discard the rest
                // 	real_number remx = fmod(position.get(0), dp);
                // 	real_number remz = fmod(position.get(1), dp);
                // 	real_number tol = 0.5 * dp * 10e-2;

                // 	if (remx > tol && remx < dp - tol)
                // 	{
                // 		++bound_box;
                // 		continue;
                // 	}
                // 	if (remz > tol && remz < dp - tol)
                // 	{
                // 		++bound_box;
                // 		continue;
                // 	}
                // }
                vd.add();

                vd.template getLastProp<vd0_type>() = BOUNDARY;
                vd.template getLastProp<vd1_rho>() = params.rho0;
                vd.template getLastProp<vd2_pressure>() = 0.0;
                vd.template getLastProp<vd3_drho>() = 0.0;

                for (int xyz = 0; xyz < DIM; xyz++)
                {
                    vd.getLastPos()[xyz] = bound_box.get().get(xyz);
                    vd.template getLastProp<vd6_force>()[xyz] = 0.0;
                    vd.template getLastProp<vd7_force_t>()[xyz] = 0.0;
                    vd.template getLastProp<vd5_velocity_t>()[xyz] = 0.0;
                    vd.template getLastProp<vd8_normal>()[xyz] = 0.0;
                    if (position.get(1) < dp / 4.0) // bottom wall
                    {
                        vd.template getLastProp<vd4_velocity>()[xyz] = params.vw_bottom[xyz];
                    }
                    else if (position.get(1) > params.length[1] - dp / 4.0) // top wall
                    {
                        vd.template getLastProp<vd4_velocity>()[xyz] = params.vw_top[xyz];
                    }
                }

                vd.template getLastProp<vd9_volume>()[0] = dp;

                ++bound_box;
            }
        }
    }
}
int FindSector(real_number th, std::vector<double> angles)
{
    // angles contain angles defining circular sectors, in increasing order
    // need to find to which sector the angle th belongs to
    // return the index beginning from 0

    // normalzie th to [0, 2pi]
    while (th < 0)
        th += 2 * M_PI;
    while (th >= 2 * M_PI)
        th -= 2 * M_PI;
    // normalize angles to [0, 2pi]
    for (int i = 0; i < angles.size(); i++)
    {
        while (angles[i] < 0)
            angles[i] += 2 * M_PI;
        while (angles[i] >= 2 * M_PI)
            angles[i] -= 2 * M_PI;
    }

    int nsectors = angles.size();
    int sector = -1;
    // Handle the last sector wrapping around to the first
    for (int i = 0; i < nsectors; i++)
    {
        double start_angle = angles[i];
        double end_angle = (i + 1 < nsectors) ? angles[i + 1] : angles[0] + 2 * M_PI;

        if (th >= start_angle && th < end_angle)
        {
            sector = i;
            break;
        }
    }

    // printf("th = %f, sector = %d\n", th * 180.0 / M_PI, sector);

    return sector;
}
void CreateParticleGeometryFlower(particles &vd, std::vector<std::pair<probe_particles, int>> &vp_vec, Obstacle *&obstacle_ptr, Parameters &params, AuxiliarParameters &auxParams)
{

    Vcluster<> &v_cl = create_vcluster();

    // Non periodic situation grid of 5 fluid particles and 3 boundary particles
    // We need a virtual grid of 5 + 2*(3+1) particles,
    // therefore the domain is discretized with 13 grid points,
    // when we use DrawParticles::DrawBox we will draw only the particles at the grid positons strictly inside the box,
    // the () repesent the recipient box, and the || represent the fluid box, we can see how this distribution places exactly 5 fluid particles inside and 3 boundary particles
    //           D-(-o--o--o-|-x--x--x--x--x--|-o-o-o-)-D
    // D: domain, o: boundary, x: fluid, --: dp distance
    // in a periodic situation we have the following
    // .....--x--x--D-|-x--x--x--x--x--|-D--x--x--......
    // therefore we need a grid of 5 + 2 particles, and the domain is discretized with 7 grid points

    // Size of the virtual cartesian grid that defines where to place the particles
    size_t sz[DIM];

    obstacle_ptr = new EmptyObstacle(params);
    Obstacle *obstacle_ptr2 = new EmptyObstacle(params);

    real_number dp = params.dp;

    // obstacle_ptr = new CurveObstacle(params.flowerA, params.flowerB, params.flowerK, params.flowerM, params.ObstacleCenter, params, params.ObstacleVelocity, params.ObstacleOmega, params.rf);
    obstacle_ptr = new EpiCycloid(params.flowerPropRadius * params.flowerA, params.flowerNlobes, params.ObstacleCenter, params, params.ObstacleVelocity, 0.0, params.rf);
    obstacle_ptr2 = new HipoCycloid(params.flowerPropRadius * params.flowerA / params.flowerPropRadius2, params.flowerNlobes, params.ObstacleCenter, params, params.ObstacleVelocity, 0.0, params.rf);

    real_number rotor_pos_r = params.flowerPropRadius * (5.0 / 3.0) * params.flowerA / 2.0;

    std::vector<Obstacle *> rotors;
    int nlobes = (int)params.flowerNlobes;
    std::vector<real_number> ths;
    for (int kk = 0; kk < nlobes; kk++)
    {
        int sign = 1;
        if (params.flowerAlternate == 1)
            sign = std::pow(-1, kk);

        real_number th = M_PI / params.flowerNlobes + 2 * M_PI * kk / params.flowerNlobes;
        ths.push_back(2 * M_PI * kk / params.flowerNlobes);
        // printf("th = %f\n", 2 * M_PI * kk / params.flowerNlobes);
        Point<DIM, real_number> rotor_pos = params.ObstacleCenter;
        Point<DIM, real_number> offset = {rotor_pos_r * cos(th), rotor_pos_r * sin(th)};
        rotor_pos = rotor_pos + offset;
        Obstacle *tmp = new CurveObstacle(params.flowerA, params.flowerB, params.flowerK, params.flowerM, rotor_pos, params, params.ObstacleVelocity, sign * params.ObstacleOmega, params.rf);
        rotors.push_back(tmp);
    }

    // Now define the iterator boxes
    // We define the boxes in terms of offstes with respect to the fluid box that goes from 0 to length
    real_number offset_domain_left[DIM] = {0.0};
    real_number offset_domain_right[DIM] = {0.0};
    real_number offset_recipient[DIM] = {0.0};
    real_number offset_periodic_fluid[DIM] = {0.0};
    real_number offset_periodic_recipient[DIM] = {0.0};

    for (int xyz = 0; xyz < DIM; xyz++)
    {
        if (params.bc[xyz] == NON_PERIODIC) // non periodic, fluid covered by boundary
        {
            sz[xyz] = params.Nfluid[xyz] + 2 * (params.Nboundary[xyz] + 1);
            offset_domain_left[xyz] = (0.5 + params.Nboundary[xyz]) * dp;
            offset_domain_right[xyz] = (0.5 + params.Nboundary[xyz]) * dp;

            if (params.BC_TYPE == NEW_NO_SLIP) // Nboundary should only be 0 or 1 if we are using the new bc
                offset_recipient[xyz] = 0.1 * params.Nboundary[xyz] * dp;
            else if (params.BC_TYPE == NO_SLIP)
                offset_recipient[xyz] = params.Nboundary[xyz] * dp;
        }
        else // periodic, open ended
        {
            sz[xyz] = params.Nfluid[xyz] + 1;

            offset_domain_left[xyz] = 0.0;
            offset_domain_right[xyz] = dp;
            offset_periodic_fluid[xyz] = 0.75 * dp;
            offset_periodic_recipient[xyz] = 0.85 * dp;
        }
    }

    // Define the boxes
    Box<DIM, real_number> domain({-offset_domain_left[0],
                                  -offset_domain_left[1]},
                                 {params.length[0] + offset_domain_right[0],
                                  params.length[1] + offset_domain_right[1]});

    Box<DIM, real_number> fluid_box({0.0,
                                     0.0},
                                    {params.length[0] + offset_periodic_fluid[0],
                                     params.length[1] + offset_periodic_fluid[1]});

    Box<DIM, real_number> recipient({-offset_recipient[0],
                                     -offset_recipient[1]},
                                    {params.length[0] + offset_recipient[0] + offset_periodic_recipient[0],
                                     params.length[1] + offset_recipient[1] + offset_periodic_recipient[1]});

    // Will only be used in the new bc
    // Box<DIM, real_number> recipient_hole({offset_recipient[0],
    //                                       offset_recipient[1]},
    //                                      {params.length[0] - offset_recipient[0] + offset_periodic_fluid[0],
    //                                       params.length[1] - offset_recipient[1] + offset_periodic_fluid[1]});

    for (int xyz = 0; xyz < DIM; xyz++) // correct length in periodic case
    {
        if (params.bc[xyz] == PERIODIC)
            params.length[xyz] += dp;
    }

    // extended boundary around the domain, and the processor domain
    Ghost<DIM, real_number> g(params.r_cut);

    // create particle object
    particles vd_loc(0, domain, params.bc, g, DEC_GRAN(128));
    // vd is argument passed as reference we want to fill with particles
    vd = vd_loc;

    // For the new bc we place the wall/obstacle particles manually
    // Add the obstacle/walls as marker particles only on processor 0
    if (params.BC_TYPE == NEW_NO_SLIP && v_cl.getProcessUnitID() == 0)
    {
        // Add obstacle
        obstacle_ptr->AddObstacle(vd);
        obstacle_ptr2->AddObstacle(vd);

        for (auto &rotor : rotors)
        {
            rotor->AddObstacle(vd);
        }
    }

    // return an iterator to the fluid particles to add to vd
    auto fluid_it = DrawParticles::DrawBox(vd, sz, domain, fluid_box);

    // for each particle inside the fluid box ...
    if (params.BC_TYPE == NEW_NO_SLIP)
    {
        while (fluid_it.isNext())
        {

            Point<DIM, real_number> iterator_position = fluid_it.get();

            if (IsInsideAll(rotors, iterator_position) || !(*obstacle_ptr).isInside(iterator_position) || (*obstacle_ptr2).isInside(iterator_position)) // if not inside the obstacles, or outside the containing obstacle
            {
                ++fluid_it;
                continue;
            }
            else // if not inside obstacle at all just add fluid particles
            {
                // ... add a particle ...
                vd.add();
                vd.template getLastProp<vd0_type>() = FLUID;
                vd.template getLastProp<vd10_omega>() = 0.0;
                for (int xyz = 0; xyz < DIM; xyz++)
                {
                    vd.template getLastProp<vd4_velocity>()[xyz] = 0.0;
                    vd.template getLastProp<vd7_force_t>()[xyz] = 0.0;
                }
            }

            // Set properties common to fluid and obstacle particles
            vd.template getLastProp<vd1_rho>() = params.rho0;
            vd.template getLastProp<vd2_pressure>() = 0.0;
            vd.template getLastProp<vd3_drho>() = 0.0;

            for (int xyz = 0; xyz < DIM; xyz++)
            {
                vd.getLastPos()[xyz] = iterator_position.get(xyz);
                vd.template getLastProp<vd5_velocity_t>()[xyz] = 0.0;
                vd.template getLastProp<vd8_normal>()[xyz] = 0.0;
                vd.template getLastProp<vd6_force>()[xyz] = 0.0;
            }

            // fluid particles carry dp in volume[0], it is useful to read it in postprocessing
            // this array is unused for fluid particles anyway
            vd.template getLastProp<vd9_volume>()[0] = dp;

            real_number theta = atan2(iterator_position.get(1) - params.ObstacleCenter[1], iterator_position.get(0) - params.ObstacleCenter[0]);

            int idx = FindSector(theta, ths);

            // if ((iterator_position.get(0) > params.ObstacleCenter[0]) && (iterator_position.get(1) > params.ObstacleCenter[1]))
            //     vd.template getLastProp<vd9_volume>()[1] = 100.0;
            // else if ((iterator_position.get(0) < params.ObstacleCenter[0]) && (iterator_position.get(1) > params.ObstacleCenter[1]))
            //     vd.template getLastProp<vd9_volume>()[1] = 200.0;
            // else if ((iterator_position.get(0) < params.ObstacleCenter[0]) && (iterator_position.get(1) < params.ObstacleCenter[1]))
            //     vd.template getLastProp<vd9_volume>()[1] = 300.0;
            // else if ((iterator_position.get(0) > params.ObstacleCenter[0]) && (iterator_position.get(1) < params.ObstacleCenter[1]))
            //     vd.template getLastProp<vd9_volume>()[1] = 400.0;
            // else
            //     vd.template getLastProp<vd9_volume>()[1] = 0.0;

            vd.template getLastProp<vd9_volume>()[1] = 100.0 * idx;

            vd.template getLastProp<vd9_volume>()[2] = 0.0;

            vd.template getLastProp<vd11_vorticity>() = 0.0;
            vd.template getLastProp<vd12_vel_red>() = 0.0;
            vd.template getLastProp<vd13_force_red_x>() = 0.0;
            vd.template getLastProp<vd14_force_red_y>() = 0.0;

            // next fluid particle
            ++fluid_it;
        }
    }
    else
    {
        while (fluid_it.isNext())
        {

            Point<DIM, real_number> iterator_position = fluid_it.get();

            if (IsInsideAll(rotors, iterator_position) || !(*obstacle_ptr).isInside(iterator_position) || (*obstacle_ptr2).isInside(iterator_position)) // if not inside the obstacles, or outside the containing obstacle
            {
                // ... add a particle ...
                vd.add();
                vd.template getLastProp<vd0_type>() = OBSTACLE;
                vd.template getLastProp<vd10_omega>() = 0.0;

                for (int xyz = 0; xyz < DIM; xyz++)
                {
                    vd.template getLastProp<vd4_velocity>()[xyz] = 0.0;
                    vd.template getLastProp<vd7_force_t>()[xyz] = 0.0;
                }

                for (auto rotor : rotors)
                {
                    if (rotor->isInside(iterator_position))
                    {
                        vd.template getLastProp<vd10_omega>() = rotor->AngularVelocity_;
                        for (int xyz = 0; xyz < DIM; xyz++)
                        {
                            vd.template getLastProp<vd4_velocity>()[xyz] = (rotor->LinearVelocity_).get(xyz);
                            vd.template getLastProp<vd7_force_t>()[xyz] = (rotor->Centre_).get(xyz);
                        }
                    }
                }
            }
            else // if not inside obstacle at all just add fluid particles
            {
                // ... add a particle ...
                vd.add();
                vd.template getLastProp<vd0_type>() = FLUID;
                vd.template getLastProp<vd10_omega>() = 0.0;
                for (int xyz = 0; xyz < DIM; xyz++)
                {
                    vd.template getLastProp<vd4_velocity>()[xyz] = 0.0;
                    vd.template getLastProp<vd7_force_t>()[xyz] = 0.0;
                }
            }

            // Set properties common to fluid and obstacle particles
            vd.template getLastProp<vd1_rho>() = params.rho0;
            vd.template getLastProp<vd2_pressure>() = 0.0;
            vd.template getLastProp<vd3_drho>() = 0.0;

            for (int xyz = 0; xyz < DIM; xyz++)
            {
                vd.getLastPos()[xyz] = iterator_position.get(xyz);
                vd.template getLastProp<vd5_velocity_t>()[xyz] = 0.0;
                vd.template getLastProp<vd8_normal>()[xyz] = 0.0;
                vd.template getLastProp<vd6_force>()[xyz] = 0.0;
            }

            // fluid particles carry dp in volume[0], it is useful to read it in postprocessing
            // this array is unused for fluid particles anyway
            vd.template getLastProp<vd9_volume>()[0] = dp;

            real_number theta = atan2(iterator_position.get(1) - params.ObstacleCenter[1], iterator_position.get(0) - params.ObstacleCenter[0]);

            int idx = FindSector(theta, ths);

            // if ((iterator_position.get(0) > params.ObstacleCenter[0]) && (iterator_position.get(1) > params.ObstacleCenter[1]))
            //     vd.template getLastProp<vd9_volume>()[1] = 100.0;
            // else if ((iterator_position.get(0) < params.ObstacleCenter[0]) && (iterator_position.get(1) > params.ObstacleCenter[1]))
            //     vd.template getLastProp<vd9_volume>()[1] = 200.0;
            // else if ((iterator_position.get(0) < params.ObstacleCenter[0]) && (iterator_position.get(1) < params.ObstacleCenter[1]))
            //     vd.template getLastProp<vd9_volume>()[1] = 300.0;
            // else if ((iterator_position.get(0) > params.ObstacleCenter[0]) && (iterator_position.get(1) < params.ObstacleCenter[1]))
            //     vd.template getLastProp<vd9_volume>()[1] = 400.0;
            // else
            //     vd.template getLastProp<vd9_volume>()[1] = 0.0;

            if (vd.template getLastProp<vd0_type>() == FLUID)
                vd.template getLastProp<vd9_volume>()[1] = 100.0 * idx;

            vd.template getLastProp<vd9_volume>()[2] = 0.0;

            vd.template getLastProp<vd11_vorticity>() = 0.0;
            vd.template getLastProp<vd12_vel_red>() = 0.0;
            vd.template getLastProp<vd13_force_red_x>() = 0.0;
            vd.template getLastProp<vd14_force_red_y>() = 0.0;

            // next fluid particle
            ++fluid_it;
        }
    }
}
