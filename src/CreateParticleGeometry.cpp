#include "CreateParticleGeometry.hpp"

void CreateParticleGeometry(particles &vd, std::vector<std::pair<probe_particles, int>> &vp_vec, Vcluster<> &v_cl, Obstacle *&obstacle_ptr, Parameters &params)
{

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
    double dp = params.dp;
    // Initialize obstacle in scenarios where needed
    if (params.SCENARIO == CYLINDER_ARRAY)
        obstacle_ptr = new CylinderObstacle(params.LengthScale, params.ObstacleCenter, params, params.ObstacleVelocity, params.ObstacleOmega, params.rf);
    else if (params.SCENARIO == CYLINDER_LATTICE)
        obstacle_ptr = new CylinderObstacle(params.LengthScale, params.ObstacleCenter, params, params.ObstacleVelocity, params.ObstacleOmega, params.rf);
    else if (params.SCENARIO == SQUARE)
        obstacle_ptr = new RectangleObstacle(params.ObstacleCenter, params, params.ObstacleBase, params.ObstacleHeight, params.ObstacleVelocity, params.ObstacleOmega, params.rf);
    else if (params.SCENARIO == TRIANGLE)
        obstacle_ptr = new TriangleObstacle(params.ObstacleCenter, params, params.ObstacleBase, params.ObstacleHeight, params.ObstacleVelocity, params.ObstacleOmega, params.rf);
    else if (params.SCENARIO == TRIANGLE_EQUILATERAL)
        obstacle_ptr = new TriangleEqui(params.ObstacleCenter, params, params.ObstacleBase, params.ObstacleVelocity, params.ObstacleOmega, params.rf);
    else if (params.SCENARIO == MOVING_OBSTACLE)
        obstacle_ptr = new TriangleEqui(params.ObstacleCenter, params, params.ObstacleBase, params.ObstacleVelocity, params.ObstacleOmega, params.rf);
    else if (params.SCENARIO == ELLIPSE)
        obstacle_ptr = new EllipticObstacle(params.ObstacleBase, params.ObstacleHeight, params.ObstacleTilt, params.ObstacleCenter, params, params.ObstacleVelocity, params.ObstacleOmega, params.rf);
    double refine_factor = params.rf;

    // Now define the iterator boxes
    // We define the boxes in terms of offstes with respect to the fluid box that goes from 0 to length
    double offset_domain_left[DIM] = {0.0};
    double offset_domain_right[DIM] = {0.0};
    double offset_recipient[DIM] = {0.0};
    double offset_periodic_fluid[DIM] = {0.0};
    double offset_periodic_recipient[DIM] = {0.0};

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
    Box<DIM, double> domain({-offset_domain_left[0],
                             -offset_domain_left[1]},
                            {params.length[0] + offset_domain_right[0],
                             params.length[1] + offset_domain_right[1]});

    Box<DIM, double> fluid_box({0.0,
                                0.0},
                               {params.length[0] + offset_periodic_fluid[0],
                                params.length[1] + offset_periodic_fluid[1]});

    Box<DIM, double> recipient({-offset_recipient[0],
                                -offset_recipient[1]},
                               {params.length[0] + offset_recipient[0] + offset_periodic_recipient[0],
                                params.length[1] + offset_recipient[1] + offset_periodic_recipient[1]});

    // Will only be used in the new bc
    Box<DIM, double> recipient_hole({offset_recipient[0],
                                     offset_recipient[1]},
                                    {params.length[0] - offset_recipient[0] + offset_periodic_fluid[0],
                                     params.length[1] - offset_recipient[1] + offset_periodic_fluid[1]});

    for (int xyz = 0; xyz < DIM; xyz++) // correct length in periodic case
    {
        if (params.bc[xyz] == PERIODIC)
            params.length[xyz] += dp;
    }

    // extended boundary around the domain, and the processor domain
    Ghost<DIM, double> g(params.r_threshold);

    // create particle object
    particles vd_loc(0, domain, params.bc, g, DEC_GRAN(512));
    // vd is argument passed as reference we want to fill with particles
    vd = vd_loc;

    // Write constants on file
    // WriteParameters(v_cl, params);

    // place probes
    if (params.PROBES_ENABLED)
    {
        // we want to place probes  in a vertical line at this locations
        Point<DIM, double> EndChannel = {params.length[0], 0.0};
        Point<DIM, double> HalfChannel = {params.length[0] / 2.0, 0.0};
        Point<DIM, double> HalfHeight = {0.0, params.length[1] / 2.0};
        Point<DIM, double> VerticalOffset = {0.0, dp};
        Point<DIM, double> HorizontalOffset = {dp, 0.0};
        int k0 = 0;
        int kendHeight = params.Nfluid[1] + 1;
        int kendWidth = params.Nfluid[0] + 1;

        std::vector<Point<DIM, double>> ProbePoints; // start points for the PlaceProbes function
        std::vector<int> ProbeComponents;            // component to measure 0 for x 1 for y
        std::vector<Point<DIM, double>> Offsets;
        std::vector<int> maxIters;

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
        }
        else
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

        for (unsigned int k = 0; k < ProbePoints.size(); k++)
        {
            // create probe object
            Ghost<DIM, double> gp(0);
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
                PlaceProbes(vp_loc, k0, maxIters[k], ProbePoints[k], Offsets[k]);
            }
            std::pair<probe_particles, int> tmp = std::make_pair(vp_loc, ProbeComponents[k]);
            vp_vec.push_back(tmp);
            params.probe_filenames.push_back("probes_" + std::to_string(k) + "_" + params.filename);
        }
    }

    // Add the obstacle/walls as marker particles only on processor 0
    if (params.BC_TYPE == NEW_NO_SLIP && v_cl.getProcessUnitID() == 0)
    {
        // Add obstacle
        obstacle_ptr->AddObstacle(vd);

        // Add walls
        if (params.bc[0] == PERIODIC && params.bc[1] == NON_PERIODIC) // Channel like scenario
        {
            double dx_wall = dp / refine_factor;
            int Nwall = ceil(params.length[0] / dx_wall);
            dx_wall = params.length[0] / Nwall;
            Point<DIM, double> X_Offset = {dx_wall, 0.0};

            Point<DIM, double> LL_corner = {0.0, 0.0};
            Point<DIM, double> UL_corner = {0.0, params.length[1]};
            // Top And Bottom Walls
            AddFlatWallNewBC(vd, 0, Nwall, LL_corner, X_Offset, dx_wall, {0.0, 0.0}, params.vw_bottom, params, 0.0);
            AddFlatWallNewBC(vd, 0, Nwall, UL_corner, X_Offset, dx_wall, {0.0, 0.0}, params.vw_top, params, 0.0);
        }
        else if (params.bc[0] == NON_PERIODIC && params.bc[1] == NON_PERIODIC) // Box like scenario
        {
            double dx_wall_x = dp / refine_factor;
            int Nwall_x = ceil(params.length[0] / dx_wall_x);
            dx_wall_x = params.length[0] / Nwall_x;
            Point<DIM, double> X_Offset = {dx_wall_x, 0.0};

            double dx_wall_y = dp / refine_factor;
            int Nwall_y = ceil(params.length[1] / dx_wall_y);
            dx_wall_y = params.length[1] / Nwall_y;
            Point<DIM, double> Y_Offset = {0.0, dx_wall_y};

            Point<DIM, double> LL_corner = {0.0, 0.0};
            Point<DIM, double> LR_corner = {params.length[0], 0.0};

            Point<DIM, double> UL_corner = {0.0, params.length[1]};

            // Top And Bottom Walls
            AddFlatWallNewBC(vd, 0, Nwall_x + 1, LL_corner, X_Offset, dx_wall_x, {0.0, 0.0}, params.vw_bottom, params, 0.0);
            AddFlatWallNewBC(vd, 0, Nwall_x + 1, UL_corner, X_Offset, dx_wall_x, {0.0, 0.0}, params.vw_top, params, 0.0);

            // Left And Right Walls
            AddFlatWallNewBC(vd, 1, Nwall_y, LL_corner, Y_Offset, dx_wall_y, {0.0, 0.0}, {0.0, 0.0}, params, 0.0);
            AddFlatWallNewBC(vd, 1, Nwall_y, LR_corner, Y_Offset, dx_wall_y, {0.0, 0.0}, {0.0, 0.0}, params, 0.0);
        }
    }

    // return an iterator to the fluid particles to add to vd
    auto fluid_it = DrawParticles::DrawBox(vd, sz, domain, fluid_box);

    // for each particle inside the fluid box ...
    while (fluid_it.isNext())
    {

        Point<DIM, double> iterator_position = fluid_it.get();

        if ((*obstacle_ptr).isInside(iterator_position)) // if inside the obstacle region
        {
            if (params.BC_TYPE == NO_SLIP) // add particle but set it as boundary
            {
                // ... add a particle ...
                vd.add();
                vd.template getLastProp<type>() = BOUNDARY;
                vd.template getLastProp<vd_omega>() = (*obstacle_ptr).AngularVelocity_;
                for (int xyz = 0; xyz < DIM; xyz++)
                {
                    vd.template getLastProp<velocity>()[xyz] = ((*obstacle_ptr).LinearVelocity_).get(xyz);
                    vd.template getLastProp<force_transport>()[xyz] = ((*obstacle_ptr).Centre_).get(xyz);
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
            vd.template getLastProp<type>() = FLUID;
            vd.template getLastProp<vd_omega>() = 0.0;
            for (int xyz = 0; xyz < DIM; xyz++)
            {
                vd.template getLastProp<velocity>()[xyz] = 0.0;
                vd.template getLastProp<force_transport>()[xyz] = 0.0;
            }
        }

        // for validating the vorticity calculation, vorticity should be 2w
        // double w = 5.0;
        // if (params.SCENARIO == HYDROSTATIC)
        // {
        //     vd.template getLastProp<velocity>()[0] = -w * iterator_position.get(1);
        //     vd.template getLastProp<velocity>()[1] = w * iterator_position.get(0);
        // }

        // Set properties
        vd.template getLastProp<rho>() = params.rho_zero;
        vd.template getLastProp<pressure>() = 0.0;
        vd.template getLastProp<drho>() = 0.0;

        for (int xyz = 0; xyz < DIM; xyz++)
        {
            vd.getLastPos()[xyz] = iterator_position.get(xyz);
            vd.template getLastProp<v_transport>()[xyz] = 0.0;
            vd.template getLastProp<normal_vector>()[xyz] = 0.0;
        }

        vd.template getLastProp<curvature_boundary>() = 0.0;
        vd.template getLastProp<arc_length>() = dp;

        // next fluid particle
        ++fluid_it;
    }

    // Now place solid walls using iterators (only for OLD BC)

    if (params.BC_TYPE == NO_SLIP)
    {

        openfpm::vector<Box<DIM, double>> holes;
        holes.add(fluid_box);
        Box<DIM, double> hole_box = holes.get(0);
        auto bound_box = DrawParticles::DrawSkin(vd, sz, domain, holes, recipient);

        if (params.bc[0] != PERIODIC || params.bc[1] != PERIODIC) // no walls in all periodic scenario
        {
            while (bound_box.isNext())
            {
                Point<DIM, double> position = bound_box.get();

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
                // 	double remx = fmod(position.get(0), dp);
                // 	double remz = fmod(position.get(1), dp);
                // 	double tol = 0.5 * dp * 10e-2;

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

                vd.template getLastProp<type>() = BOUNDARY;
                vd.template getLastProp<rho>() = params.rho_zero;
                vd.template getLastProp<pressure>() = 0.0;
                vd.template getLastProp<drho>() = 0.0;

                for (int xyz = 0; xyz < DIM; xyz++)
                {
                    vd.getLastPos()[xyz] = bound_box.get().get(xyz);
                    vd.template getLastProp<force>()[xyz] = 0.0;
                    vd.template getLastProp<force_transport>()[xyz] = 0.0;
                    vd.template getLastProp<v_transport>()[xyz] = 0.0;
                    vd.template getLastProp<normal_vector>()[xyz] = 0.0;
                    if (position.get(1) < dp / 4.0) // bottom wall
                    {
                        vd.template getLastProp<velocity>()[xyz] = params.vw_bottom.get(xyz);
                    }
                    else if (position.get(1) > params.length[1] - dp / 4.0) // top wall
                    {
                        vd.template getLastProp<velocity>()[xyz] = params.vw_top.get(xyz);
                    }
                }

                vd.template getLastProp<curvature_boundary>() = 0.0;
                vd.template getLastProp<arc_length>() = dp;

                ++bound_box;
            }
        }
    }
}
void CreateParticleGeometryTaylorCouette(particles &vd, std::vector<std::pair<probe_particles, int>> &vp_vec, Vcluster<> &v_cl, Obstacle *&obstacle_ptr, Parameters params)
{
    // Size of the virtual cartesian grid that defines where to place the particles
    size_t sz[DIM];

    double dp = params.dp;
    // We define the boxes in terms of offstes with respect to the fluid box that goes from 0 to length
    double offset_domain_left[DIM] = {0.0};
    double offset_domain_right[DIM] = {0.0};

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

    double Rin = params.Rin;
    double Rout = params.Rout;
    double Win = params.Win;
    double Wout = params.Wout;

    double a_tc = -((Rout * Rout * Rin * Rin) / (Rout * Rout - Rin * Rin)) * (Wout - Win);
    double b_tc = (Wout * Rout * Rout - Win * Rin * Rin) / (Rout * Rout - Rin * Rin);

    size_t Nbound = (params.BC_TYPE == NEW_NO_SLIP) ? 1 : 3;

    for (int dim = 0; dim < DIM; dim++)
    {
        params.length[dim] = dp * params.Nfluid[dim];
        sz[dim] = params.Nfluid[dim] + 2 * (Nbound + 1);
        offset_domain_left[dim] = (0.5 + Nbound) * dp;
        offset_domain_right[dim] = (0.5 + Nbound) * dp;
    }

    // Define the boxes
    Box<DIM, double> domain({-params.length[0] / 2.0 - offset_domain_left[0],
                             -params.length[1] / 2.0 - offset_domain_left[1]},
                            {params.length[0] / 2.0 + offset_domain_right[0],
                             params.length[1] / 2.0 + offset_domain_right[1]});
    Box<DIM, double> fluid_box({-params.length[0] / 2.0,
                                -params.length[1] / 2.0},
                               {params.length[0] / 2.0,
                                params.length[1] / 2.0});

    // extended boundary around the domain, and the processor domain
    Ghost<DIM, double> g(params.r_threshold);

    // create particle object
    particles vd_loc(0, domain, params.bc, g, DEC_GRAN(512));
    vd = vd_loc;

    // Write constants on file
    double rf = params.rf;

    // Set cylindrical object parameters
    Point<DIM, double> CylinderCentre = {0.0, 0.0};

    obstacle_ptr = new EmptyObstacle(params);

    const Point<DIM, double> vel = {0.0, 0.0};

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

    Box<DIM, double> fluid_box_aux({-Rout - 3.5 * dp,
                                    -Rout - 3.5 * dp},
                                   {Rout + 3.5 * dp,
                                    Rout + 3.5 * dp});

    // Outer Cylinder boundary particles
    auto out_it = DrawParticles::DrawBox(vd, sz, domain, fluid_box_aux);
    if (params.BC_TYPE == NO_SLIP)
    {
        while (out_it.isNext())
        {

            Point<DIM, double> iterator_position = out_it.get();
            if (!(*obstacle_ptr_out).isOutside(iterator_position) && (*obstacle_ptr_out_aux).isInside(iterator_position)) // if outside the outer cylinder and inside outer cylinder aux
            {
                if (params.BC_TYPE == NO_SLIP)
                {
                    vd.add();
                    // Set properties
                    vd.template getLastProp<type>() = BOUNDARY;
                    vd.template getLastProp<vd_omega>() = (*obstacle_ptr_out).AngularVelocity_;
                    for (int xyz = 0; xyz < DIM; xyz++)
                    {
                        vd.template getLastProp<velocity>()[xyz] = ((*obstacle_ptr_out).LinearVelocity_).get(xyz);
                        vd.template getLastProp<force_transport>()[xyz] = ((*obstacle_ptr_out).Centre_).get(xyz);
                    }
                    vd.template getLastProp<rho>() = params.rho_zero;
                    vd.template getLastProp<pressure>() = 0.0;
                    vd.template getLastProp<drho>() = 0.0;

                    for (int xyz = 0; xyz < DIM; xyz++)
                    {
                        vd.getLastPos()[xyz] = iterator_position.get(xyz);
                        vd.template getLastProp<v_transport>()[xyz] = 0.0;
                        vd.template getLastProp<normal_vector>()[xyz] = 0.0;
                    }

                    vd.template getLastProp<curvature_boundary>() = 0.0;
                    vd.template getLastProp<arc_length>() = dp;

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

        Point<DIM, double> iterator_position = fluid_it.get();
        if ((*obstacle_ptr_out).isOutside(iterator_position)) // if inside the outer cylinder
        {
            if ((*obstacle_ptr_in).isInside(iterator_position)) // if inside the inner cylinder region
            {
                if (!(*obstacle_ptr_in_aux).isOutside(iterator_position))
                {
                    if (params.BC_TYPE == NO_SLIP) // add particle but set it as boundary
                    {
                        // ... add a particle ...
                        vd.add();
                        vd.template getLastProp<type>() = BOUNDARY;
                        vd.template getLastProp<vd_omega>() = (*obstacle_ptr_in).AngularVelocity_;
                        for (int xyz = 0; xyz < DIM; xyz++)
                        {
                            vd.template getLastProp<velocity>()[xyz] = ((*obstacle_ptr_in).LinearVelocity_).get(xyz);
                            vd.template getLastProp<force_transport>()[xyz] = ((*obstacle_ptr_in).Centre_).get(xyz);
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
                vd.template getLastProp<type>() = FLUID;
                vd.template getLastProp<vd_omega>() = 0.0;

                double r = iterator_position.get(0) * iterator_position.get(0) + iterator_position.get(1) * iterator_position.get(1);
                r = sqrt(r);
                double uth = a_tc / r + b_tc * r;

                double ux = uth * (-iterator_position.get(1) / r);
                double uy = uth * (iterator_position.get(0) / r);

                vd.template getLastProp<velocity>()[0] = ux;
                vd.template getLastProp<velocity>()[1] = uy;

                for (int xyz = 0; xyz < DIM; xyz++)
                {
                    vd.template getLastProp<force_transport>()[xyz] = 0.0;
                }
            }
        }
        else // skip fluid particle
        {
            ++fluid_it;
            continue;
        }
        // Set properties
        vd.template getLastProp<rho>() = params.rho_zero;
        vd.template getLastProp<pressure>() = 0.0;
        vd.template getLastProp<drho>() = 0.0;

        for (int xyz = 0; xyz < DIM; xyz++)
        {
            vd.getLastPos()[xyz] = iterator_position.get(xyz);
            vd.template getLastProp<v_transport>()[xyz] = 0.0;
            vd.template getLastProp<normal_vector>()[xyz] = 0.0;
        }

        vd.template getLastProp<curvature_boundary>() = 0.0;
        vd.template getLastProp<arc_length>() = dp;

        // next fluid particle
        ++fluid_it;
    }
}
void CreateParticleGeometryStep(particles &vd, std::vector<std::pair<probe_particles, int>> &vp_vec, Vcluster<> &v_cl, Parameters params)
{
    // Size of the virtual cartesian grid that defines where to place the particles
    size_t sz[DIM];

    double length_small[DIM];
    double length_big[DIM];
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
    double offset_domain[DIM] = {0.0};
    double offset_recipient[DIM] = {0.0};
    double offset_periodic_fluid[DIM] = {0.0};
    double offset_periodic_recipient[DIM] = {0.0};

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

    double StepHeight = 4.9;
    params.LengthScale = StepHeight;

    Nfluid_big[0] = 275;
    Nfluid_big[1] = 31;

    Nfluid_small[0] = 100;
    Nfluid_small[1] = 16;

    size_t Nbound = (params.BC_TYPE == NEW_NO_SLIP) ? 1 : 3;
    Nboundary_big[0] = 0;
    Nboundary_big[1] = Nbound;
    // double Nboundary_small_up = Nbound;
    // double Nboundary_small_down = Nbound + Nfluid_big[1] - Nfluid_small[1];

    bc[0] = PERIODIC;
    bc[1] = NON_PERIODIC;
    params.dp = params.LengthScale / ((double)Nfluid_big[1] - Nfluid_small[1]);
    double dp = params.dp;
    params.umax = 1.4 * 1e-1;

    params.H = params.Hconst * dp;
    // r_threshold = (KERNEL == CUBIC ? 2.0 * H : 3.0 * H);
    params.r_threshold = 3.0 * params.H;
    params.Kquintic = (DIM == 3) ? 1.0 / 120.0 / M_PI / params.H / params.H / params.H : 7.0 / 478.0 / M_PI / params.H / params.H;
    params.MassFluid = params.rho_zero * (DIM == 3 ? dp * dp * dp : dp * dp);
    params.MassBound = params.rho_zero * (DIM == 3 ? dp * dp * dp : dp * dp);
    params.cbar = params.coeff_sound * params.umax;
    params.B = params.rho_zero * params.cbar * params.cbar / params.gamma_;
    params.Pbackground = params.Bfactor * params.B;
    params.eta = params.nu * params.rho_zero;
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
    Box<DIM, double> domain({-offset_domain[0],
                             -offset_domain[1]},
                            {params.length[0] + offset_domain[0],
                             params.length[1] + offset_domain[1]});

    Box<DIM, double> fluid_box_small({0.0,
                                      (Nfluid_big[1] - Nfluid_small[1]) * dp},
                                     {length_small[0],
                                      (Nfluid_big[1] - Nfluid_small[1]) * dp + length_small[1]});

    Box<DIM, double> fluid_box_big({length_small[0],
                                    0.0},
                                   {length_small[0] + length_big[0] + offset_periodic_fluid[0],
                                    length_big[1]});

    Box<DIM, double> recipient({-offset_recipient[0],
                                -offset_recipient[1]},
                               {params.length[0] + offset_recipient[0] + offset_periodic_recipient[0],
                                params.length[1] + offset_recipient[1] + offset_periodic_recipient[1]});

    // Will only be used in the new bc
    Box<DIM, double> recipient_hole_small({offset_recipient[0],
                                           (Nfluid_big[1] - Nfluid_small[1]) * dp + offset_recipient[1]},
                                          {length_small[0] - offset_recipient[0],
                                           (Nfluid_big[1] - Nfluid_small[1]) * dp + length_small[1] - offset_recipient[1]});

    Box<DIM, double> recipient_hole_big({length_small[0] + offset_recipient[0],
                                         offset_recipient[1]},
                                        {length_small[0] + length_big[0] - offset_recipient[0] + offset_periodic_fluid[0],
                                         length_big[1] - offset_recipient[1]});

    Box<DIM, double> CornerHole{{3 * dp, -3 * dp}, {(3 + Nfluid_small[0] - 6) * dp, (Nfluid_big[1] - Nfluid_small[1] - 3) * dp}};
    Box<DIM, double> CornerHole_New{{dp, -1 * dp}, {(1 + Nfluid_small[0] - 2) * dp, (Nfluid_big[1] - Nfluid_small[1]) * dp - 0.5 * dp}};

    // extended boundary around the domain, and the processor domain
    Ghost<DIM, double> g(params.r_threshold);

    // create particle object
    particles vd_loc(0, domain, bc, g, DEC_GRAN(512));
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
    Point<DIM, double> P1 = {0.75 * Nfluid_small[0] * dp, 0.0};
    Point<DIM, double> P2 = {1.2 * Nfluid_small[0] * dp, 0.0};
    Point<DIM, double> P3 = {1.4 * Nfluid_small[0] * dp, 0.0};
    Point<DIM, double> P4 = {Nfluid_small[0] * dp + Nfluid_big[0] * dp * 0.7, 0.0};

    std::vector<Point<DIM, double>> ProbePoints = {P1, P2, P3, P4};

    Point<DIM, double> VerticalOffset = {0.0, dp};
    int k0 = 0;
    int kendHeight = Nfluid_big[1];

    // place probes
    if (params.PROBES_ENABLED)
    {
        for (int k = 0; k < 4; k++)
        {
            // create probe object
            Ghost<DIM, double> gp(0);
            size_t bc_p[DIM] = {NON_PERIODIC, NON_PERIODIC};
            probe_particles vp_loc(0, domain, bc_p, gp, DEC_GRAN(512));
            openfpm::vector<std::string> names_p({"vx"});
            vp_loc.setPropNames(names_p);

            PlaceProbes(vp_loc, k0, kendHeight, ProbePoints[k], VerticalOffset);
            std::pair<probe_particles, int> tmp = std::make_pair(vp_loc, 0);
            vp_vec.push_back(tmp);
            params.probe_filenames.push_back("probes_" + std::to_string(k) + "_" + params.filename);
        }
    }

    // return an iterator to the fluid particles to add to vd
    auto fluid_it1 = DrawParticles::DrawBox(vd, sz, domain, fluid_box_big);
    auto fluid_it2 = DrawParticles::DrawBox(vd, sz, domain, fluid_box_small);

    // for each particle inside the fluid box ...
    while (fluid_it1.isNext())
    {

        Point<DIM, double> iterator_position = fluid_it1.get();

        // ... add a particle ...
        vd.add();
        vd.template getLastProp<type>() = FLUID;

        // Set properties
        vd.template getLastProp<rho>() = params.rho_zero;
        vd.template getLastProp<pressure>() = 0.0;
        vd.template getLastProp<drho>() = 0.0;

        for (int xyz = 0; xyz < DIM; xyz++)
        {
            vd.getLastPos()[xyz] = iterator_position.get(xyz);
            vd.template getLastProp<force_transport>()[xyz] = 0.0;
            vd.template getLastProp<v_transport>()[xyz] = 0.0;
            vd.template getLastProp<normal_vector>()[xyz] = 0.0;
            vd.template getLastProp<velocity>()[xyz] = 0.0;
        }

        vd.template getLastProp<curvature_boundary>() = 0.0;
        vd.template getLastProp<arc_length>() = dp;

        // next fluid particle
        ++fluid_it1;
    }
    // for each particle inside the fluid box ...
    while (fluid_it2.isNext())
    {

        Point<DIM, double> iterator_position = fluid_it2.get();

        // ... add a particle ...
        vd.add();
        vd.template getLastProp<type>() = FLUID;

        // Set properties
        vd.template getLastProp<rho>() = params.rho_zero;
        vd.template getLastProp<pressure>() = 0.0;
        vd.template getLastProp<drho>() = 0.0;

        for (int xyz = 0; xyz < DIM; xyz++)
        {
            vd.getLastPos()[xyz] = iterator_position.get(xyz);
            vd.template getLastProp<force_transport>()[xyz] = 0.0;
            vd.template getLastProp<v_transport>()[xyz] = 0.0;
            vd.template getLastProp<normal_vector>()[xyz] = 0.0;
            vd.template getLastProp<velocity>()[xyz] = 0.0;
        }

        vd.template getLastProp<curvature_boundary>() = 0.0;
        vd.template getLastProp<arc_length>() = dp;

        // next fluid particle
        ++fluid_it2;
    }

    // Now place solid walls
    openfpm::vector<Box<DIM, double>> holes;

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
    Box<DIM, double> hole_get0 = holes.get(0);
    Box<DIM, double> hole_get1 = holes.get(1);
    auto bound_box = DrawParticles::DrawSkin(vd, sz, domain, holes, recipient);

    if (bc[0] != PERIODIC || bc[1] != PERIODIC) // no walls in all periodic scenario
    {
        while (bound_box.isNext())
        {
            Point<DIM, double> position = bound_box.get();

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
                // double remx = fmod(position.get(0), dp);
                double remy = fmod(position.get(1), dp);
                double tol = 0.5 * dp * 10e-2;

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

            vd.template getLastProp<type>() = BOUNDARY;
            vd.template getLastProp<rho>() = params.rho_zero;
            vd.template getLastProp<pressure>() = 0.0;
            vd.template getLastProp<drho>() = 0.0;

            for (int xyz = 0; xyz < DIM; xyz++)
            {
                vd.getLastPos()[xyz] = bound_box.get().get(xyz);
                vd.template getLastProp<force>()[xyz] = 0.0;
                vd.template getLastProp<force_transport>()[xyz] = 0.0;
                vd.template getLastProp<v_transport>()[xyz] = 0.0;
                vd.template getLastProp<normal_vector>()[xyz] = 0.0;
                if (position.get(1) < dp / 4.0) // bottom wall
                {
                    vd.template getLastProp<velocity>()[xyz] = params.vw_bottom.get(xyz);
                }
                else if (position.get(1) > params.length[1] - dp / 4.0) // top wall
                {
                    vd.template getLastProp<velocity>()[xyz] = params.vw_top.get(xyz);
                }
            }

            vd.template getLastProp<curvature_boundary>() = 0.0;
            vd.template getLastProp<arc_length>() = dp;

            ++bound_box;
        }
    }
}