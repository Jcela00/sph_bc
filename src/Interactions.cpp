#include "Interactions.hpp"

void interact_fluid_boundary_new(particles &vd,
                                 vect_dist_key_dx fluid_key,
                                 const double &massf,
                                 const double &rhof,
                                 const double &Pf,
                                 const Point<DIM, double> &xf,
                                 const Point<DIM, double> &vf,
                                 const Point<DIM, double> &xw,
                                 const Point<DIM, double> &r_wall_to_fluid,
                                 unsigned long &boundary_key,
                                 bool accumulate_force,
                                 double &obstacle_force_x,
                                 double &obstacle_force_y,
                                 const Parameters &params)
{

    double dp = params.dp;
    // Points from fluid to wall
    Point<DIM, double> r_fluid_to_wall = -1.0 * r_wall_to_fluid;
    // double dist2marker = getVectorNorm(r_fluid_to_wall);
    Point<DIM, double> vw = vd.getProp<velocity>(boundary_key);

    double ang_vel = vd.getProp<vd_omega>(boundary_key);

    if (ang_vel != 0.0) // if the solid is rotating, we need to add the tangential velocity of the rotation to vw
    {
        // marker particles store centre of solid body in force_transport since it is unused
        // vector pointing from centre of rotation to marker particle
        const Point<DIM, double> radial_vec = {xw.get(0) - vd.getProp<force_transport>(boundary_key)[0],
                                               xw.get(1) - vd.getProp<force_transport>(boundary_key)[1]};
        const double radius = getVectorNorm(radial_vec);

        // get vector tangential to the radial vector, rotation velocity is in this direction
        const Point<DIM, double> tangential_rotation = getPerpendicularUnit2D(radial_vec);

        // Wall velocity is linear velocity + w*R*tangential
        vw.get(0) += radius * ang_vel * tangential_rotation.get(0);
        vw.get(1) += radius * ang_vel * tangential_rotation.get(1);
    }

    // Get normal and tangential vectors for velocity mirroring
    const Point<DIM, double> normal = vd.getProp<normal_vector>(boundary_key);
    const Point<DIM, double> tangential = getPerpendicularUnit2D(normal);

    // wall acceleration
    const Point<DIM, double> aw = vd.getProp<force>(boundary_key);

    // Difference between fluid transport and momentum velocity
    const Point<DIM, double> vtf = vd.getProp<v_transport>(fluid_key);
    const Point<DIM, double> vdiff_f = vtf - vf;

    // get curvature, and arc length
    // double kappa = vd.getProp<curvature_boundary>(boundary_key);

    // Project vf and vw on tangential and normal directions
    double vt = dotProduct(vf, tangential);
    double vn = dotProduct(vf, normal);
    double vwt = dotProduct(vw, tangential);

    // vertical distance from fluid particle to wall
    double lf = dotProduct(r_fluid_to_wall, normal);
    lf = (lf < 0.0 ? -1.0 * lf : lf); // absolute value

    // Get array of vectors from fluid to 3 boundary particles and its norms
    std::array<Point<DIM, double>, 3> r_boundary = getBoundaryPositions(r_fluid_to_wall, normal, dp);
    std::array<double, 3> r_boundary_norm = {getVectorNorm(r_boundary[0]), getVectorNorm(r_boundary[1]), getVectorNorm(r_boundary[2])};

    // const double dist2third = r_boundary_norm[2];
    // distance from 3 boundary particles to marker
    std::array<double, 3> lwall = {0.5 * dp, 1.5 * dp, 2.5 * dp};

    // project forces on normal direction
    double g_normal = dotProduct(params.gravity_vector, normal);
    double a_normal = dotProduct(aw, normal);

    // For gradient of Af tensor
    const std::array<Point<DIM, double>, DIM> Af = dyadicProduct(rhof * vf, vdiff_f);
    // to avoid division by zero
    lf = std::max(lf, 0.25 * dp);

    const Point<3, double> vol = vd.getProp<vd_volume>(boundary_key);

    for (int i = 0; i < 3; i++) // for the 3 boundary particles
    {
        // double rmax = sqrt(3.0 * 3.0 - (0.5 + (double)i) * (0.5 + (double)i)) * dp;
        // double rmin = (3.0 - (0.5 + (double)i)) * dp;
        // double kappa_max = 1.0 / (3.0 * dp);

        // kappa 0 gets rmax, kappa = kappa_max gets rmin
        // double r_interp = (rmin - rmax) / kappa_max * kappa + rmax;
        // double r_interp = rmin;
        // if (dist2marker < r_interp)
        // {
        const double Mass_boundary = vol[i] * params.rho_zero;
        const Point<DIM, double> v_boundary = ((vwt - vt) * (lwall[i] / lf) + vwt) * tangential + vn * normal;
        const double p_boundary = Pf - rhof * (g_normal - a_normal) * (lf + lwall[i]); //  dot(r_boundary,normal) = -(lf+lw)
        const double rho_boundary = InvEqState_particle(p_boundary, params.rho_zero, params.B, params.gamma_, params.xi);
        // if (i == 0)
        // {
        //     // set the marker particle properties equal to the first boundary particle, for visualization purposes only, they are not used in the simulation
        //     vd.getProp<pressure>(boundary_key) = p_boundary;
        //     vd.getProp<rho>(boundary_key) = rho_boundary;
        // }

        // flip sign of r_boundary to get vector pointing from boundary to fluid (Force routines use the vector pointing from b to a)
        r_boundary[i] = -1.0 * r_boundary[i];

        // Evaluate kernel gradient
        const Point<DIM, double> DW = DWab(r_boundary[i], r_boundary_norm[i], params.H, params.Kquintic);

        // Compute forces
        const Point<DIM, double> v_rel = vf - v_boundary;
        // const double Va2 = (massf / rhof) * (massf / rhof);
        // const double Vb2 = Volumes[i] * Volumes[i];
        // const double Vb2 = (Mass_boundary / rho_boundary) * (Mass_boundary / rho_boundary); // this is Vol_boundary^2*(rhob/rho0)^2 therefore considers density changes,
        const double Vb = (Mass_boundary / rho_boundary);
        // But I found no big difference with the commented out expression

        const Point<DIM, double> ViscosityTerm = Pi_physical(r_boundary[i], r_boundary_norm[i], v_rel, DW, params.eta);
        const double PressureTerm = PressureForce(rhof, rho_boundary, Pf, p_boundary); //-p_boundary - Pf;
        const Point<DIM, double> GradATerm = 0.5 * matVec(Af, DW);

        // std::cout << "PressureTerm_x: " << PressureTerm * DW.get(0) << "PressureTerm_y: " << PressureTerm * DW.get(1) << std::endl;
        // std::cout << "ViscosityTerm_x: " << ViscosityTerm.get(0) << "ViscosityTerm_y: " << ViscosityTerm.get(1) << std::endl;
        // std::cout << "GradATerm_x: " << GradATerm.get(0) << "GradATerm_y: " << GradATerm.get(1) << std::endl;
        if (accumulate_force) // we accumulate x force on cylinder ( this is just to compute drag coefficient)
        {
            obstacle_force_x += -2.0 * (Vb / rhof) * (PressureTerm * DW.get(0) + ViscosityTerm.get(0) + GradATerm.get(0));
            obstacle_force_y += -2.0 * (Vb / rhof) * (PressureTerm * DW.get(1) + ViscosityTerm.get(1) + GradATerm.get(1));
        }
        for (int xyz = 0; xyz < DIM; ++xyz)
        {
            // write to particles
            vd.getProp<force>(fluid_key)[xyz] += 2.0 * (Vb / rhof) * (PressureTerm * DW.get(xyz) + ViscosityTerm.get(xyz) + GradATerm.get(xyz));
            vd.getProp<force_transport>(fluid_key)[xyz] += -2.0 * (Vb / rhof) * (params.Pbackground) * DW.get(xyz);
            // vd.getProp<force>(fluid_key)[xyz] += (Va2 + Vb2) * (PressureTerm * DW.get(xyz) + ViscosityTerm.get(xyz) + GradATerm.get(xyz)) / massf;
            // vd.getProp<force_transport>(fluid_key)[xyz] += -1.0 * (Va2 + Vb2) * (Pbackground)*DW.get(xyz) / massf;
        }

        if (params.DENSITY_TYPE == DENSITY_DIFFERENTIAL) // this doesnt work well I havent touched in a long time and I have made may changes
        {
            vd.getProp<drho>(fluid_key) += Mass_boundary * dotProduct(v_rel, DW);
        }
        // }
    }
}

void interact_fluid_boundary_old(particles &vd,
                                 vect_dist_key_dx &fluid_key,
                                 const double &massf,
                                 const double &rhof,
                                 const double &Pf,
                                 const Point<DIM, double> &xf,
                                 const Point<DIM, double> &vf,
                                 const Point<DIM, double> &xb,
                                 const Point<DIM, double> &r_ab,
                                 const double &r2,
                                 unsigned long &boundary_key,
                                 bool accumulate_force,
                                 double &obstacle_force_x,
                                 double &obstacle_force_y,
                                 const Parameters &params)
{
    const double massb = params.MassBound;
    const double rhob = vd.getProp<rho>(boundary_key);
    const double Pb = vd.getProp<pressure>(boundary_key);
    const Point<DIM, double> vb = vd.getProp<velocity>(boundary_key);
    const Point<DIM, double> vb_noslip = vd.getProp<v_transport>(boundary_key); // here we store the extrapolated velocity for no slip BC

    const double r = sqrt(r2);

    const Point<DIM, double> v_rel = vf - vb;
    const Point<DIM, double> v_rel_aux = vf - vb_noslip;

    const Point<DIM, double> DW = DWab(r_ab, r, params.H, params.Kquintic);

    const double Va2 = (massf / rhof) * (massf / rhof);
    const double Vb2 = (massb / rhob) * (massb / rhob);

    const Point<DIM, double> ViscosityTerm = Pi_physical(r_ab, r, v_rel_aux, DW, params.eta);
    const double PressureTerm = PressureForce(rhof, rhob, Pf, Pb);

    const Point<DIM, double> vtf = vd.getProp<v_transport>(fluid_key);
    const Point<DIM, double> vdiff_f = vtf - vf;
    // Boundary particles have no transport velocity difference
    std::array<Point<DIM, double>, DIM> Af = dyadicProduct(rhof * vf, vdiff_f);
    Point<DIM, double> GradATerm = 0.5 * matVec(Af, DW);

    if (accumulate_force) // we accumulate x force on cylinder
    {
        obstacle_force_x += -1.0 * (Va2 + Vb2) * (PressureTerm * DW.get(0) + ViscosityTerm.get(0) + GradATerm.get(0)) / massf;
        obstacle_force_y += -1.0 * (Va2 + Vb2) * (PressureTerm * DW.get(1) + ViscosityTerm.get(1) + GradATerm.get(1)) / massf;
    }
    for (int xyz = 0; xyz < DIM; ++xyz)
    {
        vd.getProp<force>(fluid_key)[xyz] += (Va2 + Vb2) * (PressureTerm * DW.get(xyz) + ViscosityTerm.get(xyz) + GradATerm.get(xyz)) / massf;
        vd.getProp<force_transport>(fluid_key)[xyz] += -1.0 * (Va2 + Vb2) * params.Pbackground * DW.get(xyz) / massf;
    }

    if (params.DENSITY_TYPE == DENSITY_DIFFERENTIAL)
    {
        vd.getProp<drho>(fluid_key) += massb * dotProduct(v_rel, DW);
    }
}

void interact_fluid_fluid(particles &vd,
                          const vect_dist_key_dx &a_key,
                          const double &massa,
                          const double &rhoa,
                          const double &Pa,
                          const Point<DIM, double> &xa,
                          const Point<DIM, double> &va,
                          const Point<DIM, double> &xb,
                          const Point<DIM, double> &r_ab,
                          const double &r2,
                          const unsigned long &b_key,
                          const Parameters &params)
{

    const double massb = params.MassFluid;
    const double rhob = vd.getProp<rho>(b_key);
    const double Pb = vd.getProp<pressure>(b_key);
    const Point<DIM, double> vb = vd.getProp<velocity>(b_key);

    const double r = sqrt(r2);

    const Point<DIM, double> v_rel = va - vb;

    const Point<DIM, double> DW = DWab(r_ab, r, params.H, params.Kquintic);

    const double Va2 = (massa / rhoa) * (massa / rhoa);
    const double Vb2 = (massb / rhob) * (massb / rhob);

    const Point<DIM, double> ViscosityTerm = Pi_physical(r_ab, r, v_rel, DW, params.eta);
    const double PressureTerm = PressureForce(rhoa, rhob, Pa, Pb);

    const Point<DIM, double> vta = vd.getProp<v_transport>(a_key);
    const Point<DIM, double> vtb = vd.getProp<v_transport>(b_key);
    const Point<DIM, double> vdiff_a = vta - va;
    const Point<DIM, double> vdiff_b = vtb - vb;

    // double GradATerm = 0.5 * (rhoa * dotProduct(va, vdiff_a) + rhob * dotProduct(vb, vdiff_b));
    std::array<Point<DIM, double>, DIM> Aa = dyadicProduct(rhoa * va, vdiff_a);
    std::array<Point<DIM, double>, DIM> Ab = dyadicProduct(rhob * vb, vdiff_b);
    std::array<Point<DIM, double>, DIM> SumA;
    SumA[0] = Aa[0] + Ab[0];
    SumA[1] = Aa[1] + Ab[1];
    if (DIM == 3)
        SumA[2] = Aa[2] + Ab[2];

    Point<DIM, double> GradATerm = 0.5 * matVec(SumA, DW);

    for (int xyz = 0; xyz < DIM; ++xyz)
    {
        vd.getProp<force>(a_key)[xyz] += (Va2 + Vb2) * (PressureTerm * DW.get(xyz) + ViscosityTerm.get(xyz) + GradATerm.get(xyz)) / massa;
        vd.getProp<force_transport>(a_key)[xyz] += -1.0 * (Va2 + Vb2) * params.Pbackground * DW.get(xyz) / massa;
    }
    if (params.DENSITY_TYPE == DENSITY_DIFFERENTIAL)
    {
        vd.getProp<drho>(a_key) += massb * dotProduct(v_rel, DW);
    }
}
