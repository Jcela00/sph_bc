#include "Obstacle.hpp"

void AddFlatWallNewBC(particles &vd,
                      const int k0,
                      const int kmax,
                      const Point<DIM, real_number> Corner,
                      const Point<DIM, real_number> UnitOffset,
                      const real_number dx,
                      const Point<DIM, real_number> obstacle_centre,
                      const Point<DIM, real_number> obstacle_velocity,
                      const Parameters &arg_p,
                      const size_t particle_type,
                      const real_number obstacle_omega)
{

    for (int k = k0; k < kmax; k++)
    {
        Point<DIM, real_number> WallPos = Corner + k * UnitOffset;

        vd.add();
        vd.template getLastProp<vd0_type>() = particle_type;
        for (int xyz = 0; xyz < DIM; xyz++)
        {
            vd.getLastPos()[xyz] = WallPos.get(xyz); //+ ((real_number)rand() / RAND_MAX - 0.5) * dp;
            vd.template getLastProp<vd4_velocity>()[xyz] = obstacle_velocity.get(xyz);
            vd.template getLastProp<vd6_force>()[xyz] = 0.0;
            vd.template getLastProp<vd7_force_t>()[xyz] = obstacle_centre.get(xyz);
            vd.template getLastProp<vd5_velocity_t>()[xyz] = 0.0;
            vd.template getLastProp<vd8_normal>()[xyz] = 0.0;
        }

        vd.template getLastProp<vd2_pressure>() = 0.0;
        vd.template getLastProp<vd1_rho>() = arg_p.rho0;
        vd.template getLastProp<vd3_drho>() = 0.0;
        vd.template getLastProp<vd9_volume>()[0] = dx;
        vd.template getLastProp<vd10_omega>() = obstacle_omega;
    }
}

void AddFlatWallModNewBC(particles &vd,
                         const int k0,
                         const int kmax,
                         const Point<DIM, real_number> Corner,
                         const Point<DIM, real_number> UnitOffset,
                         const real_number dx,
                         const Point<DIM, real_number> obstacle_centre,
                         const Point<DIM, real_number> obstacle_velocity,
                         const Parameters &arg_p,
                         const size_t particle_type,
                         const Point<DIM, real_number> given_normal,
                         const real_number obstacle_omega)
{

    for (int k = k0; k < kmax; k++)
    {
        Point<DIM, real_number> WallPos = Corner + k * UnitOffset;

        vd.add();
        vd.template getLastProp<vd0_type>() = particle_type;
        for (int xyz = 0; xyz < DIM; xyz++)
        {
            vd.getLastPos()[xyz] = WallPos.get(xyz); //+ ((real_number)rand() / RAND_MAX - 0.5) * dp;
            vd.template getLastProp<vd4_velocity>()[xyz] = obstacle_velocity.get(xyz);
            vd.template getLastProp<vd6_force>()[xyz] = 0.0;
            vd.template getLastProp<vd7_force_t>()[xyz] = obstacle_centre.get(xyz);
            vd.template getLastProp<vd5_velocity_t>()[xyz] = 0.0;
            vd.template getLastProp<vd8_normal>()[xyz] = given_normal.get(xyz);
        }

        vd.template getLastProp<vd2_pressure>() = 0.0;
        vd.template getLastProp<vd1_rho>() = arg_p.rho0;
        vd.template getLastProp<vd3_drho>() = 0.0;
        vd.template getLastProp<vd9_volume>()[0] = dx;
        vd.template getLastProp<vd10_omega>() = obstacle_omega;
    }
}

bool isAvobeLine(Point<DIM, real_number> P, Point<DIM, real_number> Q, Point<DIM, real_number> EvalPoint, real_number dp)
{
    // P,Q two points forming a line, EvalPoint the point to check if it is above the line or not
    // with a small epsilon to avoid particles to be very close to the line
    real_number slope = (Q.get(1) - P.get(1)) / (Q.get(0) - P.get(0));
    // we compare y- y1 + m(x-x1)
    // =0 means the point is in the line
    // >0 means the point is above the line
    // <0 means the point is below the line
    real_number yline = P.get(1) + slope * (EvalPoint.get(0) - P.get(0));
    real_number epsilon = 0.10 * dp;
    if (EvalPoint.get(1) > yline + epsilon)
        return true;
    else
        return false;
}

bool isBelowLine(Point<DIM, real_number> P, Point<DIM, real_number> Q, Point<DIM, real_number> EvalPoint, real_number dp)
{
    // P,Q two points forming a line, EvalPoint the point to check if it is above the line or not
    // with a small epsilon to avoid particles to be very close to the line
    real_number slope = (Q.get(1) - P.get(1)) / (Q.get(0) - P.get(0));
    // we compare y- y1 + m(x-x1)
    // =0 means the point is in the line
    // >0 means the point is above the line
    // <0 means the point is below the line
    real_number yline = P.get(1) + slope * (EvalPoint.get(0) - P.get(0));
    real_number epsilon = 0.10 * dp;
    if (EvalPoint.get(1) < yline - epsilon)
        return true;
    else
        return false;
}

// Obstacle class
Obstacle::Obstacle(Point<DIM, real_number> Centre,
                   const Parameters &p,
                   Point<DIM, real_number> vel,
                   real_number omega,
                   real_number rf) : Centre_(Centre), params_(p), LinearVelocity_(vel), AngularVelocity_(omega), refine_factor(rf) {}
Obstacle::Obstacle(const Parameters &p) : params_(p) {}

// EmptyObstacle class
EmptyObstacle::EmptyObstacle(const Parameters &p) : Obstacle(p) {}
bool EmptyObstacle::isInside(Point<DIM, real_number> P)
{
    return false;
}
void EmptyObstacle::AddObstacle(particles &vd)
{
}

// CylinderObstacle class
CylinderObstacle::CylinderObstacle(real_number Radius,
                                   Point<DIM, real_number> centre,
                                   const Parameters &p,
                                   Point<DIM, real_number> vel,
                                   real_number omega,
                                   real_number rf) : Obstacle(centre, p, vel, omega, rf), Radius_(Radius) {};

bool CylinderObstacle::isInside(Point<DIM, real_number> P)
{
    bool is_inside = false;
    real_number tmp = (P.get(0) - Centre_.get(0)) * (P.get(0) - Centre_.get(0)) + (P.get(1) - Centre_.get(1)) * (P.get(1) - Centre_.get(1));
    if (tmp <= Radius_ * Radius_)
    {
        is_inside = true;
    }
    return is_inside;
}
bool CylinderObstacle::isInside_minEps(Point<DIM, real_number> P) // for outer cylinder in taylor couette
{
    bool is_outside = false;
    real_number tmp = (P.get(0) - Centre_.get(0)) * (P.get(0) - Centre_.get(0)) + (P.get(1) - Centre_.get(1)) * (P.get(1) - Centre_.get(1));
    if (tmp <= Radius_ * Radius_)
    {
        is_outside = true;
    }
    return is_outside;
}

void CylinderObstacle::AddObstacle(particles &vd)
{
    const real_number dx_self = params_.dp / refine_factor;
    const real_number perimeter = 2.0 * M_PI * Radius_;
    const int Np_cylinder = ceil(perimeter / dx_self);
    const real_number dtheta = 2.0 * M_PI / (real_number)Np_cylinder;
    const real_number dxwall = dtheta * Radius_;
    real_number theta = 0.0;

    Point<DIM, real_number> Cylinder_particle;

    for (int k = 0; k < Np_cylinder; k++)
    {

        Cylinder_particle[0] = Centre_.get(0) + Radius_ * cos(theta);
        Cylinder_particle[1] = Centre_.get(1) + Radius_ * sin(theta);

        if (DIM == 3)
        {
            Cylinder_particle[2] = Centre_.get(2);
        }

        vd.add();
        vd.template getLastProp<vd0_type>() = OBSTACLE;
        for (int xyz = 0; xyz < DIM; xyz++)
        {
            vd.getLastPos()[xyz] = Cylinder_particle.get(xyz); // + ((real_number)rand() / RAND_MAX - 0.5) * dp;
            vd.template getLastProp<vd6_force>()[xyz] = 0.0;
            vd.template getLastProp<vd7_force_t>()[xyz] = Centre_.get(xyz);
            vd.template getLastProp<vd5_velocity_t>()[xyz] = 0.0;
            vd.template getLastProp<vd8_normal>()[xyz] = 0.0;
            vd.template getLastProp<vd4_velocity>()[xyz] = LinearVelocity_.get(xyz);
        }

        vd.template getLastProp<vd2_pressure>() = 0.0;
        vd.template getLastProp<vd1_rho>() = params_.rho0;
        vd.template getLastProp<vd3_drho>() = 0.0;
        vd.template getLastProp<vd9_volume>()[0] = dxwall;
        vd.template getLastProp<vd10_omega>() = AngularVelocity_;
        theta += dtheta;
    }
}

// Ellipse class
EllipticObstacle::EllipticObstacle(real_number Major,
                                   real_number Minor,
                                   real_number tilt,
                                   Point<DIM, real_number> centre,
                                   const Parameters &p,
                                   Point<DIM, real_number> vel,
                                   real_number omega,
                                   real_number rf) : Obstacle(centre, p, vel, omega, rf), Major_(Major), Minor_(Minor), tilt_(-tilt * M_PI / 180.0) {}

bool EllipticObstacle::isInside(Point<DIM, real_number> P)
{
    // get point minus the centre
    Point<DIM, real_number> Paux = P - Centre_;
    // rotate by -tilt
    Paux = ApplyRotation(Paux, -tilt_, {0.0, 0.0});
    // check if the point is inside the ellipse
    return (Paux.get(0) * Paux.get(0) / (Major_ * Major_) + Paux.get(1) * Paux.get(1) / (Minor_ * Minor_) <= 1.0);
}

void EllipticObstacle::AddObstacle(particles &vd)
{
    // first we need to determine the perimeter of the ellipse
    // we can use c++ std::comp_ellint_2 function
    // first compute eccentricity
    const real_number e = sqrtf(1.0 - Minor_ * Minor_ / Major_ / Major_);
    // compute the perimeter
    const real_number perimeter = 4.0 * Major_ * std::comp_ellint_2(e);
    const real_number dx_self = params_.dp / refine_factor;
    const int Np_ellipse_ = ceil(perimeter / dx_self);
    const real_number dxwall = perimeter / (real_number)Np_ellipse_;

    const int M = 1e6; // number of points to discretize the ellipse
    real_number accumulated_arc = 0.0;
    real_number theta = 0.0;
    real_number dtheta = 2.0 * M_PI / (real_number)M;

    Point<DIM, real_number> P = parametricEllipse(theta);

    // add the first particle
    vd.add();
    vd.template getLastProp<vd0_type>() = OBSTACLE;
    for (int xyz = 0; xyz < DIM; xyz++)
    {
        vd.getLastPos()[xyz] = P.get(xyz);
        vd.template getLastProp<vd6_force>()[xyz] = 0.0;
        vd.template getLastProp<vd7_force_t>()[xyz] = Centre_.get(xyz);
        vd.template getLastProp<vd5_velocity_t>()[xyz] = 0.0;
        vd.template getLastProp<vd8_normal>()[xyz] = 0.0;
        vd.template getLastProp<vd4_velocity>()[xyz] = LinearVelocity_.get(xyz);
    }

    vd.template getLastProp<vd2_pressure>() = 0.0;
    vd.template getLastProp<vd1_rho>() = params_.rho0;
    vd.template getLastProp<vd3_drho>() = 0.0;
    vd.template getLastProp<vd9_volume>()[0] = dxwall;
    vd.template getLastProp<vd10_omega>() = AngularVelocity_;

    Point<DIM, real_number> P_prev = P;

    theta += dtheta;

    for (int step = 1; step < M; step++)
    {
        P = parametricEllipse(theta);
        accumulated_arc += norm(P - P_prev);

        P_prev = P;
        theta += dtheta;

        // if accumulated arc is greater than dxwall, we add a particle
        if (accumulated_arc > dxwall)
        {
            vd.add();
            vd.template getLastProp<vd0_type>() = OBSTACLE;
            for (int xyz = 0; xyz < DIM; xyz++)
            {
                vd.getLastPos()[xyz] = P.get(xyz);
                vd.template getLastProp<vd6_force>()[xyz] = 0.0;
                vd.template getLastProp<vd7_force_t>()[xyz] = Centre_.get(xyz);
                vd.template getLastProp<vd5_velocity_t>()[xyz] = 0.0;
                vd.template getLastProp<vd8_normal>()[xyz] = 0.0;
                vd.template getLastProp<vd4_velocity>()[xyz] = LinearVelocity_.get(xyz);
            }

            vd.template getLastProp<vd2_pressure>() = 0.0;
            vd.template getLastProp<vd1_rho>() = params_.rho0;
            vd.template getLastProp<vd3_drho>() = 0.0;
            vd.template getLastProp<vd9_volume>()[0] = dxwall;
            vd.template getLastProp<vd10_omega>() = AngularVelocity_;

            // reset accumulated arc
            accumulated_arc = 0.0;
        }
    }
}

Point<DIM, real_number> EllipticObstacle::parametricEllipse(real_number theta)
{
    // given the angle theta, return the point in the ellipse
    Point<DIM, real_number> P;
    P.get(0) = Centre_.get(0) + (Major_ * cos(theta) * cos(tilt_) - Minor_ * sin(theta) * sin(tilt_));
    P.get(1) = Centre_.get(1) + (Major_ * cos(theta) * sin(tilt_) + Minor_ * sin(theta) * cos(tilt_));
    return P;
}

// RectangleObstacle class
RectangleObstacle::RectangleObstacle(Point<DIM, real_number> centre,
                                     const Parameters &p,
                                     real_number BaseLength,
                                     real_number HeigthLength,
                                     Point<DIM, real_number> vel,
                                     real_number omega,
                                     real_number rf) : Obstacle(centre, p, vel, omega, rf), BaseLength_(BaseLength), HeigthLength_(HeigthLength)
{
    LowerLeft_ = Point<DIM, real_number>{Centre_.get(0) - BaseLength_ / 2.0f,
                                         Centre_.get(1) - HeigthLength_ / 2.0f};

    UpperRight_ = Point<DIM, real_number>{Centre_.get(0) + BaseLength_ / 2.0f,
                                          Centre_.get(1) + HeigthLength_ / 2.0f};

    // LowerLeft_ = Point<DIM, real_number>{Centre_.get(0) - (((real_number)BaseLength_ - 1.0) * params_.dp) / 2.0,
    //                                 Centre_.get(1) - (((real_number)HeigthLength_ - 1.0) * params_.dp) / 2.0};
    // LowerRight_ = Point<DIM, real_number>{Centre_.get(0) + (((real_number)BaseLength_ - 1.0) * params_.dp) / 2.0,
    //                                  Centre_.get(1) - (((real_number)HeigthLength_ - 1.0) * params_.dp) / 2.0};
    // UpperLeft_ = Point<DIM, real_number>{Centre_.get(0) - (((real_number)BaseLength_ - 1.0) * params_.dp) / 2.0,
    //                                 Centre_.get(1) + (((real_number)HeigthLength_ - 1.0) * params_.dp) / 2.0};
    // UpperRight_ = Point<DIM, real_number>{Centre_.get(0) + (((real_number)BaseLength_ - 1.0) * params_.dp) / 2.0,
    //                                  Centre_.get(1) + (((real_number)HeigthLength_ - 1.0) * params_.dp) / 2.0};

    Rectangle_ = Box<DIM, real_number>(LowerLeft_, UpperRight_);
}

bool RectangleObstacle::isInside(Point<DIM, real_number> P)
{
    return Rectangle_.isInside(P);
}

void RectangleObstacle::AddObstacle(particles &vd)
{

    const real_number dx_self = params_.dp / refine_factor;

    // Horizontal walls
    const int N_bottom = ceil(BaseLength_ / dx_self);
    const real_number dxwall_bottom = BaseLength_ / N_bottom;
    Point<DIM, real_number> Xoffset = {dxwall_bottom, 0.0};

    // Lower wall
    AddFlatWallNewBC(vd, 0, N_bottom + 1, LowerLeft_, Xoffset, dxwall_bottom, Centre_, LinearVelocity_, params_, OBSTACLE, AngularVelocity_);
    // Upper wall
    AddFlatWallNewBC(vd, 0, N_bottom + 1, UpperRight_, -1.0 * Xoffset, dxwall_bottom, Centre_, LinearVelocity_, params_, OBSTACLE, AngularVelocity_);

    // Vertical walls
    const int N_right = ceil(HeigthLength_ / dx_self);
    const real_number dxwall_right = HeigthLength_ / N_right;
    Point<DIM, real_number> Yoffset = {0.0, dxwall_right};

    // Left wall
    AddFlatWallNewBC(vd, 1, N_right, LowerLeft_, Yoffset, dxwall_right, Centre_, LinearVelocity_, params_, OBSTACLE, AngularVelocity_);
    // Right wall
    AddFlatWallNewBC(vd, 1, N_right, UpperRight_, -1.0 * Yoffset, dxwall_right, Centre_, LinearVelocity_, params_, OBSTACLE, AngularVelocity_);

    // WHEN using non mod flat wall this are the k and kmax values
    // // Lower wall
    // AddFlatWallModNewBC(vd, 0, N_bottom + 1, LowerLeft_, Xoffset, dxwall_bottom, Centre_, LinearVelocity_, {0.0, -1.0}, AngularVelocity_);
    // // Upper wall
    // AddFlatWallModNewBC(vd, 0, N_bottom + 1, UpperRight_, -1.0 * Xoffset, dxwall_bottom, Centre_, LinearVelocity_, {0.0, 1.0}, AngularVelocity_);

    // // Vertical walls
    // const int N_right = ceil(HeigthLength_ / dx_self);
    // const real_number dxwall_right = HeigthLength_ / N_right;
    // Point<DIM, real_number> Yoffset = {0.0, dxwall_right};

    // // Left wall
    // AddFlatWallModNewBC(vd, 1, N_right, LowerLeft_, Yoffset, dxwall_right, Centre_, LinearVelocity_, {-1.0, 0.0}, AngularVelocity_);
    // // Right wall
    // AddFlatWallModNewBC(vd, 1, N_right, UpperRight_, -1.0 * Yoffset, dxwall_right, Centre_, LinearVelocity_, {1.0, 0.0}, AngularVelocity_);
}

// TriangleObstacle class
TriangleObstacle::TriangleObstacle(Point<DIM, real_number> centre,
                                   const Parameters &p,
                                   real_number BaseLength,
                                   real_number HeigthLength,
                                   Point<DIM, real_number> vel,
                                   real_number omega,
                                   real_number rf) : Obstacle(centre, p, vel, omega, rf), BaseLength_(BaseLength), HeigthLength_(HeigthLength)
{
    LowerLeft_ = Point<DIM, real_number>{Centre_.get(0) - 2.0f * BaseLength_ / 3.0f,
                                         Centre_.get(1) - HeigthLength_ / 3.0f};

    LowerRight_ = Point<DIM, real_number>{Centre_.get(0) + BaseLength_ / 3.0f,
                                          Centre_.get(1) - HeigthLength_ / 3.0f};

    UpperRight_ = Point<DIM, real_number>{Centre_.get(0) + BaseLength_ / 3.0f,
                                          Centre_.get(1) + 2.0f * HeigthLength_ / 3.0f};

    ContainingRectangle_ = Box<DIM, real_number>(LowerLeft_, UpperRight_);
}

bool TriangleObstacle::isInside(Point<DIM, real_number> P)
{
    return (ContainingRectangle_.isInside(P) && !isAvobeLine(LowerLeft_, UpperRight_, P, params_.dp));
}

void TriangleObstacle::AddObstacle(particles &vd)
{
    const real_number dx_self = params_.dp / refine_factor;
    // Lower wall
    const int N_bottom = ceil(BaseLength_ / dx_self);
    const real_number dxwall_bottom = BaseLength_ / N_bottom;
    const Point<DIM, real_number> Xoffset = {dxwall_bottom, 0.0};
    AddFlatWallNewBC(vd, 0, N_bottom + 1, LowerLeft_, Xoffset, dxwall_bottom, Centre_, LinearVelocity_, params_, OBSTACLE, AngularVelocity_);

    // Right wall
    const int N_right = ceil(HeigthLength_ / dx_self);
    const real_number dxwall_right = HeigthLength_ / N_right;
    const Point<DIM, real_number> Yoffset = {0.0, dxwall_right};
    AddFlatWallNewBC(vd, 1, N_right + 1, LowerRight_, Yoffset, dxwall_right, Centre_, LinearVelocity_, params_, OBSTACLE, AngularVelocity_);

    //  Hypothenuse wall
    // We want particles spaced roughly by dp
    const real_number HypothenuseLength = sqrtf(BaseLength_ * BaseLength_ + HeigthLength_ * HeigthLength_);
    const int Ndiag = ceil(HypothenuseLength / dx_self);       // integer number of particles that can fit in the diagonal
    const real_number dxwall_diag = HypothenuseLength / Ndiag; // actual spacing between particles ( close to dp but not exactly)
    const real_number sin_theta = HeigthLength_ / HypothenuseLength;
    const real_number cos_theta = BaseLength_ / HypothenuseLength;
    const Point<DIM, real_number> Diagoffset{dxwall_diag * cos_theta, dxwall_diag * sin_theta};

    AddFlatWallNewBC(vd, 1, Ndiag, UpperRight_, -1.0 * Diagoffset, dxwall_diag, Centre_, LinearVelocity_, params_, OBSTACLE, AngularVelocity_);
}

TriangleEqui::TriangleEqui(Point<DIM, real_number> centre,
                           const Parameters &p,
                           real_number sidelength,
                           Point<DIM, real_number> vel,
                           real_number omega,
                           real_number rf) : Obstacle(centre, p, vel, omega, rf), SideLength_(sidelength)
{

    UpperRight_ = Point<DIM, real_number>{Centre_.get(0) + sqrtf(3.0f) * SideLength_ / 6.0f,
                                          Centre_.get(1) + SideLength_ / 2.0f};
    LowerRight_ = Point<DIM, real_number>{Centre_.get(0) + sqrtf(3.0f) * SideLength_ / 6.0f,
                                          Centre_.get(1) - SideLength_ / 2.0f};
    TriangleTip_ = Point<DIM, real_number>{Centre_.get(0) - sqrtf(3.0f) * SideLength_ / 3.0f,
                                           Centre_.get(1)};

    LowerLeft_ = Point<DIM, real_number>{TriangleTip_.get(0),
                                         LowerRight_.get(1)};

    ContainingRectangle_ = Box<DIM, real_number>(LowerLeft_, UpperRight_);
}

bool TriangleEqui::isInside(Point<DIM, real_number> P)
{
    return (ContainingRectangle_.isInside(P) && !isAvobeLine(TriangleTip_, UpperRight_, P, params_.dp) && !isBelowLine(TriangleTip_, LowerRight_, P, params_.dp));
}

void TriangleEqui::AddObstacle(particles &vd)
{
    real_number dx_self = params_.dp / refine_factor;
    const int N_wall = ceil(SideLength_ / dx_self);
    const real_number dxwall = SideLength_ / N_wall;
    Point<DIM, real_number> Yoffset = {0.0, dxwall};

    // Right wall
    AddFlatWallNewBC(vd, 0, N_wall + 1, LowerRight_, Yoffset, dxwall, Centre_, LinearVelocity_, params_, OBSTACLE, AngularVelocity_);

    // Inclined walls

    const real_number cos_theta = sqrtf(3.0) / 2.0f;
    const real_number sin_theta = 1 / 2.0f;
    const Point<DIM, real_number> Diagoffset{dxwall * cos_theta, dxwall * sin_theta};
    Point<DIM, real_number> DiagoffsetNW = Diagoffset;
    DiagoffsetNW.get(0) = -1.0 * DiagoffsetNW.get(0);
    //  Hypothenuse upper wall
    AddFlatWallNewBC(vd, 1, N_wall + 1, UpperRight_, -1.0 * Diagoffset, dxwall, Centre_, LinearVelocity_, params_, OBSTACLE, AngularVelocity_);
    // Hypothenuse lower wall
    AddFlatWallNewBC(vd, 1, N_wall, LowerRight_, DiagoffsetNW, dxwall, Centre_, LinearVelocity_, params_, OBSTACLE, AngularVelocity_);
}
