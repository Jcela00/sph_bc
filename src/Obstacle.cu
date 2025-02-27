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
        vd.template getLastProp<vd9_volume>()[1] = 0.0;
        vd.template getLastProp<vd9_volume>()[2] = 0.0;

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
                         const real_number given_kappa,
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
        vd.template getLastProp<vd9_volume>()[1] = given_kappa;
        vd.template getLastProp<vd9_volume>()[2] = 0.0;

        vd.template getLastProp<vd10_omega>() = obstacle_omega;
    }
}

void AddFlatWallModNewBC3D(particles &vd,
                           const int k0,
                           const int kmax,
                           const int j0,
                           const int jmax,
                           const Point<DIM, real_number> Corner,
                           const Point<DIM, real_number> UnitOffsetK,
                           const Point<DIM, real_number> UnitOffsetJ,
                           const real_number dx,
                           const Point<DIM, real_number> obstacle_centre,
                           const Point<DIM, real_number> obstacle_velocity,
                           const Parameters &arg_p,
                           const size_t particle_type,
                           const Point<DIM, real_number> given_normal,
                           const real_number given_kappa,
                           const real_number obstacle_omega)
{

    for (int k = k0; k < kmax; k++)
    {
        for (int j = j0; j < jmax; j++)
        {

            Point<DIM, real_number> WallPos = Corner + k * UnitOffsetK + j * UnitOffsetJ;

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
            vd.template getLastProp<vd9_volume>()[1] = given_kappa;
            vd.template getLastProp<vd9_volume>()[2] = 0.0;

            vd.template getLastProp<vd10_omega>() = obstacle_omega;
        }
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
bool EmptyObstacle::isInsidePlusTol(Point<DIM, real_number> P)
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
bool CylinderObstacle::isInsidePlusTol(Point<DIM, real_number> P)
{
    bool is_inside = false;
    real_number tol = 0.25 * params_.dp;
    real_number tmp = (P.get(0) - Centre_.get(0)) * (P.get(0) - Centre_.get(0)) + (P.get(1) - Centre_.get(1)) * (P.get(1) - Centre_.get(1));
    if (tmp <= (Radius_ + tol) * (Radius_ + tol))
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
        if (k == 0)
            vd.template getLastProp<vd3_drho>() = 10.0;
        else
            vd.template getLastProp<vd3_drho>() = 0.0;

        vd.template getLastProp<vd9_volume>()[0] = dxwall;
        vd.template getLastProp<vd9_volume>()[1] = 0.0;
        vd.template getLastProp<vd9_volume>()[2] = 0.0;

        vd.template getLastProp<vd10_omega>() = AngularVelocity_;
        theta += dtheta;
    }
}

// SphereObstacle class
SphereObstacle::SphereObstacle(real_number Radius,
                               Point<DIM, real_number> centre,
                               const Parameters &p,
                               Point<DIM, real_number> vel,
                               real_number omega,
                               real_number rf,
                               bool autoNormals) : Obstacle(centre, p, vel, omega, rf), Radius_(Radius), autoNormals_(autoNormals) {};

bool SphereObstacle::isInside(Point<DIM, real_number> P)
{
    bool is_inside = false;
    real_number tmp = (P.get(0) - Centre_.get(0)) * (P.get(0) - Centre_.get(0)) + (P.get(1) - Centre_.get(1)) * (P.get(1) - Centre_.get(1)) + (P.get(2) - Centre_.get(2)) * (P.get(2) - Centre_.get(2));
    if (tmp <= Radius_ * Radius_)
    {
        is_inside = true;
    }
    return is_inside;
}
bool SphereObstacle::isInsidePlusTol(Point<DIM, real_number> P)
{
    return false;
}

void SphereObstacle::AddObstacle(particles &vd)
{
    const real_number dx_self = params_.dp / refine_factor;
    const real_number sphere_area = 4.0 * M_PI * Radius_ * Radius_;
    const real_number area_element = dx_self * dx_self;

    const unsigned int Npoints = ceil(sphere_area / area_element);
    const real_number NpointsReal = static_cast<real_number>(Npoints);

    const real_number area_element_actual = sphere_area / NpointsReal; // not exactly area_element since sphere area need not be an exact multiple of area_element

    Point<DIM, real_number> SpherePositon;
    const real_number golden_ratio = (1.0 + sqrt(5.0)) / 2.0;

    for (int k = 1; k < Npoints + 1; k++)
    {
        real_number latitude = asin((2.0 * k - NpointsReal - 1.0) / NpointsReal);
        real_number longitude = 2.0 * M_PI * k / golden_ratio;

        SpherePositon.get(0) = Radius_ * cos(longitude) * cos(latitude);
        SpherePositon.get(1) = Radius_ * sin(longitude) * cos(latitude);
        SpherePositon.get(2) = Radius_ * sin(latitude);

        vd.add();
        vd.template getLastProp<vd0_type>() = OBSTACLE;
        for (int xyz = 0; xyz < DIM; xyz++)
        {
            vd.getLastPos()[xyz] = Centre_.get(xyz) + SpherePositon.get(xyz); // + ((real_number)rand() / RAND_MAX - 0.5) * dp;
            vd.template getLastProp<vd6_force>()[xyz] = 0.0;
            vd.template getLastProp<vd7_force_t>()[xyz] = Centre_.get(xyz);
            vd.template getLastProp<vd5_velocity_t>()[xyz] = 0.0;
            if (autoNormals_)
                vd.template getLastProp<vd8_normal>()[xyz] = 100.0 * SpherePositon.get(xyz);
            else
                vd.template getLastProp<vd8_normal>()[xyz] = 0.0;

            vd.template getLastProp<vd4_velocity>()[xyz] = LinearVelocity_.get(xyz);
        }

        vd.template getLastProp<vd2_pressure>() = 0.0;
        vd.template getLastProp<vd1_rho>() = params_.rho0;
        vd.template getLastProp<vd3_drho>() = 0.0;
        vd.template getLastProp<vd9_volume>()[0] = area_element_actual;
        vd.template getLastProp<vd9_volume>()[1] = 0.0;
        vd.template getLastProp<vd9_volume>()[2] = 0.0;
        vd.template getLastProp<vd10_omega>() = AngularVelocity_;
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
bool EllipticObstacle::isInsidePlusTol(Point<DIM, real_number> P)
{
    return false;
}

void EllipticObstacle::AddObstacle(particles &vd)
{
    // first we need to determine the perimeter of the ellipse
    // we can use c++ std::comp_ellint_2 function
    // first compute eccentricity
    const real_number e = sqrt(1.0 - Minor_ * Minor_ / Major_ / Major_);
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
    vd.template getLastProp<vd9_volume>()[1] = 0.0;
    vd.template getLastProp<vd9_volume>()[2] = 0.0;
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
            vd.template getLastProp<vd9_volume>()[1] = 0.0;
            vd.template getLastProp<vd9_volume>()[2] = 0.0;
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
    LowerLeft_ = Point<DIM, real_number>{Centre_.get(0) - BaseLength_ / 2.0,
                                         Centre_.get(1) - HeigthLength_ / 2.0};

    UpperRight_ = Point<DIM, real_number>{Centre_.get(0) + BaseLength_ / 2.0,
                                          Centre_.get(1) + HeigthLength_ / 2.0};

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
bool RectangleObstacle::isInsidePlusTol(Point<DIM, real_number> P)
{
    return false;
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
    LowerLeft_ = Point<DIM, real_number>{Centre_.get(0) - 2.0 * BaseLength_ / 3.0,
                                         Centre_.get(1) - HeigthLength_ / 3.0};

    LowerRight_ = Point<DIM, real_number>{Centre_.get(0) + BaseLength_ / 3.0,
                                          Centre_.get(1) - HeigthLength_ / 3.0};

    UpperRight_ = Point<DIM, real_number>{Centre_.get(0) + BaseLength_ / 3.0,
                                          Centre_.get(1) + 2.0 * HeigthLength_ / 3.0};

    ContainingRectangle_ = Box<DIM, real_number>(LowerLeft_, UpperRight_);
}

bool TriangleObstacle::isInside(Point<DIM, real_number> P)
{
    return (ContainingRectangle_.isInside(P) && !isAvobeLine(LowerLeft_, UpperRight_, P, params_.dp));
}
bool TriangleObstacle::isInsidePlusTol(Point<DIM, real_number> P)
{
    return false;
}

void TriangleObstacle::AddObstacle(particles &vd)
{
    const real_number dx_self = params_.dp / refine_factor;
    // Lower wall
    const int N_bottom = ceil(BaseLength_ / dx_self);
    const real_number dxwall_bottom = BaseLength_ / N_bottom;
    const Point<DIM, real_number> Xoffset = {dxwall_bottom, 0.0};
    AddFlatWallModNewBC(vd, 0, 1, LowerLeft_, Xoffset, dxwall_bottom, Centre_, LinearVelocity_, params_, OBSTACLE, {-1.0, 0.0}, 0.0, AngularVelocity_);
    AddFlatWallNewBC(vd, 1, N_bottom + 1, LowerLeft_, Xoffset, dxwall_bottom, Centre_, LinearVelocity_, params_, OBSTACLE, AngularVelocity_);

    // Right wall
    const int N_right = ceil(HeigthLength_ / dx_self);
    const real_number dxwall_right = HeigthLength_ / N_right;
    const Point<DIM, real_number> Yoffset = {0.0, dxwall_right};
    AddFlatWallNewBC(vd, 1, N_right + 1, LowerRight_, Yoffset, dxwall_right, Centre_, LinearVelocity_, params_, OBSTACLE, AngularVelocity_);

    //  Hypothenuse wall
    // We want particles spaced roughly by dp
    const real_number HypothenuseLength = sqrt(BaseLength_ * BaseLength_ + HeigthLength_ * HeigthLength_);
    const int Ndiag = ceil(HypothenuseLength / dx_self);       // integer number of particles that can fit in the diagonal
    const real_number dxwall_diag = HypothenuseLength / Ndiag; // actual spacing between particles ( close to dp but not exactly)
    const real_number sin_theta = HeigthLength_ / HypothenuseLength;
    const real_number cos_theta = BaseLength_ / HypothenuseLength;
    const Point<DIM, real_number> Diagoffset{dxwall_diag * cos_theta, dxwall_diag * sin_theta};

    AddFlatWallNewBC(vd, 1, Ndiag, UpperRight_, -1.0 * Diagoffset, dxwall_diag, Centre_, LinearVelocity_, params_, OBSTACLE, AngularVelocity_);
}

// TriangleTestObstacle class
TriangleTestObstacle::TriangleTestObstacle(Point<DIM, real_number> centre,
                                           const Parameters &p,
                                           real_number BaseLength,
                                           real_number HeigthLength,
                                           Point<DIM, real_number> vel,
                                           real_number omega,
                                           real_number rf) : Obstacle(centre, p, vel, omega, rf), BaseLength_(BaseLength), HeigthLength_(HeigthLength)
{
    LowerLeft_ = Point<DIM, real_number>{Centre_.get(0) - 2.0 * BaseLength_ / 3.0,
                                         Centre_.get(1) - HeigthLength_ / 3.0};

    LowerRight_ = Point<DIM, real_number>{Centre_.get(0) + BaseLength_ / 3.0,
                                          Centre_.get(1) - HeigthLength_ / 3.0};

    UpperRight_ = Point<DIM, real_number>{Centre_.get(0) + BaseLength_ / 3.0,
                                          Centre_.get(1) + 2.0 * HeigthLength_ / 3.0};

    ContainingRectangle_ = Box<DIM, real_number>(LowerLeft_, UpperRight_);
}

bool TriangleTestObstacle::isInside(Point<DIM, real_number> P)
{
    return (ContainingRectangle_.isInside(P) && !isAvobeLine(LowerLeft_, UpperRight_, P, params_.dp));
}
bool TriangleTestObstacle::isInsidePlusTol(Point<DIM, real_number> P)
{
    return false;
}

void TriangleTestObstacle::AddObstacle(particles &vd)
{

    // OPTION TWO
    const real_number dx_self = params_.dp / refine_factor;
    // Lower wall
    const int N_bottom = ceil(BaseLength_ / dx_self);
    const real_number dxwall_bottom = BaseLength_ / N_bottom;
    const Point<DIM, real_number> Xoffset = {dxwall_bottom, 0.0};
    // Lower walls except corners
    AddFlatWallModNewBC(vd, 1, N_bottom, LowerLeft_, Xoffset, dxwall_bottom, Centre_, LinearVelocity_, params_, OBSTACLE, {0.0, -10.0}, 0.0, AngularVelocity_);
    // Lower right corner
    AddFlatWallModNewBC(vd, N_bottom, N_bottom + 1, LowerLeft_, Xoffset, dxwall_bottom, Centre_, LinearVelocity_, params_, OBSTACLE, {0.5, -0.5}, 0.0, AngularVelocity_);

    // Right wall except corners
    const int N_right = ceil(HeigthLength_ / dx_self);
    const real_number dxwall_right = HeigthLength_ / N_right;
    const Point<DIM, real_number> Yoffset = {0.0, dxwall_right};
    AddFlatWallModNewBC(vd, 1, N_right, LowerRight_, Yoffset, dxwall_right, Centre_, LinearVelocity_, params_, OBSTACLE, {10.0, 0.0}, 0.0, AngularVelocity_);

    // Hypothenuse wall
    // We want particles spaced roughly by dp
    const real_number HypothenuseLength = sqrt(BaseLength_ * BaseLength_ + HeigthLength_ * HeigthLength_);
    const int Ndiag = ceil(HypothenuseLength / dx_self);       // integer number of particles that can fit in the diagonal
    const real_number dxwall_diag = HypothenuseLength / Ndiag; // actual spacing between particles ( close to dp but not exactly)
    const real_number sin_theta = HeigthLength_ / HypothenuseLength;
    const real_number cos_theta = BaseLength_ / HypothenuseLength;
    const Point<DIM, real_number> Diagoffset{dxwall_diag * cos_theta, dxwall_diag * sin_theta};
    const Point<DIM, real_number> DiagNormal{-10.0 * sin_theta, 10.0 * cos_theta};
    AddFlatWallModNewBC(vd, 1, Ndiag, UpperRight_, -1.0 * Diagoffset, dxwall_diag, Centre_, LinearVelocity_, params_, OBSTACLE, DiagNormal, 0.0, AngularVelocity_);

    // const Point<DIM, real_number> NormalURCorner = {cos_theta, sin_theta};

    AddFlatWallModNewBC(vd, 0, 1, UpperRight_, -1.0 * Diagoffset, dxwall_diag, Centre_, LinearVelocity_, params_, OBSTACLE, {0.5, 0.5}, 0.0, AngularVelocity_);
    AddFlatWallModNewBC(vd, Ndiag, Ndiag + 1, UpperRight_, -1.0 * Diagoffset, dxwall_diag, Centre_, LinearVelocity_, params_, OBSTACLE, {-0.5, 0.0}, 0.0, AngularVelocity_);
}

// Triangle equi class
TriangleEqui::TriangleEqui(Point<DIM, real_number> centre,
                           const Parameters &p,
                           real_number sidelength,
                           Point<DIM, real_number> vel,
                           real_number omega,
                           real_number rf) : Obstacle(centre, p, vel, omega, rf), SideLength_(sidelength)
{

    UpperRight_ = Point<DIM, real_number>{Centre_.get(0) + sqrt(3.0) * SideLength_ / 6.0,
                                          Centre_.get(1) + SideLength_ / 2.0};
    LowerRight_ = Point<DIM, real_number>{Centre_.get(0) + sqrt(3.0) * SideLength_ / 6.0,
                                          Centre_.get(1) - SideLength_ / 2.0};
    TriangleTip_ = Point<DIM, real_number>{Centre_.get(0) - sqrt(3.0) * SideLength_ / 3.0,
                                           Centre_.get(1)};

    LowerLeft_ = Point<DIM, real_number>{TriangleTip_.get(0),
                                         LowerRight_.get(1)};

    ContainingRectangle_ = Box<DIM, real_number>(LowerLeft_, UpperRight_);
}

bool TriangleEqui::isInside(Point<DIM, real_number> P)
{
    return (ContainingRectangle_.isInside(P) && !isAvobeLine(TriangleTip_, UpperRight_, P, params_.dp) && !isBelowLine(TriangleTip_, LowerRight_, P, params_.dp));
}
bool TriangleEqui::isInsidePlusTol(Point<DIM, real_number> P)
{
    return false;
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

    const real_number cos_theta = sqrt(3.0) / 2.0;
    const real_number sin_theta = 1 / 2.0;
    const Point<DIM, real_number> Diagoffset{dxwall * cos_theta, dxwall * sin_theta};
    Point<DIM, real_number> DiagoffsetNW = Diagoffset;
    DiagoffsetNW.get(0) = -1.0 * DiagoffsetNW.get(0);
    //  Hypothenuse upper wall
    AddFlatWallNewBC(vd, 1, N_wall + 1, UpperRight_, -1.0 * Diagoffset, dxwall, Centre_, LinearVelocity_, params_, OBSTACLE, AngularVelocity_);
    // Hypothenuse lower wall
    AddFlatWallNewBC(vd, 1, N_wall, LowerRight_, DiagoffsetNW, dxwall, Centre_, LinearVelocity_, params_, OBSTACLE, AngularVelocity_);
}

// CurveObstacle class
CurveObstacle::CurveObstacle(real_number a,
                             real_number b,
                             real_number k,
                             real_number m,
                             Point<DIM, real_number> centre,
                             const Parameters &p,
                             Point<DIM, real_number> vel,
                             real_number omega,
                             real_number rf) : Obstacle(centre, p, vel, omega, rf), a_(a), b_(b), k_(k), m_(m) {}

void CurveObstacle::AddObstacle(particles &vd)
{
    // first compute total perimeter of the curve

    real_number nsteps = 1e7;
    real_number dtheta = 2.0 * M_PI / nsteps;
    real_number arc_length = 0.0;
    real_number theta = 0.0;

    real_number R = parametricRadius(theta);
    real_number dR = parametricRadiusDerivative(theta);

    real_number f0 = sqrt(R * R + dR * dR);
    real_number fprev = f0;
    theta += dtheta;

    for (int i = 0; i < nsteps; i++)
    {
        real_number R = parametricRadius(theta);
        real_number dR = parametricRadiusDerivative(theta);
        real_number f = sqrt(R * R + dR * dR);
        real_number ds = 0.5 * (f + fprev) * dtheta;
        fprev = f;
        arc_length += ds;
        theta += dtheta;
    }

    printf("perimeter of the curve is %f\n", arc_length);

    // now we want to discretize the curve with particles
    const real_number dx_self = params_.dp / refine_factor;
    const int Npoints = ceil(arc_length / dx_self);

    const real_number darc = arc_length / (real_number)Npoints;
    arc_length = 0.0;

    R = parametricRadius(0.0);
    // add first particle at theta = 0
    vd.add();

    vd.template getLastProp<vd0_type>() = OBSTACLE;
    vd.getLastPos()[0] = Centre_.get(0) + R * cos(0.0);
    vd.getLastPos()[1] = Centre_.get(1) + R * sin(0.0);

    for (int xyz = 0; xyz < DIM; xyz++)
    {
        vd.template getLastProp<vd6_force>()[xyz] = 0.0;
        vd.template getLastProp<vd7_force_t>()[xyz] = Centre_.get(xyz);
        vd.template getLastProp<vd5_velocity_t>()[xyz] = 0.0;
        vd.template getLastProp<vd8_normal>()[xyz] = 0.0;
        vd.template getLastProp<vd4_velocity>()[xyz] = LinearVelocity_.get(xyz);
    }

    vd.template getLastProp<vd2_pressure>() = 0.0;
    vd.template getLastProp<vd1_rho>() = params_.rho0;
    vd.template getLastProp<vd3_drho>() = 0.0;
    vd.template getLastProp<vd9_volume>()[0] = darc;
    vd.template getLastProp<vd9_volume>()[1] = 0.0;
    vd.template getLastProp<vd9_volume>()[2] = 0.0;
    vd.template getLastProp<vd10_omega>() = AngularVelocity_;

    // add the rest of the particles
    theta = dtheta;
    fprev = f0;
    real_number accumulated_arc = 0.0;

    for (int i = 0; i < nsteps; i++)
    {
        real_number R = parametricRadius(theta);
        real_number dR = parametricRadiusDerivative(theta);
        real_number f = sqrt(R * R + dR * dR);
        real_number ds = 0.5 * (f + fprev) * dtheta;
        fprev = f;
        accumulated_arc += ds;
        theta += dtheta;

        if (accumulated_arc > darc)
        {
            vd.add();
            vd.template getLastProp<vd0_type>() = OBSTACLE;
            vd.getLastPos()[0] = Centre_.get(0) + R * cos(theta);
            vd.getLastPos()[1] = Centre_.get(1) + R * sin(theta);

            for (int xyz = 0; xyz < DIM; xyz++)
            {
                vd.template getLastProp<vd6_force>()[xyz] = 0.0;
                vd.template getLastProp<vd7_force_t>()[xyz] = Centre_.get(xyz);
                vd.template getLastProp<vd5_velocity_t>()[xyz] = 0.0;
                vd.template getLastProp<vd8_normal>()[xyz] = 0.0;
                vd.template getLastProp<vd4_velocity>()[xyz] = LinearVelocity_.get(xyz);
            }

            vd.template getLastProp<vd2_pressure>() = 0.0;
            vd.template getLastProp<vd1_rho>() = params_.rho0;
            vd.template getLastProp<vd3_drho>() = 0.0;
            vd.template getLastProp<vd9_volume>()[0] = darc;
            vd.template getLastProp<vd9_volume>()[1] = 0.0;
            vd.template getLastProp<vd9_volume>()[2] = 0.0;
            vd.template getLastProp<vd10_omega>() = AngularVelocity_;

            accumulated_arc = 0.0;
        }
    }
}

bool CurveObstacle::isInside(Point<DIM, real_number> P)
{
    // get point minus the centre
    Point<DIM, real_number> Paux = P - Centre_;

    real_number Rpoint = sqrt(Paux.get(0) * Paux.get(0) + Paux.get(1) * Paux.get(1));
    real_number TPoint = atan2(Paux.get(1), Paux.get(0));

    real_number Rcurve = parametricRadius(TPoint);

    return (Rpoint <= Rcurve);
}

bool CurveObstacle::isInsidePlusTol(Point<DIM, real_number> P)
{
    return false;
}

real_number CurveObstacle::parametricRadius(real_number theta)
{
    return a_ * (1 + b_ * std::pow(cos(k_ * theta), m_));
}

real_number CurveObstacle::parametricRadiusDerivative(real_number theta)
{
    return -a_ * b_ * m_ * k_ * std::pow(cos(k_ * theta), m_ - 1) * sin(k_ * theta);
}

// EpiCycloid Class
EpiCycloid::EpiCycloid(real_number R,
                       real_number k,
                       Point<DIM, real_number> centre,
                       const Parameters &p,
                       Point<DIM, real_number> vel,
                       real_number omega,
                       real_number rf) : Obstacle(centre, p, vel, omega, rf), R_(R), k_(k)
{
    r_ = R_ / k_;
}

void EpiCycloid::AddObstacle(particles &vd)
{
    // first compute total perimeter of the curve

    real_number nsteps = 1e7;
    real_number dtheta = 2.0 * M_PI / nsteps;
    real_number arc_length = 0.0;
    real_number theta = 0.0;

    Point<DIM, real_number> dP = ParametricCoordsDerivative(theta);

    real_number f0 = sqrt(dP.get(0) * dP.get(0) + dP.get(1) * dP.get(1));
    real_number fprev = f0;
    theta += dtheta;

    for (int i = 0; i < nsteps; i++)
    {
        Point<DIM, real_number> dP = ParametricCoordsDerivative(theta);
        real_number f = sqrt(dP.get(0) * dP.get(0) + dP.get(1) * dP.get(1));
        real_number ds = 0.5 * (f + fprev) * dtheta;
        fprev = f;
        arc_length += ds;
        theta += dtheta;
    }
    printf("perimeter of the Epicycloid is %f\n", arc_length);

    // now we want to discretize the curve with particles
    const real_number dx_self = params_.dp / refine_factor;
    const int Npoints = k_ * ceil(arc_length / k_ / dx_self);

    const real_number darc = arc_length / (real_number)Npoints;

    AddLobes(vd, Npoints, darc, nsteps);
}

void EpiCycloid::AddLobes(particles &vd, const int Npoints, const real_number darc, const int nsteps)
{
    real_number dtheta = (2 * M_PI / k_) / nsteps;

    // place lobes one by one

    for (int kk = 0; kk < k_; kk++)
    {
        // add the first particle manually
        real_number theta0 = kk * 2 * M_PI / k_;
        Point<DIM, real_number> P0 = ParametricCoords(theta0);
        Point<DIM, real_number> dP0 = ParametricCoordsDerivative(theta0);
        real_number f0 = sqrt(dP0.get(0) * dP0.get(0) + dP0.get(1) * dP0.get(1));

        vd.add();

        vd.template getLastProp<vd0_type>() = OBSTACLE;
        vd.getLastPos()[0] = Centre_.get(0) + P0.get(0);
        vd.getLastPos()[1] = Centre_.get(1) + P0.get(1);

        // normal in - radial direction
        real_number nx = R_ * cos(theta0);
        real_number ny = R_ * sin(theta0);
        vd.template getLastProp<vd8_normal>()[0] = -1000.0 * nx; // -0.5 * 1-;
        vd.template getLastProp<vd8_normal>()[1] = -1000.0 * ny; // 0.5 * dx;-y
        for (int xyz = 0; xyz < DIM; xyz++)
        {
            vd.template getLastProp<vd6_force>()[xyz] = 0.0;
            vd.template getLastProp<vd7_force_t>()[xyz] = Centre_.get(xyz);
            vd.template getLastProp<vd5_velocity_t>()[xyz] = 0.0;
            vd.template getLastProp<vd4_velocity>()[xyz] = LinearVelocity_.get(xyz);
        }

        vd.template getLastProp<vd2_pressure>() = 0.0;
        vd.template getLastProp<vd1_rho>() = params_.rho0;
        vd.template getLastProp<vd3_drho>() = 0.0;
        vd.template getLastProp<vd9_volume>()[0] = darc;
        vd.template getLastProp<vd9_volume>()[1] = 0.0;
        vd.template getLastProp<vd9_volume>()[2] = 0.0;
        vd.template getLastProp<vd10_omega>() = AngularVelocity_;

        // add the rest of the particles
        real_number theta = theta0 + dtheta;
        real_number fprev = f0;
        real_number accumulated_arc = 0.0;

        for (int i = 0; i < nsteps; i++)
        {
            // take a step in theta
            Point<DIM, real_number> P = ParametricCoords(theta);
            Point<DIM, real_number> dP = ParametricCoordsDerivative(theta);

            // integrate arc length
            real_number f = sqrt(dP.get(0) * dP.get(0) + dP.get(1) * dP.get(1));
            real_number ds = 0.5 * (f + fprev) * dtheta;

            fprev = f;
            accumulated_arc += ds;
            theta += dtheta;

            if (accumulated_arc > darc)
            {
                vd.add();
                vd.template getLastProp<vd0_type>() = OBSTACLE;
                vd.getLastPos()[0] = Centre_.get(0) + P.get(0);
                vd.getLastPos()[1] = Centre_.get(1) + P.get(1);

                // real_number dx = -r_ * (k_ + 1) * sin(theta) + r_ * (k_ + 1) * sin((k_ + 1) * theta);
                // real_number dy = r_ * (k_ + 1) * cos(theta) - r_ * (k_ + 1) * cos((k_ + 1) * theta);

                vd.template getLastProp<vd8_normal>()[0] = 0.0; // P.get(0) - Centre_.get(0);
                vd.template getLastProp<vd8_normal>()[1] = 0.0; // P.get(1) - Centre_.get(1);

                for (int xyz = 0; xyz < DIM; xyz++)
                {
                    vd.template getLastProp<vd6_force>()[xyz] = 0.0;
                    vd.template getLastProp<vd7_force_t>()[xyz] = Centre_.get(xyz);
                    vd.template getLastProp<vd5_velocity_t>()[xyz] = 0.0;
                    vd.template getLastProp<vd4_velocity>()[xyz] = LinearVelocity_.get(xyz);
                }

                vd.template getLastProp<vd2_pressure>() = 0.0;
                vd.template getLastProp<vd1_rho>() = params_.rho0;
                vd.template getLastProp<vd3_drho>() = 0.0;
                vd.template getLastProp<vd9_volume>()[0] = darc;
                vd.template getLastProp<vd9_volume>()[1] = 0.0;
                vd.template getLastProp<vd9_volume>()[2] = 0.0;
                vd.template getLastProp<vd10_omega>() = AngularVelocity_;

                accumulated_arc = 0.0;
            }
        }
    }
}

bool EpiCycloid::isInside(Point<DIM, real_number> P)
{

    Point<DIM, real_number> Paux = P - Centre_;
    real_number Rpoint = sqrt(Paux.get(0) * Paux.get(0) + Paux.get(1) * Paux.get(1));

    real_number Rmax = (5.0 / 3.0) * R_;
    real_number Rmin = R_;
    if (Rpoint > Rmax)
    {
        return false;
    }
    else if (Rpoint < Rmin)
    {
        return true;
    }
    else // more expensive check, winding number
    {
        real_number ngrid = 1e3;
        real_number dtheta = 2.0 * M_PI / ngrid;

        real_number accumulated_angle = 0.0;

        for (int i = 0; i < ngrid; i++)
        {
            real_number theta = i * dtheta;
            Point<DIM, real_number> Pcurve = ParametricCoords(theta);
            Point<DIM, real_number> Pcurve_next = ParametricCoords(theta + dtheta);

            // angle = atan2(ycurvenext-ypoint,xcurvenext-xpoint) - atan2(ycurve-ypoint,xcurve-xpoint)
            real_number angle = atan2(Pcurve_next.get(1) - Paux.get(1), Pcurve_next.get(0) - Paux.get(0)) - atan2(Pcurve.get(1) - Paux.get(1), Pcurve.get(0) - Paux.get(0));
            // normalize to -pi,pi
            if (angle > M_PI)
            {
                angle -= 2.0 * M_PI;
            }
            else if (angle < -M_PI)
            {
                angle += 2.0 * M_PI;
            }
            accumulated_angle += angle;
        }
        real_number tol = 0.1 * 2.0 * M_PI;

        if (abs(accumulated_angle) > tol)
        {
            return true;
        }
        else
        {
            return false;
        }
    }
}

bool EpiCycloid::isInsidePlusTol(Point<DIM, real_number> P)
{
    return false;
}

Point<DIM, real_number> EpiCycloid::ParametricCoords(real_number theta)
{
    real_number x = r_ * (k_ + 1) * cos(theta) - r_ * cos((k_ + 1) * theta);
    real_number y = r_ * (k_ + 1) * sin(theta) - r_ * sin((k_ + 1) * theta);
    return Point<DIM, real_number>{x, y};
}

Point<DIM, real_number> EpiCycloid::ParametricCoordsDerivative(real_number theta)
{
    real_number dx = -r_ * (k_ + 1) * sin(theta) + r_ * (k_ + 1) * sin((k_ + 1) * theta);
    real_number dy = r_ * (k_ + 1) * cos(theta) - r_ * (k_ + 1) * cos((k_ + 1) * theta);
    return Point<DIM, real_number>{dx, dy};
}

// HipoCycloid Class
HipoCycloid::HipoCycloid(real_number R,
                         real_number k,
                         Point<DIM, real_number> centre,
                         const Parameters &p,
                         Point<DIM, real_number> vel,
                         real_number omega,
                         real_number rf) : Obstacle(centre, p, vel, omega, rf), R_(R), k_(k)
{
    r_ = R_ / k_;
}

void HipoCycloid::AddObstacle(particles &vd)
{
    // first compute total perimeter of the curve

    real_number nsteps = 1e7;
    real_number dtheta = 2.0 * M_PI / nsteps;
    real_number arc_length = 0.0;
    real_number theta = 0.0;

    Point<DIM, real_number> dP = ParametricCoordsDerivative(theta);

    real_number f0 = sqrt(dP.get(0) * dP.get(0) + dP.get(1) * dP.get(1));
    real_number fprev = f0;
    theta += dtheta;

    for (int i = 0; i < nsteps; i++)
    {
        Point<DIM, real_number> dP = ParametricCoordsDerivative(theta);
        real_number f = sqrt(dP.get(0) * dP.get(0) + dP.get(1) * dP.get(1));
        real_number ds = 0.5 * (f + fprev) * dtheta;
        fprev = f;
        arc_length += ds;
        theta += dtheta;
    }
    printf("perimeter of the HipoCycloid is %f\n", arc_length);

    // now we want to discretize the curve with particles
    const real_number dx_self = params_.dp / refine_factor;
    const int Npoints = k_ * ceil(arc_length / k_ / dx_self);

    const real_number darc = arc_length / (real_number)Npoints;

    AddLobes(vd, Npoints, darc, nsteps);
}

void HipoCycloid::AddLobes(particles &vd, const int Npoints, const real_number darc, const int nsteps)
{
    real_number dtheta = (2 * M_PI / k_) / nsteps;

    // place lobes one by one

    for (int kk = 0; kk < k_; kk++)
    {
        // add the first particle manually
        real_number theta0 = kk * 2 * M_PI / k_;
        Point<DIM, real_number> P0 = ParametricCoords(theta0);
        Point<DIM, real_number> dP0 = ParametricCoordsDerivative(theta0);
        real_number f0 = sqrt(dP0.get(0) * dP0.get(0) + dP0.get(1) * dP0.get(1));

        vd.add();

        vd.template getLastProp<vd0_type>() = OBSTACLE;
        vd.getLastPos()[0] = Centre_.get(0) + P0.get(0);
        vd.getLastPos()[1] = Centre_.get(1) + P0.get(1);

        // normal in - radial direction
        real_number nx = R_ * cos(theta0);
        real_number ny = R_ * sin(theta0);
        vd.template getLastProp<vd8_normal>()[0] = 1000.0 * nx; // -0.5 * 1-;
        vd.template getLastProp<vd8_normal>()[1] = 1000.0 * ny; // 0.5 * dx;-y
        for (int xyz = 0; xyz < DIM; xyz++)
        {
            vd.template getLastProp<vd6_force>()[xyz] = 0.0;
            vd.template getLastProp<vd7_force_t>()[xyz] = Centre_.get(xyz);
            vd.template getLastProp<vd5_velocity_t>()[xyz] = 0.0;
            vd.template getLastProp<vd4_velocity>()[xyz] = LinearVelocity_.get(xyz);
        }

        vd.template getLastProp<vd2_pressure>() = 0.0;
        vd.template getLastProp<vd1_rho>() = params_.rho0;
        vd.template getLastProp<vd3_drho>() = 0.0;
        vd.template getLastProp<vd9_volume>()[0] = darc;
        vd.template getLastProp<vd9_volume>()[1] = 0.0;
        vd.template getLastProp<vd9_volume>()[2] = 0.0;
        vd.template getLastProp<vd10_omega>() = AngularVelocity_;

        // add the rest of the particles
        real_number theta = theta0 + dtheta;
        real_number fprev = f0;
        real_number accumulated_arc = 0.0;

        for (int i = 0; i < nsteps; i++)
        {
            // take a step in theta
            Point<DIM, real_number> P = ParametricCoords(theta);
            Point<DIM, real_number> dP = ParametricCoordsDerivative(theta);

            // integrate arc length
            real_number f = sqrt(dP.get(0) * dP.get(0) + dP.get(1) * dP.get(1));
            real_number ds = 0.5 * (f + fprev) * dtheta;

            fprev = f;
            accumulated_arc += ds;
            theta += dtheta;

            if (accumulated_arc > darc)
            {
                vd.add();
                vd.template getLastProp<vd0_type>() = OBSTACLE;
                vd.getLastPos()[0] = Centre_.get(0) + P.get(0);
                vd.getLastPos()[1] = Centre_.get(1) + P.get(1);

                // real_number dx = -r_ * (k_ + 1) * sin(theta) + r_ * (k_ + 1) * sin((k_ + 1) * theta);
                // real_number dy = r_ * (k_ + 1) * cos(theta) - r_ * (k_ + 1) * cos((k_ + 1) * theta);
                // real_number nx = R_ * cos(theta0 + M_PI / k_);
                // real_number ny = R_ * sin(theta0 + M_PI / k_);
                real_number nx = R_ * cos(theta);
                real_number ny = R_ * sin(theta);
                vd.template getLastProp<vd8_normal>()[0] = 0.5 * nx; // P.get(0) - Centre_.get(0);
                vd.template getLastProp<vd8_normal>()[1] = 0.5 * ny; // P.get(1) - Centre_.get(1);

                for (int xyz = 0; xyz < DIM; xyz++)
                {
                    vd.template getLastProp<vd6_force>()[xyz] = 0.0;
                    vd.template getLastProp<vd7_force_t>()[xyz] = Centre_.get(xyz);
                    vd.template getLastProp<vd5_velocity_t>()[xyz] = 0.0;
                    vd.template getLastProp<vd4_velocity>()[xyz] = LinearVelocity_.get(xyz);
                }

                vd.template getLastProp<vd2_pressure>() = 0.0;
                vd.template getLastProp<vd1_rho>() = params_.rho0;
                vd.template getLastProp<vd3_drho>() = 0.0;
                vd.template getLastProp<vd9_volume>()[0] = darc;
                vd.template getLastProp<vd9_volume>()[1] = 0.0;
                vd.template getLastProp<vd9_volume>()[2] = 0.0;
                vd.template getLastProp<vd10_omega>() = AngularVelocity_;

                accumulated_arc = 0.0;
            }
        }
    }
}

bool HipoCycloid::isInside(Point<DIM, real_number> P)
{

    Point<DIM, real_number> Paux = P - Centre_;
    real_number Rpoint = sqrt(Paux.get(0) * Paux.get(0) + Paux.get(1) * Paux.get(1));

    real_number Rmax = R_;
    if (Rpoint > Rmax)
    {
        return false;
    }
    else // more expensive check, winding number
    {
        real_number ngrid = 1000;
        real_number dtheta = 2.0 * M_PI / ngrid;

        real_number accumulated_angle = 0.0;

        for (int i = 0; i < ngrid; i++)
        {
            real_number theta = i * dtheta;
            Point<DIM, real_number> Pcurve = ParametricCoords(theta);
            Point<DIM, real_number> Pcurve_next = ParametricCoords(theta + dtheta);

            // angle = atan2(ycurvenext-ypoint,xcurvenext-xpoint) - atan2(ycurve-ypoint,xcurve-xpoint)
            real_number angle = atan2(Pcurve_next.get(1) - Paux.get(1), Pcurve_next.get(0) - Paux.get(0)) - atan2(Pcurve.get(1) - Paux.get(1), Pcurve.get(0) - Paux.get(0));
            // normalize to -pi,pi
            if (angle > M_PI)
            {
                angle -= 2.0 * M_PI;
            }
            else if (angle < -M_PI)
            {
                angle += 2.0 * M_PI;
            }
            accumulated_angle += angle;
        }
        real_number tol = 0.1 * 2.0 * M_PI;

        if (abs(accumulated_angle) > tol)
        {
            return true;
        }
        else
        {
            return false;
        }
    }
}

bool HipoCycloid::isInsidePlusTol(Point<DIM, real_number> P)
{
    return false;
}

Point<DIM, real_number> HipoCycloid::ParametricCoords(real_number theta)
{
    real_number x = r_ * (k_ - 1) * cos(theta) + r_ * cos((k_ - 1) * theta);
    real_number y = r_ * (k_ - 1) * sin(theta) - r_ * sin((k_ - 1) * theta);
    return Point<DIM, real_number>{x, y};
}

Point<DIM, real_number> HipoCycloid::ParametricCoordsDerivative(real_number theta)
{
    real_number dx = -r_ * (k_ - 1) * sin(theta) - r_ * (k_ - 1) * sin((k_ - 1) * theta);
    real_number dy = r_ * (k_ - 1) * cos(theta) - r_ * (k_ - 1) * cos((k_ - 1) * theta);
    return Point<DIM, real_number>{dx, dy};
}

bool IsInsideAll(const std::vector<Obstacle *> &obstacles, const Point<DIM, real_number> &point)
{
    for (const auto &obstacle : obstacles)
    {
        if (obstacle->isInside(point))
        {
            return true;
        }
    }
    return false; // If not inside any obstacle
}