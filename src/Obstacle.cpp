#include "Obstacle.hpp"

void AddFlatWallNewBC(particles &vd,
                      const int k0,
                      const int kmax,
                      const Point<DIM, double> Corner,
                      const Point<DIM, double> UnitOffset,
                      const double dx,
                      const Point<DIM, double> obstacle_centre,
                      const Point<DIM, double> obstacle_velocity,
                      const Parameters &arg_p,
                      const double obstacle_omega)
{

    for (int k = k0; k < kmax; k++)
    {
        Point<DIM, double> WallPos = Corner + k * UnitOffset;

        vd.add();
        vd.template getLastProp<type>() = BOUNDARY;
        for (int xyz = 0; xyz < DIM; xyz++)
        {
            vd.getLastPos()[xyz] = WallPos.get(xyz); //+ ((double)rand() / RAND_MAX - 0.5) * dp;
            vd.template getLastProp<velocity>()[xyz] = obstacle_velocity.get(xyz);
            vd.template getLastProp<force>()[xyz] = 0.0;
            vd.template getLastProp<force_transport>()[xyz] = obstacle_centre.get(xyz);
            vd.template getLastProp<v_transport>()[xyz] = 0.0;
            vd.template getLastProp<normal_vector>()[xyz] = 0.0;
        }

        vd.template getLastProp<pressure>() = 0.0;
        vd.template getLastProp<rho>() = arg_p.rho_zero;
        vd.template getLastProp<drho>() = 0.0;
        vd.template getLastProp<curvature_boundary>() = 0.0;
        vd.template getLastProp<arc_length>() = dx;
        vd.template getLastProp<vd_omega>() = obstacle_omega;
    }
}

void AddFlatWallModNewBC(particles &vd,
                         const int k0,
                         const int kmax,
                         const Point<DIM, double> Corner,
                         const Point<DIM, double> UnitOffset,
                         const double dx,
                         const Point<DIM, double> obstacle_centre,
                         const Point<DIM, double> obstacle_velocity,
                         const Point<DIM, double> given_normal,
                         const double obstacle_omega)
{

    for (int k = k0; k < kmax; k++)
    {
        Point<DIM, double> WallPos = Corner + (k + 0.5) * UnitOffset;

        vd.add();
        vd.template getLastProp<type>() = BOUNDARY;
        for (int xyz = 0; xyz < DIM; xyz++)
        {
            vd.getLastPos()[xyz] = WallPos.get(xyz); //+ ((double)rand() / RAND_MAX - 0.5) * dp;
            vd.template getLastProp<velocity>()[xyz] = obstacle_velocity.get(xyz);
            vd.template getLastProp<force>()[xyz] = 0.0;
            vd.template getLastProp<force_transport>()[xyz] = obstacle_centre.get(xyz);
            vd.template getLastProp<v_transport>()[xyz] = 0.0;
            vd.template getLastProp<normal_vector>()[xyz] = given_normal.get(xyz);
        }

        vd.template getLastProp<pressure>() = 0.0;
        vd.template getLastProp<rho>() = 0.0;
        vd.template getLastProp<drho>() = 0.0;
        vd.template getLastProp<curvature_boundary>() = 0.0;
        vd.template getLastProp<arc_length>() = dx;
        vd.template getLastProp<vd_omega>() = obstacle_omega;
    }
}

bool isAvobeLine(Point<DIM, double> P, Point<DIM, double> Q, Point<DIM, double> EvalPoint, double dp)
{
    // P,Q two points forming a line, EvalPoint the point to check if it is above the line or not
    // with a small epsilon to avoid particles to be very close to the line
    double slope = (Q.get(1) - P.get(1)) / (Q.get(0) - P.get(0));
    // we compare y- y1 + m(x-x1)
    // =0 means the point is in the line
    // >0 means the point is above the line
    // <0 means the point is below the line
    double yline = P.get(1) + slope * (EvalPoint.get(0) - P.get(0));
    double epsilon = 0.10 * dp;
    if (EvalPoint.get(1) > yline + epsilon)
        return true;
    else
        return false;
}

bool isBelowLine(Point<DIM, double> P, Point<DIM, double> Q, Point<DIM, double> EvalPoint, double dp)
{
    // P,Q two points forming a line, EvalPoint the point to check if it is above the line or not
    // with a small epsilon to avoid particles to be very close to the line
    double slope = (Q.get(1) - P.get(1)) / (Q.get(0) - P.get(0));
    // we compare y- y1 + m(x-x1)
    // =0 means the point is in the line
    // >0 means the point is above the line
    // <0 means the point is below the line
    double yline = P.get(1) + slope * (EvalPoint.get(0) - P.get(0));
    double epsilon = 0.10 * dp;
    if (EvalPoint.get(1) < yline - epsilon)
        return true;
    else
        return false;
}

// Obstacle class
Obstacle::Obstacle(Point<DIM, double> Centre,
                   const Parameters &p,
                   Point<DIM, double> vel,
                   double omega,
                   double rf) : Centre_(Centre), params_(p), LinearVelocity_(vel), AngularVelocity_(omega), refine_factor(rf) {}
Obstacle::Obstacle(const Parameters &p) : params_(p) {}

// EmptyObstacle class
EmptyObstacle::EmptyObstacle(const Parameters &p) : Obstacle(p) {}
bool EmptyObstacle::isInside(Point<DIM, double> P)
{
    return false;
}
void EmptyObstacle::AddObstacle(particles &vd)
{
}

// CylinderObstacle class
CylinderObstacle::CylinderObstacle(double Radius,
                                   Point<DIM, double> centre,
                                   const Parameters &p,
                                   Point<DIM, double> vel,
                                   double omega,
                                   double rf) : Obstacle(centre, p, vel, omega, rf), Radius_(Radius), Cylinder_(centre, Radius) {}

bool CylinderObstacle::isInside(Point<DIM, double> P)
{
    double radius_aux = Radius_;
    // if (params_.BC_TYPE == NEW_NO_SLIP)
    //     radius_aux = Radius_ + 0.1 * params_.dp;
    Sphere<DIM, double> Cylinderaux(Centre_, radius_aux);
    return Cylinderaux.isInside(P);
}
bool CylinderObstacle::isOutside(Point<DIM, double> P) // for outer cylinder in taylor couette
{
    double radius_aux = Radius_;
    // if (params_.BC_TYPE == NEW_NO_SLIP)
    //     radius_aux = Radius_ - 0.1 * params_.dp;
    Sphere<DIM, double> Cylinderaux(Centre_, radius_aux);
    return Cylinderaux.isInside(P);
}

void CylinderObstacle::AddObstacle(particles &vd)
{
    const double dx_self = params_.dp / refine_factor;
    const double perimeter = 2.0 * M_PI * Radius_;
    const int Np_cylinder = ceil(perimeter / dx_self);
    const double dtheta = 2.0 * M_PI / (double)Np_cylinder;
    const double dxwall = dtheta * Radius_;
    double theta = 0.0;

    Point<DIM, double> Cylinder_particle;

    for (int k = 0; k < Np_cylinder; k++)
    {

        Cylinder_particle[0] = Centre_.get(0) + Radius_ * cos(theta);
        Cylinder_particle[1] = Centre_.get(1) + Radius_ * sin(theta);

        if (DIM == 3)
        {
            Cylinder_particle[2] = Centre_.get(2);
        }

        vd.add();
        vd.template getLastProp<type>() = BOUNDARY;
        for (int xyz = 0; xyz < DIM; xyz++)
        {
            vd.getLastPos()[xyz] = Cylinder_particle.get(xyz); // + ((double)rand() / RAND_MAX - 0.5) * dp;
            vd.template getLastProp<force>()[xyz] = 0.0;
            vd.template getLastProp<force_transport>()[xyz] = Centre_.get(xyz);
            vd.template getLastProp<v_transport>()[xyz] = 0.0;
            vd.template getLastProp<normal_vector>()[xyz] = 0.0;
            vd.template getLastProp<velocity>()[xyz] = LinearVelocity_.get(xyz);
        }

        vd.template getLastProp<pressure>() = 0.0;
        vd.template getLastProp<rho>() = params_.rho_zero;
        vd.template getLastProp<drho>() = 0.0;
        vd.template getLastProp<curvature_boundary>() = 0.0; // 1.0 / Radius_;
        vd.template getLastProp<arc_length>() = dxwall;
        vd.template getLastProp<vd_omega>() = AngularVelocity_;
        theta += dtheta;
    }
}

// Ellipse class
EllipticObstacle::EllipticObstacle(double Major,
                                   double Minor,
                                   double tilt,
                                   Point<DIM, double> centre,
                                   const Parameters &p,
                                   Point<DIM, double> vel,
                                   double omega,
                                   double rf) : Obstacle(centre, p, vel, omega, rf), Major_(Major), Minor_(Minor), tilt_(-tilt * M_PI / 180.0) {}

bool EllipticObstacle::isInside(Point<DIM, double> P)
{
    // get point minus the centre
    Point<DIM, double> Paux = P - Centre_;
    // rotate by -tilt
    ApplyRotation(Paux, -tilt_, {0.0, 0.0});
    // check if the point is inside the ellipse
    return (Paux.get(0) * Paux.get(0) / (Major_ * Major_) + Paux.get(1) * Paux.get(1) / (Minor_ * Minor_) <= 1.0);
}

void EllipticObstacle::AddObstacle(particles &vd)
{
    // first we need to determine the perimeter of the ellipse
    // we can use c++ std::comp_ellint_2 function
    // first compute eccentricity
    const double e = sqrt(1.0 - Minor_ * Minor_ / Major_ / Major_);
    // compute the perimeter
    const double perimeter = 4.0 * Major_ * std::comp_ellint_2(e);
    const double dx_self = params_.dp / refine_factor;
    const int Np_ellipse_ = ceil(perimeter / dx_self);
    const double dxwall = perimeter / (double)Np_ellipse_;

    const int M = 1e6; // number of points to discretize the ellipse
    double accumulated_arc = 0.0;
    double theta = 0.0;
    double dtheta = 2.0 * M_PI / (double)M;

    Point<DIM, double> P = parametricEllipse(theta);

    // add the first particle
    vd.add();
    vd.template getLastProp<type>() = BOUNDARY;
    for (int xyz = 0; xyz < DIM; xyz++)
    {
        vd.getLastPos()[xyz] = P.get(xyz);
        vd.template getLastProp<force>()[xyz] = 0.0;
        vd.template getLastProp<force_transport>()[xyz] = Centre_.get(xyz);
        vd.template getLastProp<v_transport>()[xyz] = 0.0;
        vd.template getLastProp<normal_vector>()[xyz] = 0.0;
        vd.template getLastProp<velocity>()[xyz] = LinearVelocity_.get(xyz);
    }

    vd.template getLastProp<pressure>() = 0.0;
    vd.template getLastProp<rho>() = params_.rho_zero;
    vd.template getLastProp<drho>() = 0.0;
    vd.template getLastProp<curvature_boundary>() = 0.0;
    vd.template getLastProp<arc_length>() = dxwall;
    vd.template getLastProp<vd_omega>() = AngularVelocity_;

    Point<DIM, double> P_prev = P;

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
            vd.template getLastProp<type>() = BOUNDARY;
            for (int xyz = 0; xyz < DIM; xyz++)
            {
                vd.getLastPos()[xyz] = P.get(xyz);
                vd.template getLastProp<force>()[xyz] = 0.0;
                vd.template getLastProp<force_transport>()[xyz] = Centre_.get(xyz);
                vd.template getLastProp<v_transport>()[xyz] = 0.0;
                vd.template getLastProp<normal_vector>()[xyz] = 0.0;
                vd.template getLastProp<velocity>()[xyz] = LinearVelocity_.get(xyz);
            }

            vd.template getLastProp<pressure>() = 0.0;
            vd.template getLastProp<rho>() = params_.rho_zero;
            vd.template getLastProp<drho>() = 0.0;
            vd.template getLastProp<curvature_boundary>() = 0.0;
            vd.template getLastProp<arc_length>() = dxwall;
            vd.template getLastProp<vd_omega>() = AngularVelocity_;

            // reset accumulated arc
            accumulated_arc = 0.0;
        }
    }
}

Point<DIM, double> EllipticObstacle::parametricEllipse(double theta)
{
    // given the angle theta, return the point in the ellipse
    Point<DIM, double> P;
    P.get(0) = Centre_.get(0) + (Major_ * cos(theta) * cos(tilt_) - Minor_ * sin(theta) * sin(tilt_));
    P.get(1) = Centre_.get(1) + (Major_ * cos(theta) * sin(tilt_) + Minor_ * sin(theta) * cos(tilt_));
    return P;
}

// RectangleObstacle class
RectangleObstacle::RectangleObstacle(Point<DIM, double> centre,
                                     const Parameters &p,
                                     double BaseLength,
                                     double HeigthLength,
                                     Point<DIM, double> vel,
                                     double omega,
                                     double rf) : Obstacle(centre, p, vel, omega, rf), BaseLength_(BaseLength), HeigthLength_(HeigthLength)
{
    LowerLeft_ = Point<DIM, double>{Centre_.get(0) - BaseLength_ / 2.0,
                                    Centre_.get(1) - HeigthLength_ / 2.0};

    UpperRight_ = Point<DIM, double>{Centre_.get(0) + BaseLength_ / 2.0,
                                     Centre_.get(1) + HeigthLength_ / 2.0};

    // LowerLeft_ = Point<DIM, double>{Centre_.get(0) - (((double)BaseLength_ - 1.0) * params_.dp) / 2.0,
    //                                 Centre_.get(1) - (((double)HeigthLength_ - 1.0) * params_.dp) / 2.0};
    // LowerRight_ = Point<DIM, double>{Centre_.get(0) + (((double)BaseLength_ - 1.0) * params_.dp) / 2.0,
    //                                  Centre_.get(1) - (((double)HeigthLength_ - 1.0) * params_.dp) / 2.0};
    // UpperLeft_ = Point<DIM, double>{Centre_.get(0) - (((double)BaseLength_ - 1.0) * params_.dp) / 2.0,
    //                                 Centre_.get(1) + (((double)HeigthLength_ - 1.0) * params_.dp) / 2.0};
    // UpperRight_ = Point<DIM, double>{Centre_.get(0) + (((double)BaseLength_ - 1.0) * params_.dp) / 2.0,
    //                                  Centre_.get(1) + (((double)HeigthLength_ - 1.0) * params_.dp) / 2.0};

    Rectangle_ = Box<DIM, double>(LowerLeft_, UpperRight_);
}

bool RectangleObstacle::isInside(Point<DIM, double> P)
{
    return Rectangle_.isInside(P);
}

void RectangleObstacle::AddObstacle(particles &vd)
{

    const double dx_self = params_.dp / refine_factor;

    // Horizontal walls
    const int N_bottom = ceil(BaseLength_ / dx_self);
    const double dxwall_bottom = BaseLength_ / N_bottom;
    Point<DIM, double> Xoffset = {dxwall_bottom, 0.0};

    // Lower wall
    AddFlatWallNewBC(vd, 0, N_bottom + 1, LowerLeft_, Xoffset, dxwall_bottom, Centre_, LinearVelocity_, params_, AngularVelocity_);
    // Upper wall
    AddFlatWallNewBC(vd, 0, N_bottom + 1, UpperRight_, -1.0 * Xoffset, dxwall_bottom, Centre_, LinearVelocity_, params_, AngularVelocity_);

    // Vertical walls
    const int N_right = ceil(HeigthLength_ / dx_self);
    const double dxwall_right = HeigthLength_ / N_right;
    Point<DIM, double> Yoffset = {0.0, dxwall_right};

    // Left wall
    AddFlatWallNewBC(vd, 1, N_right, LowerLeft_, Yoffset, dxwall_right, Centre_, LinearVelocity_, params_, AngularVelocity_);
    // Right wall
    AddFlatWallNewBC(vd, 1, N_right, UpperRight_, -1.0 * Yoffset, dxwall_right, Centre_, LinearVelocity_, params_, AngularVelocity_);

    // WHEN using non mod flat wall this are the k and kmax values
    // // Lower wall
    // AddFlatWallModNewBC(vd, 0, N_bottom + 1, LowerLeft_, Xoffset, dxwall_bottom, Centre_, LinearVelocity_, {0.0, -1.0}, AngularVelocity_);
    // // Upper wall
    // AddFlatWallModNewBC(vd, 0, N_bottom + 1, UpperRight_, -1.0 * Xoffset, dxwall_bottom, Centre_, LinearVelocity_, {0.0, 1.0}, AngularVelocity_);

    // // Vertical walls
    // const int N_right = ceil(HeigthLength_ / dx_self);
    // const double dxwall_right = HeigthLength_ / N_right;
    // Point<DIM, double> Yoffset = {0.0, dxwall_right};

    // // Left wall
    // AddFlatWallModNewBC(vd, 1, N_right, LowerLeft_, Yoffset, dxwall_right, Centre_, LinearVelocity_, {-1.0, 0.0}, AngularVelocity_);
    // // Right wall
    // AddFlatWallModNewBC(vd, 1, N_right, UpperRight_, -1.0 * Yoffset, dxwall_right, Centre_, LinearVelocity_, {1.0, 0.0}, AngularVelocity_);
}

// TriangleObstacle class
TriangleObstacle::TriangleObstacle(Point<DIM, double> centre,
                                   const Parameters &p,
                                   double BaseLength,
                                   double HeigthLength,
                                   Point<DIM, double> vel,
                                   double omega,
                                   double rf) : Obstacle(centre, p, vel, omega, rf), BaseLength_(BaseLength), HeigthLength_(HeigthLength)
{
    LowerLeft_ = Point<DIM, double>{Centre_.get(0) - 2.0 * BaseLength_ / 3.0,
                                    Centre_.get(1) - HeigthLength_ / 3.0};

    LowerRight_ = Point<DIM, double>{Centre_.get(0) + BaseLength_ / 3.0,
                                     Centre_.get(1) - HeigthLength_ / 3.0};

    UpperRight_ = Point<DIM, double>{Centre_.get(0) + BaseLength_ / 3.0,
                                     Centre_.get(1) + 2.0 * HeigthLength_ / 3.0};

    ContainingRectangle_ = Box<DIM, double>(LowerLeft_, UpperRight_);
}

bool TriangleObstacle::isInside(Point<DIM, double> P)
{
    return (ContainingRectangle_.isInside(P) && !isAvobeLine(LowerLeft_, UpperRight_, P, params_.dp));
}

void TriangleObstacle::AddObstacle(particles &vd)
{
    const double dx_self = params_.dp / refine_factor;
    // Lower wall
    const int N_bottom = ceil(BaseLength_ / dx_self);
    const double dxwall_bottom = BaseLength_ / N_bottom;
    const Point<DIM, double> Xoffset = {dxwall_bottom, 0.0};
    AddFlatWallNewBC(vd, 0, N_bottom + 1, LowerLeft_, Xoffset, dxwall_bottom, Centre_, LinearVelocity_, params_, AngularVelocity_);

    // Right wall
    const int N_right = ceil(HeigthLength_ / dx_self);
    const double dxwall_right = HeigthLength_ / N_right;
    const Point<DIM, double> Yoffset = {0.0, dxwall_right};
    AddFlatWallNewBC(vd, 1, N_right + 1, LowerRight_, Yoffset, dxwall_right, Centre_, LinearVelocity_, params_, AngularVelocity_);

    //  Hypothenuse wall
    // We want particles spaced roughly by dp
    const double HypothenuseLength = sqrt(BaseLength_ * BaseLength_ + HeigthLength_ * HeigthLength_);
    const int Ndiag = ceil(HypothenuseLength / dx_self);  // integer number of particles that can fit in the diagonal
    const double dxwall_diag = HypothenuseLength / Ndiag; // actual spacing between particles ( close to dp but not exactly)
    const double sin_theta = HeigthLength_ / HypothenuseLength;
    const double cos_theta = BaseLength_ / HypothenuseLength;
    const Point<DIM, double> Diagoffset{dxwall_diag * cos_theta, dxwall_diag * sin_theta};

    AddFlatWallNewBC(vd, 1, Ndiag, UpperRight_, -1.0 * Diagoffset, dxwall_diag, Centre_, LinearVelocity_, params_, AngularVelocity_);
}

TriangleEqui::TriangleEqui(Point<DIM, double> centre,
                           const Parameters &p,
                           double sidelength,
                           Point<DIM, double> vel,
                           double omega,
                           double rf) : Obstacle(centre, p, vel, omega, rf), SideLength_(sidelength)
{

    UpperRight_ = Point<DIM, double>{Centre_.get(0) + sqrt(3) * SideLength_ / 6.0,
                                     Centre_.get(1) + SideLength_ / 2.0};
    LowerRight_ = Point<DIM, double>{Centre_.get(0) + sqrt(3) * SideLength_ / 6.0,
                                     Centre_.get(1) - SideLength_ / 2.0};
    TriangleTip_ = Point<DIM, double>{Centre_.get(0) - sqrt(3.0) * SideLength_ / 3.0,
                                      Centre_.get(1)};

    LowerLeft_ = Point<DIM, double>{TriangleTip_.get(0),
                                    LowerRight_.get(1)};

    ContainingRectangle_ = Box<DIM, double>(LowerLeft_, UpperRight_);
}

bool TriangleEqui::isInside(Point<DIM, double> P)
{
    return (ContainingRectangle_.isInside(P) && !isAvobeLine(TriangleTip_, UpperRight_, P, params_.dp) && !isBelowLine(TriangleTip_, LowerRight_, P, params_.dp));
}

void TriangleEqui::AddObstacle(particles &vd)
{
    double dx_self = params_.dp / refine_factor;
    const int N_wall = ceil(SideLength_ / dx_self);
    const double dxwall = SideLength_ / N_wall;
    Point<DIM, double> Yoffset = {0.0, dxwall};

    // Right wall
    AddFlatWallNewBC(vd, 0, N_wall + 1, LowerRight_, Yoffset, dxwall, Centre_, LinearVelocity_, params_, AngularVelocity_);

    // Inclined walls

    const double cos_theta = sqrt(3.0) / 2.0;
    const double sin_theta = 1 / 2.0;
    const Point<DIM, double> Diagoffset{dxwall * cos_theta, dxwall * sin_theta};
    Point<DIM, double> DiagoffsetNW = Diagoffset;
    DiagoffsetNW.get(0) = -1.0 * DiagoffsetNW.get(0);
    //  Hypothenuse upper wall
    AddFlatWallNewBC(vd, 1, N_wall + 1, UpperRight_, -1.0 * Diagoffset, dxwall, Centre_, LinearVelocity_, params_, AngularVelocity_);
    // Hypothenuse lower wall
    AddFlatWallNewBC(vd, 1, N_wall, LowerRight_, DiagoffsetNW, dxwall, Centre_, LinearVelocity_, params_, AngularVelocity_);
}
