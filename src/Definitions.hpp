#ifndef DEFINITIONS_H
#define DEFINITIONS_H

#include "Vector/vector_dist.hpp"
#include "Draw/DrawParticles.hpp"

// Number of dimensions, changes size of arrays,
// changes the normalizations of the kernels
#define DIM 2

// Type of particles
#define BOUNDARY 0
#define FLUID 1

// Type of boundary condition
#define NO_SLIP 0
#define NEW_NO_SLIP 1

// Type of scenario
#define POISEUILLE 0
#define COUETTE 1
#define HYDROSTATIC 2
#define CYLINDER_ARRAY 3
#define CYLINDER_LATTICE 4
#define SQUARE 5
#define TRIANGLE 6
#define TRIANGLE_EQUILATERAL 7
#define CAVITY 8
#define STEP 9
#define TAYLOR_COUETTE 10
#define MOVING_OBSTACLE 11

// Type of density calculation
#define DENSITY_SUMMATION 0
#define DENSITY_DIFFERENTIAL 1

// Alias of the vector_dist fields
// FLUID or BOUNDARY
const size_t type = 0;
// Density
const int rho = 1;
// pressure
const int pressure = 2;
// Delta rho ( for differential density )
const int drho = 3;
// calculated force
const int force = 4;
// velocity
const int velocity = 5;
// background pressure force
const int force_transport = 6;
// transport velocity
const int v_transport = 7;
// normal vector
const int normal_vector = 8;
// curvature
const int curvature_boundary = 9;
// arc length wall
const int arc_length = 10;
// volume
const int vd_volume = 11;
// angular velocity
const int vd_omega = 12;

typedef vector_dist<DIM, double, aggregate<size_t, double, double, double, double[DIM], double[DIM], double[DIM], double[DIM], double[DIM], double, double, double[3], double>> particles;
//                                          |         |     |        |        |            |           |		     |           |	        |       |          | 		   |
//                                          |         |     |        |        |            |           |		     |	     	 |	        |       |	       |	       |
//                                       type        rho  pressure delta   force          vel         F_t           vel_t    normal curvature   arc_length     vol       omega
//                                                                 density

typedef vector_dist<DIM, double, aggregate<double>> probe_particles;

class Parameters
{
public:
    double dp;

    int SCENARIO;
    int BC_TYPE;
    int DENSITY_TYPE;
    int WRITER;             // VTK_WRITER or CSV_WRITER
    int PROBES_ENABLED = 0; // 0 for disabled, 1 for enabled

    //////// DECLARATION OF GLOBAL PARAMETERS /////////////////////////////////////////////////////////////
    // Output file name
    std::string filename;
    std::vector<std::string> probe_filenames;
    // Physical size of the fluid domain, it goes from (0,0,0) to (length[0],length[1],length[2])
    // First fluid particle will be placed at (dp/2,dp/2,dp/2) and the last particle will be placed at (length[0]-dp/2,length[1]-dp/2,length[2]-dp/2)
    double length[DIM];

    // Boundary conditions
    size_t bc[DIM];
    // Number of boundary particles in each direction
    size_t Nboundary[DIM];
    // Number of fluid particles in each direction
    size_t Nfluid[DIM];

    // problem specific length scale to compute the Reynolds number
    double LengthScale;
    double rf;
    // Factor relating H (smoothing length) and dp (particle spacing)
    const double Hconst = 1.0;
    // Smoothing length
    double H;
    // Radius of the kernel support
    double r_threshold;
    // Normalization constant for the kernels
    // double Kcubic;
    double Kquintic;
    // Reynolds number
    double Re;
    // maximum velocity
    double umax;
    // Reference density
    double rho_zero;
    // Gamma in eq of state
    double gamma_;
    // Constant used for the sound speed, number of times the max velocity
    const double coeff_sound = 10.0;
    // Sound speed
    double cbar;
    // Eq of state constant ( p0 )
    double B;
    // background pressure in eq of state
    const double xi = 0.0;
    // Gravity vector and magnitude
    Point<DIM, double> gravity_vector = {0.0};
    double gravity = 0.0;
    // wall velocity
    Point<DIM, double> vw_top = {0.0};
    Point<DIM, double> vw_bottom = {0.0};
    // Mass of the fluid and boundary particles
    double MassFluid;
    double MassBound;
    // Kinematic viscosity
    double nu;
    // Dynamic viscosity
    double eta;
    // Factor relating Pbackground and B, Pbackground = Bfactor*B
    double Bfactor;
    // Background pressure in the transport force
    double Pbackground;
    // End simulation time
    double t_end;
    // Constant used to define time integration
    const double CFLnumber = 0.1;
    // Controls otput file frequency, 1 means 1 file per time unit, 10 means 10 files per time unit, etc.
    double write_const;
    // iteration counter
    size_t cnt = 0;

    Point<DIM, double> ObstacleCenter;
    Point<DIM, double> ObstacleVelocity;
    double ObstacleOmega;

    // Taylor Couette
    double Rin;
    double Rout;
    double Win;
    double Wout;

    // Obstacle
    double ObstacleBase;
    double ObstacleHeight;
    // // custom string
    // std::string custom_string = "";
};

#endif // DEFINITIONS_H
