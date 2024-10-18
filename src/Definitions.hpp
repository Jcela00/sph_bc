#ifndef DEFINITIONS_H
#define DEFINITIONS_H

#include "Vector/vector_dist.hpp"
#include "Draw/DrawParticles.hpp"

// Type of floating point number
typedef float real_number;

// Number of dimensions, changes size of arrays, changes the normalizations of the kernels
#define DIM 2

// Type of particles
#define BOUNDARY 0
#define OBSTACLE 1
#define FLUID 2

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
#define ELLIPSE 12
#define DAM_BREAK 13
#define TRIANGLE_TEST 14

// Type of density calculation
#define DENSITY_SUMMATION 0
#define DENSITY_DIFFERENTIAL 1

// Alias of the vector_dist fields
const int vd0_type = 0; // can be BOUNDARY, OBSTACLE, FLUID
// Density
const int vd1_rho = 1;
// pressure
const int vd2_pressure = 2;
// Delta rho ( for differential density )
const int vd3_drho = 3;
// velocity
const int vd4_velocity = 4;
// transport velocity
const int vd5_velocity_t = 5;
// calculated force
const int vd6_force = 6;
// background pressure force
const int vd7_force_t = 7;
// normal vector
const int vd8_normal = 8;
// volume
const int vd9_volume = 9;
// angular velocity
const int vd10_omega = 10;
// vorticity
const int vd11_vorticity = 11;
// reduction
const int vd12_vel_red = 12;
// force_red
const int vd13_force_red_x = 13;
const int vd14_force_red_y = 14;

typedef vector_dist_gpu<DIM, real_number, aggregate<unsigned int, real_number, real_number, real_number, real_number[DIM], real_number[DIM], real_number[DIM], real_number[DIM], real_number[DIM], real_number[3], real_number, real_number, real_number, real_number, real_number>> particles;
//                                                      |       |            |            |            |                 |                 |                 |                 |                 |                 |            |           |           |            |
//                                                      |       |            |            |            |                 |                 |                 |                 |                 |                 |            |           |           |            |
//                                                     type    rho        pressure       drho       velocity       v_transport           force        force_transport        normal_vector      volume            omega        vorticity   red       force_redx    force_redy

typedef vector_dist_gpu<DIM, real_number, aggregate<real_number>> probe_particles;

// Class to store parameters that we can not (or dont want) to
// copy into the GPU __constant__ memory
class AuxiliarParameters
{
public:
    std::string filename;
    std::vector<std::string> probe_filenames;

    // iteration counter
    size_t cnt = 0;
};

class Parameters
{
public:
    real_number dp;

    int SCENARIO;
    int BC_TYPE;
    int DENSITY_TYPE;
    int WRITER;         // VTK_WRITER or CSV_WRITER
    int PROBES_ENABLED; // 0 for disabled, 1 for enabled

    // Physical size of the fluid domain, it goes from (0,0,0)
    // to (length[0],length[1],length[2])
    // First fluid particle will be placed at (dp/2,dp/2,dp/2) and the
    // last particle will be placed at (length[0]-dp/2,length[1]-dp/2,length[2]-dp/2)
    real_number length[DIM];

    // Boundary conditions
    size_t bc[DIM];

    // Number of boundary particles in each direction
    size_t Nboundary[DIM];

    // Number of fluid particles in each direction
    size_t Nfluid[DIM];

    // LengthScale (problem dependent)
    real_number LengthScale;

    // Refinement factor of boundary in new bc
    real_number rf;

    // Factor relating H (smoothing length) and dp (particle spacing)
    real_number Hconst; // = 1.0;

    // Smoothing length
    real_number H;

    // Radius of the kernel support
    real_number r_cut;
    real_number r_cut2;

    // Normalization constant for the kernels
    real_number Kquintic;

    // Reynolds number
    real_number Re;

    // maximum velocity
    real_number umax;

    // Reference density
    real_number rho0;

    // Gamma in eq of state
    real_number gamma;

    // Constant used for the sound speed, number of times the max velocity
    real_number coeff_sound; // = 10.0;

    // Sound speed
    real_number cbar;

    // Eq of state constant ( p0 )
    real_number B;

    // background pressure in eq of state
    real_number xi; //= 0.0;

    // Gravity vector and magnitude
    real_number gravity_vector[DIM];
    real_number gravity;

    // wall velocity
    real_number vw_top[DIM];
    real_number vw_bottom[DIM];

    // Mass of the fluid and boundary particles
    real_number MassFluid;
    real_number MassBound;

    // Kinematic viscosity
    real_number nu;

    // Dynamic viscosity
    real_number eta;

    // Factor relating Pbackground and B, Pbackground = Bfactor*B
    real_number Bfactor;

    // Background pressure in the transport force
    real_number Pbackground;

    // End simulation time
    real_number t_end;

    // Constant used to define time integration
    real_number CFLnumber; // = 0.1;

    // Controls otput file frequency, 1 means 1 file per time unit, 10 means 10 files per time unit, etc.
    real_number write_const;

    // Obstacle properties
    real_number ObstacleCenter[DIM];
    real_number ObstacleVelocity[DIM];
    real_number ObstacleOmega;

    // Obstacle
    real_number ObstacleBase;
    real_number ObstacleHeight;
    real_number ObstacleTilt;

    // Taylor Couette specific parameters
    real_number Rin;
    real_number Rout;
    real_number Win;
    real_number Wout;

    // Dam break specific parameters
    real_number waterH;
    real_number waterB;

    // // constructor
    // Parameters()
    // {
    // }

    // void copy_params_host_to_device(Parameters &DeviceParameters, const Parameters HostParams)
    // {
    //     // only to be called for the device parameters
    //     cudaMemcpyToSymbol(DeviceParameters.dp, &(HostParams.dp), sizeof(real_number));
    //     cudaMemcpyToSymbol(DeviceParameters.SCENARIO, &(HostParams.SCENARIO), sizeof(int));
    //     cudaMemcpyToSymbol(DeviceParameters.BC_TYPE, &(HostParams.BC_TYPE), sizeof(int));
    //     cudaMemcpyToSymbol(DeviceParameters.DENSITY_TYPE, &(HostParams.DENSITY_TYPE), sizeof(int));
    //     cudaMemcpyToSymbol(DeviceParameters.WRITER, &(HostParams.WRITER), sizeof(int));
    //     cudaMemcpyToSymbol(DeviceParameters.PROBES_ENABLED, &(HostParams.PROBES_ENABLED), sizeof(int));
    //     cudaMemcpyToSymbol(DeviceParameters.length, &(HostParams.length), sizeof(real_number) * DIM);
    //     cudaMemcpyToSymbol(DeviceParameters.bc, &(HostParams.bc), sizeof(size_t) * DIM);
    //     cudaMemcpyToSymbol(DeviceParameters.Nboundary, &(HostParams.Nboundary), sizeof(size_t) * DIM);
    //     cudaMemcpyToSymbol(DeviceParameters.Nfluid, &(HostParams.Nfluid), sizeof(size_t) * DIM);
    //     cudaMemcpyToSymbol(DeviceParameters.LengthScale, &(HostParams.LengthScale), sizeof(real_number));
    //     cudaMemcpyToSymbol(DeviceParameters.rf, &(HostParams.rf), sizeof(real_number));
    //     cudaMemcpyToSymbol(DeviceParameters.Hconst, &(HostParams.Hconst), sizeof(real_number));
    //     cudaMemcpyToSymbol(DeviceParameters.H, &(HostParams.H), sizeof(real_number));
    //     cudaMemcpyToSymbol(DeviceParameters.r_cut, &(HostParams.r_cut), sizeof(real_number));
    //     cudaMemcpyToSymbol(DeviceParameters.Kquintic, &(HostParams.Kquintic), sizeof(real_number));
    //     cudaMemcpyToSymbol(DeviceParameters.Re, &(HostParams.Re), sizeof(real_number));
    //     cudaMemcpyToSymbol(DeviceParameters.umax, &(HostParams.umax), sizeof(real_number));
    //     cudaMemcpyToSymbol(DeviceParameters.rho0, &(HostParams.rho0), sizeof(real_number));
    //     cudaMemcpyToSymbol(DeviceParameters.gamma, &(HostParams.gamma), sizeof(real_number));
    //     cudaMemcpyToSymbol(DeviceParameters.coeff_sound, &(HostParams.coeff_sound), sizeof(real_number));
    //     cudaMemcpyToSymbol(DeviceParameters.cbar, &(HostParams.cbar), sizeof(real_number));
    //     cudaMemcpyToSymbol(DeviceParameters.B, &(HostParams.B), sizeof(real_number));
    //     cudaMemcpyToSymbol(DeviceParameters.xi, &(HostParams.xi), sizeof(real_number));
    //     cudaMemcpyToSymbol(DeviceParameters.gravity_vector, &(HostParams.gravity_vector), sizeof(real_number) * DIM);
    //     cudaMemcpyToSymbol(DeviceParameters.gravity, &(HostParams.gravity), sizeof(real_number));
    //     cudaMemcpyToSymbol(DeviceParameters.vw_top, &(HostParams.vw_top), sizeof(real_number) * DIM);
    //     cudaMemcpyToSymbol(DeviceParameters.vw_bottom, &(HostParams.vw_bottom), sizeof(real_number) * DIM);
    //     cudaMemcpyToSymbol(DeviceParameters.MassFluid, &(HostParams.MassFluid), sizeof(real_number));
    //     cudaMemcpyToSymbol(DeviceParameters.MassBound, &(HostParams.MassBound), sizeof(real_number));
    //     cudaMemcpyToSymbol(DeviceParameters.nu, &(HostParams.nu), sizeof(real_number));
    //     cudaMemcpyToSymbol(DeviceParameters.eta, &(HostParams.eta), sizeof(real_number));
    //     cudaMemcpyToSymbol(DeviceParameters.Bfactor, &(HostParams.Bfactor), sizeof(real_number));
    //     cudaMemcpyToSymbol(DeviceParameters.Pbackground, &(HostParams.Pbackground), sizeof(real_number));
    //     cudaMemcpyToSymbol(DeviceParameters.t_end, &(HostParams.t_end), sizeof(real_number));
    //     cudaMemcpyToSymbol(DeviceParameters.CFLnumber, &(HostParams.CFLnumber), sizeof(real_number));
    //     cudaMemcpyToSymbol(DeviceParameters.write_const, &(HostParams.write_const), sizeof(real_number));
    //     cudaMemcpyToSymbol(DeviceParameters.ObstacleCenter, &(HostParams.ObstacleCenter), sizeof(real_number) * DIM);
    //     cudaMemcpyToSymbol(DeviceParameters.ObstacleVelocity, &(HostParams.ObstacleVelocity), sizeof(real_number) * DIM);
    //     cudaMemcpyToSymbol(DeviceParameters.ObstacleOmega, &(HostParams.ObstacleOmega), sizeof(real_number));
    //     cudaMemcpyToSymbol(DeviceParameters.Rin, &(HostParams.Rin), sizeof(real_number));
    //     cudaMemcpyToSymbol(DeviceParameters.Rout, &(HostParams.Rout), sizeof(real_number));
    //     cudaMemcpyToSymbol(DeviceParameters.Win, &(HostParams.Win), sizeof(real_number));
    //     cudaMemcpyToSymbol(DeviceParameters.Wout, &(HostParams.Wout), sizeof(real_number));
    //     cudaMemcpyToSymbol(DeviceParameters.ObstacleBase, &(HostParams.ObstacleBase), sizeof(real_number));
    //     cudaMemcpyToSymbol(DeviceParameters.ObstacleHeight, &(HostParams.ObstacleHeight), sizeof(real_number));
    //     cudaMemcpyToSymbol(DeviceParameters.ObstacleTilt, &(HostParams.ObstacleTilt), sizeof(real_number));
    // }

    // void WriteParameters(Parameters params, std::string filename)
    // {
    //     // write all the parameters to a .txt file for debugging and making sure the parameters are read correctly
    //     std::ofstream file;
    //     file.open(filename);
    //     file << "SCENARIO: " << params.SCENARIO << std::endl;
    //     file << "BC_TYPE: " << params.BC_TYPE << std::endl;
    //     file << "DENSITY_TYPE: " << params.DENSITY_TYPE << std::endl;
    //     file << "WRITER: " << params.WRITER << std::endl;
    //     file << "length={" << params.length[0] << "," << params.length[1] << "}" << std::endl;
    //     file << "bc={" << params.bc[0] << "," << params.bc[1] << "}" << std::endl;
    //     file << "Nfluid={" << params.Nfluid[0] << "," << params.Nfluid[1] << "}" << std::endl;
    //     file << "Nboundary={" << params.Nboundary[0] << "," << params.Nboundary[1] << "}" << std::endl;
    //     file << "LengthScale: " << params.LengthScale << std::endl;
    //     file << "rf: " << params.rf << std::endl;
    //     file << "H: " << params.H << std::endl;
    //     file << "r_cut: " << params.r_cut << std::endl;
    //     file << "Kquintic: " << params.Kquintic << std::endl;
    //     file << "Re: " << params.Re << std::endl;
    //     file << "umax: " << params.umax << std::endl;
    //     file << "rho0: " << params.rho0 << std::endl;
    //     file << "gamma: " << params.gamma << std::endl;
    //     file << "coeff_sound: " << params.coeff_sound << std::endl;
    //     file << "cbar: " << params.cbar << std::endl;
    //     file << "B: " << params.B << std::endl;
    //     file << "xi: " << params.xi << std::endl;
    //     file << "gravity_vector: {" << params.gravity_vector[0] << "," << params.gravity_vector[1] << "," << params.gravity_vector[2] << "}" << std::endl;
    //     file << "gravity: " << params.gravity << std::endl;
    //     file << "vw_top: {" << params.vw_top[0] << "," << params.vw_top[1] << "," << params.vw_top[2] << "}" << std::endl;
    //     file << "vw_bottom: {" << params.vw_bottom[0] << "," << params.vw_bottom[1] << "," << params.vw_bottom[2] << "}" << std::endl;
    //     file << "MassFluid: " << params.MassFluid << std::endl;
    //     file << "MassBound: " << params.MassBound << std::endl;
    //     file << "nu: " << params.nu << std::endl;
    //     file << "eta: " << params.eta << std::endl;
    //     file << "Bfactor: " << params.Bfactor << std::endl;
    //     file << "Pbackground: " << params.Pbackground << std::endl;
    //     file << "t_end: " << params.t_end << std::endl;
    //     file << "CFLnumber: " << params.CFLnumber << std::endl;
    //     file << "write_const: " << params.write_const << std::endl;
    //     file << "ObstacleCenter: {" << params.ObstacleCenter[0] << "," << params.ObstacleCenter[1] << "," << params.ObstacleCenter[2] << "}" << std::endl;
    //     file << "ObstacleVelocity: {" << params.ObstacleVelocity[0] << "," << params.ObstacleVelocity[1] << "," << params.ObstacleVelocity[2] << "}" << std::endl;
    //     file << "ObstacleOmega: " << params.ObstacleOmega << std::endl;
    //     file << "Rin: " << params.Rin << std::endl;
    //     file << "Rout: " << params.Rout << std::endl;
    //     file << "Win: " << params.Win << std::endl;
    //     file << "Wout: " << params.Wout << std::endl;
    //     file << "ObstacleBase: " << params.ObstacleBase << std::endl;
    //     file << "ObstacleHeight: " << params.ObstacleHeight << std::endl;
    //     file << "ObstacleTilt: " << params.ObstacleTilt << std::endl;
    // }

    // void copy_params_host_to_device(const Parameters HostParams)
    // {
    //     // only to be called for the device parameters
    //     cudaMemcpyToSymbol(dp, &(HostParams.dp), sizeof(real_number));
    //     cudaMemcpyToSymbol(SCENARIO, &(HostParams.SCENARIO), sizeof(int));
    //     cudaMemcpyToSymbol(BC_TYPE, &(HostParams.BC_TYPE), sizeof(int));
    //     cudaMemcpyToSymbol(DENSITY_TYPE, &(HostParams.DENSITY_TYPE), sizeof(int));
    //     cudaMemcpyToSymbol(WRITER, &(HostParams.WRITER), sizeof(int));
    //     cudaMemcpyToSymbol(PROBES_ENABLED, &(HostParams.PROBES_ENABLED), sizeof(int));
    //     cudaMemcpyToSymbol(length, &(HostParams.length), sizeof(real_number) * DIM);
    //     cudaMemcpyToSymbol(bc, &(HostParams.bc), sizeof(size_t) * DIM);
    //     cudaMemcpyToSymbol(Nboundary, &(HostParams.Nboundary), sizeof(size_t) * DIM);
    //     cudaMemcpyToSymbol(Nfluid, &(HostParams.Nfluid), sizeof(size_t) * DIM);
    //     cudaMemcpyToSymbol(LengthScale, &(HostParams.LengthScale), sizeof(real_number));
    //     cudaMemcpyToSymbol(rf, &(HostParams.rf), sizeof(real_number));
    //     cudaMemcpyToSymbol(Hconst, &(HostParams.Hconst), sizeof(real_number));
    //     cudaMemcpyToSymbol(H, &(HostParams.H), sizeof(real_number));
    //     cudaMemcpyToSymbol(r_cut, &(HostParams.r_cut), sizeof(real_number));
    //     cudaMemcpyToSymbol(Kquintic, &(HostParams.Kquintic), sizeof(real_number));
    //     cudaMemcpyToSymbol(Re, &(HostParams.Re), sizeof(real_number));
    //     cudaMemcpyToSymbol(umax, &(HostParams.umax), sizeof(real_number));
    //     cudaMemcpyToSymbol(rho0, &(HostParams.rho0), sizeof(real_number));
    //     cudaMemcpyToSymbol(gamma, &(HostParams.gamma), sizeof(real_number));
    //     cudaMemcpyToSymbol(coeff_sound, &(HostParams.coeff_sound), sizeof(real_number));
    //     cudaMemcpyToSymbol(cbar, &(HostParams.cbar), sizeof(real_number));
    //     cudaMemcpyToSymbol(B, &(HostParams.B), sizeof(real_number));
    //     cudaMemcpyToSymbol(xi, &(HostParams.xi), sizeof(real_number));
    //     cudaMemcpyToSymbol(gravity_vector, &(HostParams.gravity_vector), sizeof(real_number) * DIM);
    //     cudaMemcpyToSymbol(gravity, &(HostParams.gravity), sizeof(real_number));
    //     cudaMemcpyToSymbol(vw_top, &(HostParams.vw_top), sizeof(real_number) * DIM);
    //     cudaMemcpyToSymbol(vw_bottom, &(HostParams.vw_bottom), sizeof(real_number) * DIM);
    //     cudaMemcpyToSymbol(MassFluid, &(HostParams.MassFluid), sizeof(real_number));
    //     cudaMemcpyToSymbol(MassBound, &(HostParams.MassBound), sizeof(real_number));
    //     cudaMemcpyToSymbol(nu, &(HostParams.nu), sizeof(real_number));
    //     cudaMemcpyToSymbol(eta, &(HostParams.eta), sizeof(real_number));
    //     cudaMemcpyToSymbol(Bfactor, &(HostParams.Bfactor), sizeof(real_number));
    //     cudaMemcpyToSymbol(Pbackground, &(HostParams.Pbackground), sizeof(real_number));
    //     cudaMemcpyToSymbol(t_end, &(HostParams.t_end), sizeof(real_number));
    //     cudaMemcpyToSymbol(CFLnumber, &(HostParams.CFLnumber), sizeof(real_number));
    //     cudaMemcpyToSymbol(write_const, &(HostParams.write_const), sizeof(real_number));
    //     cudaMemcpyToSymbol(ObstacleCenter, &(HostParams.ObstacleCenter), sizeof(real_number) * DIM);
    //     cudaMemcpyToSymbol(ObstacleVelocity, &(HostParams.ObstacleVelocity), sizeof(real_number) * DIM);
    //     cudaMemcpyToSymbol(ObstacleOmega, &(HostParams.ObstacleOmega), sizeof(real_number));
    //     cudaMemcpyToSymbol(Rin, &(HostParams.Rin), sizeof(real_number));
    //     cudaMemcpyToSymbol(Rout, &(HostParams.Rout), sizeof(real_number));
    //     cudaMemcpyToSymbol(Win, &(HostParams.Win), sizeof(real_number));
    //     cudaMemcpyToSymbol(Wout, &(HostParams.Wout), sizeof(real_number));
    //     cudaMemcpyToSymbol(ObstacleBase, &(HostParams.ObstacleBase), sizeof(real_number));
    //     cudaMemcpyToSymbol(ObstacleHeight, &(HostParams.ObstacleHeight), sizeof(real_number));
    //     cudaMemcpyToSymbol(ObstacleTilt, &(HostParams.ObstacleTilt), sizeof(real_number));
    // }

    void WriteParameters(const std::string filename)
    {
        // write all the parameters to a .txt file for debugging and making sure the parameters are read correctly
        std::ofstream file;
        file.open(filename);
        file << "SCENARIO: " << SCENARIO << std::endl;
        file << "BC_TYPE: " << BC_TYPE << std::endl;
        file << "DENSITY_TYPE: " << DENSITY_TYPE << std::endl;
        file << "WRITER: " << WRITER << std::endl;
        file << "length={" << length[0] << "," << length[1] << "}" << std::endl;
        file << "bc={" << bc[0] << "," << bc[1] << "}" << std::endl;
        file << "Nfluid={" << Nfluid[0] << "," << Nfluid[1] << "}" << std::endl;
        file << "Nboundary={" << Nboundary[0] << "," << Nboundary[1] << "}" << std::endl;
        file << "LengthScale: " << LengthScale << std::endl;
        file << "rf: " << rf << std::endl;
        file << "H: " << H << std::endl;
        file << "r_cut: " << r_cut << std::endl;
        file << "r_cut2: " << r_cut2 << std::endl;
        file << "Kquintic: " << Kquintic << std::endl;
        file << "Re: " << Re << std::endl;
        file << "umax: " << umax << std::endl;
        file << "rho0: " << rho0 << std::endl;
        file << "gamma: " << gamma << std::endl;
        file << "coeff_sound: " << coeff_sound << std::endl;
        file << "cbar: " << cbar << std::endl;
        file << "B: " << B << std::endl;
        file << "xi: " << xi << std::endl;
        file << "gravity_vector: {" << gravity_vector[0] << "," << gravity_vector[1] << "}" << std::endl;
        file << "gravity: " << gravity << std::endl;
        file << "vw_top: {" << vw_top[0] << "," << vw_top[1] << "}" << std::endl;
        file << "vw_bottom: {" << vw_bottom[0] << "," << vw_bottom[1] << "}" << std::endl;
        file << "MassFluid: " << MassFluid << std::endl;
        file << "MassBound: " << MassBound << std::endl;
        file << "nu: " << nu << std::endl;
        file << "eta: " << eta << std::endl;
        file << "Bfactor: " << Bfactor << std::endl;
        file << "Pbackground: " << Pbackground << std::endl;
        file << "t_end: " << t_end << std::endl;
        file << "CFLnumber: " << CFLnumber << std::endl;
        file << "write_const: " << write_const << std::endl;
        file << "ObstacleCenter: {" << ObstacleCenter[0] << "," << ObstacleCenter[1] << "}" << std::endl;
        file << "ObstacleVelocity: {" << ObstacleVelocity[0] << "," << ObstacleVelocity[1] << "}" << std::endl;
        file << "ObstacleOmega: " << ObstacleOmega << std::endl;
        file << "Rin: " << Rin << std::endl;
        file << "Rout: " << Rout << std::endl;
        file << "Win: " << Win << std::endl;
        file << "Wout: " << Wout << std::endl;
        file << "ObstacleBase: " << ObstacleBase << std::endl;
        file << "ObstacleHeight: " << ObstacleHeight << std::endl;
        file << "ObstacleTilt: " << ObstacleTilt << std::endl;
    }
};

#ifdef __CUDACC__
__constant__ Parameters _params_gpu_;
#else
extern Parameters _params_gpu_;
#endif

#endif // DEFINITIONS_H
