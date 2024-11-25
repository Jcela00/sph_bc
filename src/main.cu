#define PRINT_STACKTRACE
#define OPENMPI

#include <math.h>
#include <sys/stat.h>
#include <cmath>

#include "CreateParticleGeometry.hpp"
#include "DragLiftCoeficient.hpp"
#include "InitializeParameters.hpp"
#include "VectorUtilities.hpp"
#include "Obstacle.hpp"
#include "Kernel.hpp"
#include "Probes.hpp"
#include "Physics.hpp"
#include "Calculations.hpp"
#include "TimeIntegration.hpp"
#include "CalcForces.hpp"
#include "ExtrapolateVelocity.hpp"
#include "Definitions.hpp"

#ifdef __CUDACC__
#else
Parameters _params_gpu_;
#endif

int main(int argc, char *argv[])
{

	// initialize the library
	openfpm_init(&argc, &argv);

	// check if DIM is 2 or 3, if not, throw an error
	if (DIM != 2 && DIM != 3)
	{
		std::cerr << "Error: DIM must be 2 or 3" << std::endl;
		return -1;
	}

#if !defined(CUDA_ON_CPU) && !defined(__HIP__)
	cudaDeviceSetCacheConfig(cudaFuncCachePreferL1);
#endif

#ifdef CUDIFY_USE_CUDA
	cudaDeviceSetCacheConfig(cudaFuncCachePreferL1);
#endif

	// Get imput xml filename
	std::string tmp = argv[1];

	// Initialize parameter objects
	Parameters MainParameters;
	AuxiliarParameters AuxParameters;

	// Read xml file and initialize simulation parameters
	ParseXMLFile(tmp, MainParameters, AuxParameters);
	InitializeConstants(MainParameters, AuxParameters);

	// Write parameters to file
	MainParameters.WriteParameters("Parameters.txt");

// Copy parameters object to gpu __constant__ memory to be accesible in device code
#ifdef __CUDACC__
	cudaError_t err = cudaMemcpyToSymbol(_params_gpu_, &MainParameters, sizeof(Parameters));
	if (err != cudaSuccess)
	{
		std::cerr << "Error copying to constant memory: " << cudaGetErrorString(err) << std::endl;
	}
	cudaDeviceSynchronize();
#else
	_params_gpu_ = MainParameters;
#endif

	// Write GPU parameters to file to check if correctly copied
	Parameters host_params_check;
	cudaMemcpyFromSymbol(&host_params_check, _params_gpu_, sizeof(Parameters));
	cudaDeviceSynchronize();
	host_params_check.WriteParameters(("check_params.txt"));

	// Create the particle data structure
	particles vd;

	// Set names for particle properties
	openfpm::vector<std::string> names({"type", "rho", "pressure", "drho", "velocity", "velocity_t", "force",
										"force_t", "normal", "volume", "omega", "vorticity", "red_vel", "red_fx", "red_fy"});
	vd.setPropNames(names);

	// Create vector of pairs of probe particles data structure and integer that represents measured velocity component, 0 for x, 1 for y
	// Vector may contain multiple probe objects that write to different files, for example for cavity scenario we place a probe at x=0.5 and another at y=0.5
	std::vector<std::pair<probe_particles, int>> vp;

	// Pointer to obstacle class
	Obstacle *obstacle_ptr = nullptr;

	// Initialize problem geometry, more complex scenarios have specialized functions to create the geometry
	if (MainParameters.SCENARIO == STEP)
	{
		CreateParticleGeometryStep(vd, vp, MainParameters, AuxParameters);
	}
	else if (MainParameters.SCENARIO == TAYLOR_COUETTE)
	{
		CreateParticleGeometryTaylorCouette(vd, vp, obstacle_ptr, MainParameters, AuxParameters);
	}
	else if (MainParameters.SCENARIO == DAM_BREAK)
	{
		CreateParticleGeometryDamBreak(vd, vp, MainParameters, AuxParameters);
	}
	else if (MainParameters.SCENARIO == DAM_BREAK_ADJ)
	{
		CreateParticleGeometryDamBreakAdj(vd, vp, MainParameters, AuxParameters);
	}
	else if (MainParameters.SCENARIO == CAVITY)
	{
		CreateParticleGeometryCavity(vd, vp, obstacle_ptr, MainParameters, AuxParameters);
	}
	else if (MainParameters.SCENARIO == SPHERE)
	{
		CreateParticleGeometrySphere(vd, vp, obstacle_ptr, MainParameters, AuxParameters);
	}
	else
	{
		CreateParticleGeometry(vd, vp, obstacle_ptr, MainParameters, AuxParameters);
	}
	vd.map();

	for (size_t k = 0; k < vp.size(); k++)
	{
		(vp[k].first).map();
	}

	// For load balancing in MPI case, no longer needed
	// Now that we fill the vector with particles, decompose the domain
	// ModelCustom md;
	// vd.addComputationCosts(md);
	// vd.getDecomposition().decompose();
	// vd.map();

	// for (int k = 0; k < vp.size(); k++)
	// {
	// 	(vp[k].first).getDecomposition().decompose();
	// 	(vp[k].first).map();
	// }

	vd.ghost_get<vd0_type, vd1_rho, vd2_pressure, vd4_velocity, vd5_velocity_t, vd6_force, vd7_force_t, vd8_normal, vd9_volume, vd10_omega, vd11_vorticity>();

	// Create a Cell list object for CPU calculations
	auto NN = vd.getCellList(MainParameters.r_cut);
	vd.ghost_get<vd0_type, vd1_rho, vd2_pressure, vd4_velocity, vd5_velocity_t, vd6_force, vd7_force_t, vd8_normal, vd9_volume, vd10_omega, vd11_vorticity>();
	vd.updateCellList(NN);

	// Set up fluid vector, normal vector, curvature and volume of the boundary particles
	if (MainParameters.BC_TYPE == NEW_NO_SLIP)
	{
		CalcFluidVec(vd, NN, MainParameters);
		vd.ghost_get<vd0_type, vd1_rho, vd2_pressure, vd4_velocity, vd5_velocity_t, vd6_force, vd7_force_t, vd8_normal, vd9_volume, vd10_omega, vd11_vorticity>();
		// vd.write_frame(AuxParameters.filename + "FluidVec", 0, MainParameters.WRITER);

		CalcNormalVec(vd, NN, MainParameters);
		vd.ghost_get<vd0_type, vd1_rho, vd2_pressure, vd4_velocity, vd5_velocity_t, vd6_force, vd7_force_t, vd8_normal, vd9_volume, vd10_omega, vd11_vorticity>();
		// vd.write_frame(AuxParameters.filename + "NormalVec", 0, MainParameters.WRITER);

		CalcCurvature(vd, NN, MainParameters);
		vd.ghost_get<vd0_type, vd1_rho, vd2_pressure, vd4_velocity, vd5_velocity_t, vd6_force, vd7_force_t, vd8_normal, vd9_volume, vd10_omega, vd11_vorticity>();
		vd.write_frame(AuxParameters.filename + "_Curvature", 0, static_cast<double>(0.0), MainParameters.WRITER);

		CalcVolume(vd, MainParameters.dp);
		vd.ghost_get<vd0_type, vd1_rho, vd2_pressure, vd4_velocity, vd5_velocity_t, vd6_force, vd7_force_t, vd8_normal, vd9_volume, vd10_omega, vd11_vorticity>();
	}
	vd.map();

	// Write to file
	vd.write_frame(AuxParameters.filename, 0, static_cast<double>(0.0), MainParameters.WRITER);

	// Move particle data structure to gpu, and create gpu cell list
	vd.hostToDevicePos();
	vd.template hostToDeviceProp<vd0_type, vd1_rho, vd2_pressure, vd3_drho, vd4_velocity, vd5_velocity_t, vd6_force, vd7_force_t, vd8_normal, vd9_volume, vd10_omega, vd11_vorticity, vd12_vel_red, vd13_force_red_x, vd14_force_red_y>();

	vd.map(RUN_ON_DEVICE);
	vd.ghost_get<vd0_type, vd1_rho, vd2_pressure, vd3_drho, vd4_velocity, vd5_velocity_t, vd6_force, vd7_force_t, vd8_normal, vd9_volume, vd10_omega, vd11_vorticity, vd12_vel_red, vd13_force_red_x, vd14_force_red_y>(RUN_ON_DEVICE);
	auto NN_gpu = vd.getCellListGPU(MainParameters.r_cut / 2.0, CL_NON_SYMMETRIC | CL_GPU_REORDER, 2);

	// Time and write output counters
	real_number t = 0.0;
	size_t write = 0;

	// Last write counter value
	size_t last_write = static_cast<size_t>(MainParameters.write_const * MainParameters.t_end);

	// Stream to output drag and lift data
	std::ofstream DragFile(AuxParameters.filename + "_DragLift.csv");

	// Time simulation
	timer tot_sim;
	tot_sim.start();

	// Like mpi communicator
	Vcluster<> &v_cl = create_vcluster();

	// Time loop
	while (t <= MainParameters.t_end)
	{

		// Dynamic load balancing for MPI version, no longer used in CUDA
		//// Do rebalancing every 200 timesteps
		// it_reb++;
		// if (it_reb == 500)
		// {
		// 	vd.map();

		// 	it_reb = 0;
		// 	ModelCustom md;
		// 	vd.addComputationCosts(md);
		// 	vd.getDecomposition().decompose();
		// 	vd.map();

		// 	if (v_cl.getProcessUnitID() == 0)
		// 		std::cout << "REBALANCED " << std::endl;
		// }

		// Calculate dt for time stepping
		real_number dt = calc_deltaT(vd, MainParameters);

		// Integrate one time step
		kick_drift_int(vd, NN_gpu, dt, t, MainParameters, AuxParameters);

		// Increment time
		t += dt;

		// Write drag and lift to a file
		if (write < t * MainParameters.write_const * 10)
		{
			// make necessary reductions for computing drag and lift
			Point<3, real_number> VxDragLift = ComputeDragLift(vd, MainParameters);
			CalcDragLift(t, DragFile, VxDragLift, MainParameters, write);
		}

		// Write the particle configuration to file
		if (write < t * MainParameters.write_const)
		{
			// Send data from GPU to CPU for writing to file
			vd.deviceToHostPos();
			vd.deviceToHostProp<vd0_type, vd1_rho, vd2_pressure, vd3_drho, vd4_velocity, vd5_velocity_t, vd6_force, vd7_force_t, vd8_normal, vd9_volume, vd10_omega, vd11_vorticity>();

			// Update CPU cell list for computing vorticity and probes
			vd.updateCellList(NN);

			if constexpr (DIM == 2)
			{
				CalcVorticity(vd, NN, MainParameters);
			}
			else if constexpr (DIM == 3)
			{
				// vorticity 3d
			}

			// Write probes to file (only last time step, we are only interested in the lattice and cavity case)
			// Probe calculation requires ghost and updated cell-list
			if (MainParameters.PROBES_ENABLED && write == last_write)
			{
				for (size_t k = 0; k < vp.size(); k++) // for each probe
				{
					probe_particles &probe = vp[k].first;
					int &component = vp[k].second;
					sensor_velocity_comp(vd, probe, NN, component, obstacle_ptr, MainParameters);
					probe.write_frame(AuxParameters.probe_filenames[k], write, MainParameters.WRITER);
				}
			}

			vd.deleteGhost();
			vd.write_frame(AuxParameters.filename, write, static_cast<double>(t), MainParameters.WRITER);
			vd.ghost_get<vd0_type, vd1_rho, vd2_pressure, vd4_velocity, vd5_velocity_t, vd8_normal, vd9_volume, vd10_omega>(RUN_ON_DEVICE);

			write++;

			if (v_cl.getProcessUnitID() == 0)
			{
				std::cout << "TIME: " << t << " step " << AuxParameters.cnt << std::endl;
			}
		}
	}
	tot_sim.stop();
	std::cout << "Time to complete: " << tot_sim.getwct() << " seconds" << std::endl;

	// Calculate final density statistics to know density fluctuation percentage
	vd.deviceToHostProp<vd0_type, vd1_rho>();

	real_number max_rho, min_rho, avg_rho;
	unsigned int count = 0;
	max_rho = -1.0e10;
	min_rho = 1.0e10;
	auto part = vd.getDomainIterator();
	while (part.isNext())
	{
		auto a = part.get();

		if (vd.template getProp<vd0_type>(a) == FLUID)
		{
			real_number rho = vd.template getProp<vd1_rho>(a);
			avg_rho += rho;
			count++;

			if (rho > max_rho)
			{
				max_rho = rho;
			}
			if (rho < min_rho)
			{
				min_rho = rho;
			}
		}
		++part;
	}

	v_cl.max(max_rho);
	v_cl.min(min_rho);
	v_cl.sum(avg_rho);
	v_cl.sum(count);
	v_cl.execute();
	avg_rho /= count;

	real_number deltaRho = max_rho - min_rho;
	real_number percentFluctuation = deltaRho / MainParameters.rho0 * 100.0;
	std::cout << "Max density: " << max_rho << " Min density: " << min_rho << "average density: " << avg_rho << " Delta rho: " << deltaRho << " % Fluctuation: " << percentFluctuation << std::endl;

	delete obstacle_ptr;
	openfpm_finalize();
}
