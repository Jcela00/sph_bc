#ifdef __NVCC__
#define PRINT_STACKTRACE
#define OPENMPI

#include <math.h>
#include <sys/stat.h>
#include <cmath>

#include "Definitions.hpp"
#include "CreateParticleGeometry.hpp"
#include "DragLiftCoeficient.hpp"
#include "InitializeParameters.hpp"
#include "VectorUtilities.hpp"
#include "Kernel.hpp"
#include "Physics.hpp"
#include "Calculations.hpp"
#include "Obstacle.hpp"
#include "Probes.hpp"
#include "TimeIntegration.hpp"
#include "CalcForces.hpp"
#include "ExtrapolateVelocity.hpp"

int main(int argc, char *argv[])
{

	// initialize the library
	openfpm_init(&argc, &argv);

#if !defined(CUDA_ON_CPU) && !defined(__HIP__)
	cudaDeviceSetCacheConfig(cudaFuncCachePreferL1);
#endif

	std::string tmp = argv[1];

	Parameters MainParameters;
	AuxiliarParameters AuxParameters;

	// create a Vcluster object ( like MPI communicator )
	Vcluster<> &v_cl = create_vcluster();

	ParseXMLFile(tmp, MainParameters, AuxParameters);
	InitializeConstants(v_cl, MainParameters, AuxParameters);

	MainParameters.WriteParameters("Parameters.txt");

	cudaError_t err = cudaMemcpyToSymbol(_params_gpu_, &MainParameters, sizeof(Parameters));
	if (err != cudaSuccess)
	{
		std::cerr << "Error copying to constant memory: " << cudaGetErrorString(err) << std::endl;
	}
	cudaDeviceSynchronize();

	// Parameters host_params_check;
	// cudaMemcpyFromSymbol(&host_params_check, _params_gpu_, sizeof(Parameters));
	// cudaDeviceSynchronize();
	// host_params_check.WriteParameters(("check_params.txt"));

	// Create a particle vector
	particles vd;
	// set names
	openfpm::vector<std::string> names({"type", "rho", "pressure", "drho", "force", "velocity", "force_transport",
										"v_transport", "normal", "curvature", "arc_length", "volume", "vd_omega", "vorticity"});
	vd.setPropNames(names);

	std::vector<std::pair<probe_particles, int>> vp;

	Obstacle *obstacle_ptr = nullptr;
	if (MainParameters.SCENARIO == STEP)
	{
		CreateParticleGeometryStep(vd, vp, v_cl, MainParameters, AuxParameters);
	}
	else if (MainParameters.SCENARIO == TAYLOR_COUETTE)
	{
		CreateParticleGeometryTaylorCouette(vd, vp, v_cl, obstacle_ptr, MainParameters, AuxParameters);
	}
	else
	{
		CreateParticleGeometry(vd, vp, v_cl, obstacle_ptr, MainParameters, AuxParameters);
	}
	vd.map();

	// for (size_t k = 0; k < vp.size(); k++)
	// {
	// 	(vp[k].first).map();
	// }

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

	vd.ghost_get<type, rho, pressure, force, velocity, force_transport, v_transport, normal_vector, curvature_boundary, arc_length, vd_volume, vd_omega>();

	auto NN = vd.getCellList(MainParameters.r_cut);
	vd.ghost_get<type, rho, pressure, force, velocity, force_transport, v_transport, normal_vector, curvature_boundary, arc_length, vd_volume, vd_omega>();
	vd.updateCellList(NN);

	if (MainParameters.BC_TYPE == NEW_NO_SLIP) // Set up fluid vector and normal vector of the boundary particles
	{
		CalcFluidVec(vd, NN, MainParameters);
		vd.ghost_get<type, rho, pressure, force, velocity, force_transport, v_transport, normal_vector, curvature_boundary, arc_length, vd_volume, vd_omega>();

		CalcNormalVec(vd, NN, MainParameters);
		vd.ghost_get<type, rho, pressure, force, velocity, force_transport, v_transport, normal_vector, curvature_boundary, arc_length, vd_volume, vd_omega>();

		CalcCurvature(vd, NN, MainParameters);
		vd.ghost_get<type, rho, pressure, force, velocity, force_transport, v_transport, normal_vector, curvature_boundary, arc_length, vd_volume, vd_omega>();

		CalcVolume(vd, MainParameters.dp);
		vd.ghost_get<type, rho, pressure, force, velocity, force_transport, v_transport, normal_vector, curvature_boundary, arc_length, vd_volume, vd_omega>();
	}
	vd.map();
	vd.write_frame("Initialization", 0, MainParameters.WRITER);

	// // regularize configuration
	// vd.write_frame("reg", 0, MainParameters.WRITER);
	// real_number dt_tmp = 10 * calc_deltaT(vd, v_cl, MainParameters);
	// size_t Niters = 1e4;
	// size_t iter = 0;
	// while (iter < Niters)
	// {
	// 	kick_drift_int_regularize(vd, NN, dt_tmp, v_cl, iter, MainParameters);
	// 	if (iter % 100 == 0)
	// 		vd.write_frame("reg", iter, MainParameters.WRITER);
	// }
	// vd.write_frame("reg", iter, MainParameters.WRITER);

	// move to gpu

	vd.hostToDevicePos();
	vd.template hostToDeviceProp<type, rho, pressure, force, velocity, force_transport, v_transport, normal_vector, curvature_boundary, arc_length, vd_volume, vd_omega, vd_vorticity>();
	vd.map(RUN_ON_DEVICE);
	vd.ghost_get<type, rho, pressure, force, velocity, force_transport, v_transport, normal_vector, curvature_boundary, arc_length, vd_volume, vd_omega, vd_vorticity>(RUN_ON_DEVICE);
	auto NN_gpu = vd.getCellListGPU(MainParameters.r_cut / 2.0f, CL_NON_SYMMETRIC | CL_GPU_REORDER_POSITION | CL_GPU_REORDER_PROPERTY | CL_GPU_RESTORE_PROPERTY, 2);

	if (MainParameters.BC_TYPE == NO_SLIP) // set up boundary particle velocity for the first iteration
	{
		ExtrapolateVelocity(vd, NN_gpu);
		vd.ghost_get<type, rho, pressure, force, velocity, force_transport, v_transport, normal_vector, curvature_boundary, arc_length, vd_volume, vd_omega>(RUN_ON_DEVICE);
	}

	// Evolve
	size_t write = 0;
	real_number t = 0.0;
	// std::ofstream avgvelstream(MainParameters.filename + "_DragLift.csv");
	// bool calc_drag = false;
	// real_number obstacle_force_x = 0.0;
	// real_number obstacle_force_y = 0.0;

	timer tot_sim;
	tot_sim.start();
	while (t <= MainParameters.t_end)
	{

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
		real_number dt = calc_deltaT(MainParameters);

		// in general we dont compute drag coeficient every time step
		// calc_drag = false;
		// if (write < (t + dt) * MainParameters.write_const)
		// {
		// 	calc_drag = true;
		// }

		// Integrate one time step
		kick_drift_int(vd, NN_gpu, dt, t, MainParameters, AuxParameters);

		// increment time
		t += dt;
		if (write < t * MainParameters.write_const)
		{
			// // // sensor calculation require ghost and update cell-list
			// if (MainParameters.PROBES_ENABLED)
			// {
			// 	vd.map();
			// 	vd.ghost_get<type, rho, pressure, force, velocity, force_transport, v_transport, normal_vector, curvature_boundary, arc_length, vd_volume, vd_omega>();

			// 	vd.updateCellList(NN);
			// 	for (size_t k = 0; k < vp.size(); k++)
			// 	{
			// 		probe_particles &probe = vp[k].first;
			// 		// probe.map();
			// 		int &component = vp[k].second;
			// 		sensor_velocity_comp(vd, probe, v_cl, NN, component, obstacle_ptr, MainParameters);
			// 		probe.write_frame(MainParameters.probe_filenames[k], write, MainParameters.WRITER);
			// 	}
			// }

			// send data from GPU to CPU
			vd.deviceToHostPos();
			vd.deviceToHostProp<type, rho, pressure, force, velocity, force_transport, v_transport, normal_vector, curvature_boundary, arc_length, vd_volume, vd_omega, vd_vorticity>();
			// Update CPU cell list for computing vorticity
			vd.updateCellList(NN);
			CalcVorticity(vd, NN, MainParameters);
			vd.deleteGhost();
			vd.write_frame(AuxParameters.filename, write, MainParameters.WRITER);
			vd.ghost_get<type, rho, pressure, force, velocity, force_transport, v_transport, normal_vector, curvature_boundary, arc_length, vd_volume, vd_omega, vd_vorticity>(RUN_ON_DEVICE);
			// CalcDragLift(vd, v_cl, t, avgvelstream, obstacle_force_x, obstacle_force_y, MainParameters, write);

			write++;

			if (v_cl.getProcessUnitID() == 0)
			{
				std::cout << "TIME: " << t << "  write " << AuxParameters.cnt << std::endl;
			}
		}
	}

	tot_sim.stop();
	std::cout << "Time to complete: " << tot_sim.getwct() << " seconds" << std::endl;

	delete obstacle_ptr;
	openfpm_finalize();
}

#else

int main(int argc, char *argv[])
{
	return 0;
}

#endif