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

	std::string tmp = argv[1];

	Parameters MainParameters;
	AuxiliarParameters AuxParameters;

	ParseXMLFile(tmp, MainParameters, AuxParameters);
	InitializeConstants(MainParameters, AuxParameters);

	MainParameters.WriteParameters("Parameters.txt");

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

	//////////////////////// CHECK CORRECT GPU PARAMS COPY /////////////////////
	Parameters host_params_check;
	cudaMemcpyFromSymbol(&host_params_check, _params_gpu_, sizeof(Parameters));
	cudaDeviceSynchronize();
	host_params_check.WriteParameters(("check_params.txt"));
	////////////////////////////////////////////////////////////////////////////

	// Create a particle vector
	particles vd;

	// set names
	openfpm::vector<std::string> names({"type", "rho", "pressure", "drho", "velocity", "velocity_t", "force",
										"force_t", "normal", "volume", "omega", "vorticity", "red_vel", "red_fx", "red_fy"});
	vd.setPropNames(names);

	std::vector<std::pair<probe_particles, int>> vp;

	Obstacle *obstacle_ptr = nullptr;
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
	vd.write_frame("After Create Geometry", 0, MainParameters.WRITER);

	for (size_t k = 0; k < vp.size(); k++)
	{
		(vp[k].first).map();
	}

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

	auto NN = vd.getCellList(MainParameters.r_cut);
	vd.ghost_get<vd0_type, vd1_rho, vd2_pressure, vd4_velocity, vd5_velocity_t, vd6_force, vd7_force_t, vd8_normal, vd9_volume, vd10_omega, vd11_vorticity>();
	vd.updateCellList(NN);

	if (MainParameters.BC_TYPE == NEW_NO_SLIP) // Set up fluid vector and normal vector of the boundary particles
	{
		CalcFluidVec(vd, NN, MainParameters);
		vd.ghost_get<vd0_type, vd1_rho, vd2_pressure, vd4_velocity, vd5_velocity_t, vd6_force, vd7_force_t, vd8_normal, vd9_volume, vd10_omega, vd11_vorticity>();
		// vd.write_frame("After Fluid Vec", 0, MainParameters.WRITER);

		CalcNormalVec(vd, NN, MainParameters);
		vd.ghost_get<vd0_type, vd1_rho, vd2_pressure, vd4_velocity, vd5_velocity_t, vd6_force, vd7_force_t, vd8_normal, vd9_volume, vd10_omega, vd11_vorticity>();
		// vd.write_frame("After Normal Vec", 0, MainParameters.WRITER);

		CalcCurvature(vd, NN, MainParameters);
		vd.ghost_get<vd0_type, vd1_rho, vd2_pressure, vd4_velocity, vd5_velocity_t, vd6_force, vd7_force_t, vd8_normal, vd9_volume, vd10_omega, vd11_vorticity>();
		vd.write_frame("Curvature", 0, MainParameters.WRITER);

		CalcVolume(vd, MainParameters.dp);
		vd.ghost_get<vd0_type, vd1_rho, vd2_pressure, vd4_velocity, vd5_velocity_t, vd6_force, vd7_force_t, vd8_normal, vd9_volume, vd10_omega, vd11_vorticity>();
	}
	vd.map();
	vd.write_frame("Initialization", 0, MainParameters.WRITER);
	vd.write_frame(AuxParameters.filename, 0, static_cast<double>(0.0), MainParameters.WRITER);

	// move to gpu
	vd.hostToDevicePos();
	vd.template hostToDeviceProp<vd0_type, vd1_rho, vd2_pressure, vd3_drho, vd4_velocity, vd5_velocity_t, vd6_force, vd7_force_t, vd8_normal, vd9_volume, vd10_omega, vd11_vorticity, vd12_vel_red, vd13_force_red_x, vd14_force_red_y>();

	vd.map(RUN_ON_DEVICE);
	vd.ghost_get<vd0_type, vd1_rho, vd2_pressure, vd3_drho, vd4_velocity, vd5_velocity_t, vd6_force, vd7_force_t, vd8_normal, vd9_volume, vd10_omega, vd11_vorticity, vd12_vel_red, vd13_force_red_x, vd14_force_red_y>(RUN_ON_DEVICE);
	auto NN_gpu = vd.getCellListGPU(MainParameters.r_cut / 2.0, CL_NON_SYMMETRIC | CL_GPU_REORDER, 2);

	if (MainParameters.BC_TYPE == NO_SLIP) // set up boundary particle velocity for the first iteration
	{
		ExtrapolateVelocity(vd, NN_gpu);
		vd.ghost_get<vd0_type, vd1_rho, vd2_pressure, vd4_velocity, vd5_velocity_t, vd6_force, vd7_force_t, vd8_normal, vd9_volume, vd10_omega, vd11_vorticity>(RUN_ON_DEVICE);
	}

	// Evolve
	size_t write = 0;
	real_number t = 0.0;
	std::ofstream avgvelstream(AuxParameters.filename + "_DragLift.csv");
	timer tot_sim, tot_step;
	tot_sim.start();

	Vcluster<> &v_cl = create_vcluster();
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
		real_number dt = calc_deltaT(vd, MainParameters);

		// Integrate one time step
		kick_drift_int(vd, NN_gpu, dt, t, MainParameters, AuxParameters);

		// increment time
		t += dt;

		// write the configuration
		if (write < t * MainParameters.write_const)
		{
			std::string timestring = "time" + std::to_string(t);

			// make necessary reductions for computing drag and lift
			Point<3, real_number> VxDragLift = ComputeDragLift(vd, MainParameters);
			// send data from GPU to CPU for writing to file
			// vd.map(RUN_ON_DEVICE);
			vd.deviceToHostPos();
			vd.deviceToHostProp<vd0_type, vd1_rho, vd2_pressure, vd3_drho, vd4_velocity, vd5_velocity_t, vd6_force, vd7_force_t, vd8_normal, vd9_volume, vd10_omega, vd11_vorticity>();

			// Update CPU cell list for computing vorticity ( and probes )
			if (MainParameters.SCENARIO == ELLIPSE || MainParameters.PROBES_ENABLED)
			{
				vd.updateCellList(NN);
			}

			// compute vorticity
			if (MainParameters.SCENARIO == ELLIPSE)
			{
				CalcVorticity(vd, NN, MainParameters);
			}

			// sensor calculation requires ghost and updated cell-list
			if (MainParameters.PROBES_ENABLED)
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
			CalcDragLift(t, avgvelstream, VxDragLift, MainParameters, write);

			write++;

			if (v_cl.getProcessUnitID() == 0)
			{
				std::cout << "TIME: " << t << " step " << AuxParameters.cnt << std::endl;
			}
		}
	}
	tot_sim.stop();
	std::cout << "Time to complete: " << tot_sim.getwct() << " seconds" << std::endl;

	vd.deviceToHostProp<vd0_type, vd1_rho>();

	// get max density
	// max_rho = reduce_local<vd0_type, _max_>(vd);

	// // change sign of density for fluid particles to get min using max reduction
	real_number max_rho, min_rho, avg_rho;
	int count = 0;
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

	// min_rho = reduce_local<vd1_rho, max>(vd);
	v_cl.max(max_rho);
	v_cl.min(min_rho);
	v_cl.sum(avg_rho);
	v_cl.sum(count);
	v_cl.execute();
	avg_rho /= count;
	// min_rho = -min_rho;

	real_number deltaRho = max_rho - min_rho;
	real_number percentFluctuation = deltaRho / MainParameters.rho0 * 100.0;
	std::cout << "Max density: " << max_rho << " Min density: " << min_rho << "average density: " << avg_rho << " Delta rho: " << deltaRho << " % Fluctuation: " << percentFluctuation << std::endl;

	delete obstacle_ptr;
	openfpm_finalize();
}
