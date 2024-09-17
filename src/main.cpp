#include "Vector/vector_dist.hpp"
#include "Draw/DrawParticles.hpp"
#include <math.h>
#include <sys/stat.h>
#include <cmath>

#include "Definitions.hpp"
#include "VectorUtilities.hpp"
#include "Kernel.hpp"
#include "Physics.hpp"
#include "Calculations.hpp"
#include "Interactions.hpp"
#include "CalcForces.hpp"
#include "Obstacle.hpp"
#include "TimeIntegration.hpp"
#include "Probes.hpp"
#include "CreateParticleGeometry.hpp"
#include "InitializeParameters.hpp"
#include "DragLiftCoeficient.hpp"

int main(int argc, char *argv[])
{

	// initialize the library
	openfpm_init(&argc, &argv);

	std::string tmp = argv[1];

	Parameters MainParameters;

	// create a Vcluster object ( like MPI communicator )
	Vcluster<> &v_cl = create_vcluster();

	ParseXMLFile(tmp, MainParameters);
	InitializeConstants(v_cl, MainParameters);

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
		CreateParticleGeometryStep(vd, vp, v_cl, MainParameters);
	}
	else if (MainParameters.SCENARIO == TAYLOR_COUETTE)
	{
		CreateParticleGeometryTaylorCouette(vd, vp, v_cl, obstacle_ptr, MainParameters);
	}
	else
	{
		CreateParticleGeometry(vd, vp, v_cl, obstacle_ptr, MainParameters);
	}
	vd.map();
	// vd.write_frame(MainParameters.filename, 0, MainParameters.WRITER);

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

	vd.ghost_get<type, rho, pressure, force, velocity, force_transport, v_transport, normal_vector, curvature_boundary, arc_length, vd_volume, vd_omega>();

	auto NN = vd.getCellList(MainParameters.r_threshold);
	vd.ghost_get<type, rho, pressure, force, velocity, force_transport, v_transport, normal_vector, curvature_boundary, arc_length, vd_volume, vd_omega>();
	vd.updateCellList(NN);

	if (MainParameters.BC_TYPE == NO_SLIP) // set up boundary particle velocity for the first iteration
	{
		ExtrapolateVelocity(vd, NN, MainParameters);
		vd.ghost_get<type, rho, pressure, force, velocity, force_transport, v_transport, normal_vector, curvature_boundary, arc_length, vd_volume, vd_omega>();
	}
	else if (MainParameters.BC_TYPE == NEW_NO_SLIP) // Set up fluid vector and normal vector of the boundary particles
	{
		CalcFluidVec(vd, NN, MainParameters);
		vd.ghost_get<type, rho, pressure, force, velocity, force_transport, v_transport, normal_vector, curvature_boundary, arc_length, vd_volume, vd_omega>();

		CalcNormalVec(vd, NN, MainParameters);
		vd.ghost_get<type, rho, pressure, force, velocity, force_transport, v_transport, normal_vector, curvature_boundary, arc_length, vd_volume, vd_omega>();

		CalcCurvature(vd, NN, v_cl, MainParameters);
		vd.ghost_get<type, rho, pressure, force, velocity, force_transport, v_transport, normal_vector, curvature_boundary, arc_length, vd_volume, vd_omega>();

		CalcVolume(vd, MainParameters.dp);
		vd.ghost_get<type, rho, pressure, force, velocity, force_transport, v_transport, normal_vector, curvature_boundary, arc_length, vd_volume, vd_omega>();
	}

	// Evolve
	size_t write = 0;

	double t = 0.0;
	std::ofstream avgvelstream(MainParameters.filename + "_DragLift.csv");
	bool calc_drag = false;
	double obstacle_force_x = 0.0;
	double obstacle_force_y = 0.0;

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
		double dt = calc_deltaT(vd, v_cl, MainParameters);

		// in general we dont compute drag coeficient every time step
		calc_drag = false;
		if (write < (t + dt) * MainParameters.write_const)
		{
			calc_drag = true;
		}

		// Integrate one time step
		kick_drift_int(vd, NN, dt, v_cl, obstacle_force_x, obstacle_force_y, calc_drag, t, MainParameters);

		// increment time
		t += dt;
		if (write < t * MainParameters.write_const)
		{
			// // sensor calculation require ghost and update cell-list
			if (MainParameters.PROBES_ENABLED)
			{
				vd.map();
				vd.ghost_get<type, rho, pressure, force, velocity, force_transport, v_transport, normal_vector, curvature_boundary, arc_length, vd_volume, vd_omega>();

				vd.updateCellList(NN);
				for (size_t k = 0; k < vp.size(); k++)
				{
					probe_particles &probe = vp[k].first;
					// probe.map();
					int &component = vp[k].second;
					sensor_velocity_comp(vd, probe, v_cl, NN, component, obstacle_ptr, MainParameters);
					probe.write_frame(MainParameters.probe_filenames[k], write, MainParameters.WRITER);
				}
			}

			CalcVorticity(vd, NN, MainParameters);
			vd.deleteGhost();
			vd.write_frame(MainParameters.filename, write, MainParameters.WRITER);
			vd.ghost_get<type, rho, pressure, force, velocity, force_transport, v_transport, normal_vector, curvature_boundary, arc_length, vd_volume, vd_omega>();

			CalcDragLift(vd, v_cl, t, avgvelstream, obstacle_force_x, obstacle_force_y, MainParameters, write);
			write++;
			if (v_cl.getProcessUnitID() == 0)
			{
				std::cout << "TIME: " << t << "  write " << MainParameters.cnt << std::endl;
			}
		}
		// else
		// {
		// 	if (v_cl.getProcessUnitID() == 0)
		// 	{
		// 		std::cout << "TIME: " << t << "   " << cnt << std::endl;
		// 	}
		// }
	}
	delete obstacle_ptr;
	openfpm_finalize();
}
