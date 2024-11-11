#ifndef PROBES_HPP
#define PROBES_HPP

#include "Definitions.hpp"
#include "VectorUtilities.hpp"
#include "Kernel.hpp"
#include "Obstacle.hpp"

// void interact_probe_boundary_new(particles &vd,
// 								 vect_dist_key_dx probekey,
// 								 const Point<DIM, real_number> &r_wall_to_probe,
// 								 unsigned long &boundary_key,
// 								 const int component,
// 								 real_number &W_sum,
// 								 real_number &magnitude_tmp,
// 								 const Parameters &params);

void PlaceProbes(probe_particles &probes,
				 const int k0,
				 const int kmax,
				 const Point<DIM, real_number> Corner,
				 const Point<DIM, real_number> UnitOffset,
				 Obstacle *obstacle_ptr,
				 const std::vector<int> FixedProbeIndices = {},
				 const std::vector<real_number> FixedProbeValues = {});

template <typename CellList>
inline void sensor_velocity_comp(particles &vd,
								 probe_particles &probes,
								 CellList &NN,
								 const int component,
								 Obstacle *obstacle_ptr,
								 const Parameters &params)
{

	Vcluster<> &v_cl = create_vcluster();

	auto probe_it = probes.getDomainIterator();
	while (probe_it.isNext()) // iterate all probes
	{
		auto probekey = probe_it.get();

		if (probes.template getProp<probe0_type>(probekey) == FIXED_PROBE) // if probe is fixed value, skip
		{
			++probe_it;
			continue;
		}

		Point<DIM, real_number> xp = probes.getPos(probekey);
		real_number magnitude_tmp = 0.0;
		real_number W_sum = 0.0;

		// get the iterator over the neighbohood particles of the probe position
		auto itg = NN.getNNIteratorBox(NN.getCell(xp));

		while (itg.isNext())
		{
			// get key of the particle
			auto q = itg.get();

			// Get the position of the neighborhood particle q
			Point<DIM, real_number> xq = vd.getPos(q);
			Point<DIM, real_number> dr = xp - xq;
			// Calculate the distance
			real_number r2 = norm2(dr);

			if (r2 < params.r_cut2) // if inside the kernel
			{
				real_number r = sqrt(r2);

				if (vd.template getProp<vd0_type>(q) != FLUID)
				{
					if (r < 0.1 * params.dp) // if probe is placed on top of a wall particle
					{
						W_sum = 1.0;
						magnitude_tmp = vd.template getProp<vd4_velocity>(q)[component];
						break;
					}

					// if (params.BC_TYPE == NEW_NO_SLIP)
					// {
					// 	interact_probe_boundary_new(vd, probekey, dr, q, component, W_sum, magnitude_tmp);
					// }
					// else if (params.BC_TYPE == NO_SLIP)
					// {
					// 	real_number ker = Wab(r, params.H, params.Kquintic) * (params.MassBound / params.rho0);
					// 	W_sum += ker;
					// 	magnitude_tmp += vd.template getProp<vd5_velocity_t>(q)[component] * ker;
					// }
				}
				else
				{
					real_number rhoq = vd.template getProp<vd1_rho>(q);
					real_number ker = Wab(r, params.H, params.Kquintic) * (params.MassFluid / rhoq);

					// Add the total magnitude contribution
					magnitude_tmp += vd.template getProp<vd4_velocity>(q)[component] * ker;

					// Also keep track of the calculation of the summed kernel
					W_sum += ker;
				}
			}
			// next neighborhood particle
			++itg;
		}

		// We calculate the magnitude normalizing with W_sum
		if (W_sum == 0.0)
			magnitude_tmp = 0.0;
		else
			magnitude_tmp = magnitude_tmp / W_sum;

		probes.getProp<probe1_quantity>(probekey) = magnitude_tmp;
		++probe_it;
	}
}

#endif // PROBES_HPP
