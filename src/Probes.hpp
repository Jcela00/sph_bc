#ifndef PROBES_HPP
#define PROBES_HPP

#include "Definitions.hpp"
#include "VectorUtilities.hpp"
#include "Kernel.hpp"
#include "Obstacle.hpp"

void interact_probe_boundary_new(particles &vd,
								 vect_dist_key_dx probekey,
								 const Point<DIM, real_number> &r_wall_to_probe,
								 unsigned long &boundary_key,
								 const int component,
								 real_number &W_sum,
								 real_number &magnitude_tmp,
								 const Parameters &params);

void PlaceProbes(probe_particles &probes,
				 const int k0,
				 const int kmax,
				 const Point<DIM, real_number> Corner,
				 const Point<DIM, real_number> UnitOffset);

template <typename CellList>
inline void sensor_velocity_comp(particles &vd,
								 probe_particles &probes,
								 Vcluster<> &v_cl,
								 CellList &NN,
								 const int component,
								 Obstacle *obstacle_ptr,
								 const Parameters &params)
{

	auto it = probes.getDomainIterator();
	while (it.isNext()) // iterate all probes
	{

		auto probekey = it.get();
		Point<DIM, real_number> xp = probes.getPos(probekey);

		real_number magnitude_tmp = 0.0;
		real_number W_sum = 0.0;

		if (obstacle_ptr->isInside(xp)) // if probe is inside the obstacle it gets a 0 value
		{
			probes.getProp<0>(probekey) = 0.0;
			++it;
			continue;
		}
		else
		{
			// get the iterator over the neighbohood particles of the probes position
			auto itg = NN.getNNIteratorBox(NN.getCell(xp));
			while (itg.isNext())
			{
				// get key of the particle
				auto q = itg.get();

				// Get the position of the neighborhood particle q
				Point<DIM, real_number> xq = vd.getPos(q);
				Point<DIM, real_number> dr = xp - xq;
				// Calculate the distance
				real_number r = getVectorNorm(dr);

				if (vd.template getProp<type>(q) == BOUNDARY)
				{
					if (r < 0.01 * params.dp) // if probe is placed on top of a wall particle
					{
						W_sum = 1.0;
						magnitude_tmp = vd.template getProp<velocity>(q)[component];
						break;
					}

					// if (BC_TYPE == NEW_NO_SLIP)
					// {
					// 	// interact_probe_boundary_new(vd, probekey, dr, q, component, W_sum, magnitude_tmp);
					// }
					else if (params.BC_TYPE == NO_SLIP)
					{
						real_number ker = Wab(r, params.H, params.Kquintic) * (params.MassBound / params.rho0);
						W_sum += ker;
						magnitude_tmp += vd.template getProp<v_transport>(q)[component] * ker;
					}
				}
				else if (vd.template getProp<type>(q) == FLUID)
				{
					real_number rhoq = vd.template getProp<rho>(q);
					real_number ker = Wab(r, params.H, params.Kquintic) * (params.MassFluid / rhoq);

					// Also keep track of the calculation of the summed kernel
					W_sum += ker;

					// Add the total magnitude contribution
					magnitude_tmp += vd.template getProp<velocity>(q)[component] * ker;
				}
				// next neighborhood particle
				++itg;
			}

			// We calculate the magnitude normalizing with W_sum
			if (W_sum == 0.0)
				magnitude_tmp = 0.0;
			else
				magnitude_tmp = 1.0 / W_sum * magnitude_tmp;

			probes.getProp<0>(probekey) = magnitude_tmp;
			++it;
		}
	}
}

#endif // PROBES_HPP
