#include "Probes.hpp"

// void interact_probe_boundary_new(particles &vd,
//                                  vect_dist_key_dx probekey,
//                                  const Point<DIM, real_number> &r_wall_to_probe,
//                                  unsigned long &boundary_key,
//                                  const int component,
//                                  real_number &W_sum,
//                                  real_number &magnitude_tmp,
//                                  const Parameters &params)
// {
//     Point<DIM, real_number> r_probe_to_wall = -1.0 * r_wall_to_probe;

//     Point<DIM, real_number> normal = vd.getProp<normal_vector>(boundary_key);
//     // get curvature, and arc length
//     real_number kappa = vd.getProp<curvature_boundary>(boundary_key);
//     real_number dxwall = vd.getProp<arc_length>(boundary_key);

//     std::array<Point<DIM, real_number>, 3> r_boundary = getBoundaryPositions(r_probe_to_wall, normal, params.dp);
//     std::array<real_number, 3> r_boundary_norm = {getVectorNorm(r_boundary[0]), getVectorNorm(r_boundary[1]), getVectorNorm(r_boundary[2])};
//     std::array<real_number, 3> Volume_boundary;
//     // std::array<real_number, 3> Mass_boundary;

//     for (int i = 0; i < 3; i++)
//     {
//         if (r_boundary_norm[i] < params.r_cut)
//         {
//             // compute volume of boundary particle, this gives density and pressure
//             Volume_boundary[i] = 0.5 * (2.0 * params.dp + params.dp * params.dp * kappa - 2.0 * (i + 1.0) * params.dp * params.dp * kappa) * dxwall;
//             real_number ker = Wab(r_boundary_norm[i], params.H, params.Kquintic) * Volume_boundary[i];
//             W_sum += ker;
//             magnitude_tmp += vd.getProp<velocity>(boundary_key)[component] * ker;
//         }
//     }
// }

void PlaceProbes(probe_particles &probes,
                 const int k0,
                 const int kmax,
                 const Point<DIM, real_number> Corner,
                 const Point<DIM, real_number> UnitOffset)
{
    for (int k = k0; k < kmax; k++)
    {
        Point<DIM, real_number> probe_position = Corner + k * UnitOffset;
        probes.add();
        for (int xyz = 0; xyz < DIM; xyz++)
        {
            probes.getLastPos()[xyz] = probe_position.get(xyz);
        }
        probes.template getLastProp<0>() = 0.0;
    }
}