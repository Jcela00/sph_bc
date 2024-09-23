#ifndef CREATEPARTICLEGEOMETRY_HPP
#define CREATEPARTICLEGEOMETRY_HPP

#include "Definitions.hpp"
#include "Obstacle.hpp"
#include "Probes.hpp"
#include "VectorUtilities.hpp"
#include "InitializeParameters.hpp"

void CreateParticleGeometry(particles &vd,
                            std::vector<std::pair<probe_particles, int>> &vp_vec,
                            Vcluster<> &v_cl,
                            Obstacle *&obstacle_ptr,
                            Parameters &params);

void CreateParticleGeometryTaylorCouette(particles &vd,
                                         std::vector<std::pair<probe_particles, int>> &vp_vec,
                                         Vcluster<> &v_cl,
                                         Obstacle *&obstacle_ptr,
                                         Parameters params);

void CreateParticleGeometryStep(particles &vd,
                                std::vector<std::pair<probe_particles, int>> &vp_vec,
                                Vcluster<> &v_cl,
                                Parameters params);

#endif // CREATEPARTICLEGEOMETRY_HPP
