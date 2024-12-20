#ifndef CREATEPARTICLEGEOMETRY_HPP
#define CREATEPARTICLEGEOMETRY_HPP

#include "Calculations.hpp"
#include "Obstacle.hpp"
#include "Probes.hpp"
#include "InitializeParameters.hpp"

void CreateParticleGeometry(particles &vd,
                            std::vector<std::pair<probe_particles, int>> &vp_vec,
                            Obstacle *&obstacle_ptr,
                            Parameters &params,
                            AuxiliarParameters &auxParams);

void CreateParticleGeometryPoiseuilleTank(particles &vd, Parameters &params, AuxiliarParameters &auxParams);

void CreateParticleGeometryTaylorCouette(particles &vd,
                                         std::vector<std::pair<probe_particles, int>> &vp_vec,
                                         Obstacle *&obstacle_ptr,
                                         Parameters params,
                                         AuxiliarParameters &auxParams);

void CreateParticleGeometryStep(particles &vd,
                                std::vector<std::pair<probe_particles, int>> &vp_vec,
                                Parameters params,
                                AuxiliarParameters &auxParams);

void CreateParticleGeometryDamBreak(particles &vd,
                                    std::vector<std::pair<probe_particles, int>> &vp_vec,
                                    Parameters &params,
                                    AuxiliarParameters &auxParams);

void CreateParticleGeometryDamBreakAdj(particles &vd,
                                       std::vector<std::pair<probe_particles, int>> &vp_vec,
                                       Parameters &params,
                                       AuxiliarParameters &auxParams);

void CreateParticleGeometryCavity(particles &vd,
                                  std::vector<std::pair<probe_particles, int>> &vp_vec,
                                  Obstacle *&obstacle_ptr,
                                  Parameters &params,
                                  AuxiliarParameters &auxParams);

void CreateParticleGeometrySphere(particles &vd,
                                  std::vector<std::pair<probe_particles, int>> &vp_vec,
                                  Obstacle *&obstacle_ptr,
                                  Parameters &params,
                                  AuxiliarParameters &auxParams);

void CreateParticleGeometryFlower(particles &vd,
                                  std::vector<std::pair<probe_particles, int>> &vp_vec,
                                  Obstacle *&obstacle_ptr,
                                  Parameters &params,
                                  AuxiliarParameters &auxParams);

#endif // CREATEPARTICLEGEOMETRY_HPP
