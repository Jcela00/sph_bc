#ifndef INTERACTIONS_H
#define INTERACTIONS_H

#include "Definitions.hpp"
#include "VectorUtilities.hpp"
#include "Physics.hpp"
#include "Kernel.hpp"

void interact_fluid_boundary_new(particles &vd,
                                 vect_dist_key_dx fluid_key,
                                 const real_number &massf,
                                 const real_number &rhof,
                                 const real_number &Pf,
                                 const Point<DIM, real_number> &xf,
                                 const Point<DIM, real_number> &vf,
                                 const Point<DIM, real_number> &xw,
                                 const Point<DIM, real_number> &r_wall_to_fluid,
                                 unsigned long &boundary_key,
                                 bool accumulate_force,
                                 real_number &obstacle_force_x,
                                 real_number &obstacle_force_y,
                                 const Parameters &params);

void interact_fluid_boundary_old(particles &vd,
                                 vect_dist_key_dx &fluid_key,
                                 const real_number &massf,
                                 const real_number &rhof,
                                 const real_number &Pf,
                                 const Point<DIM, real_number> &xf,
                                 const Point<DIM, real_number> &vf,
                                 const Point<DIM, real_number> &xb,
                                 const Point<DIM, real_number> &r_ab,
                                 const real_number &r2,
                                 unsigned long &boundary_key,
                                 bool accumulate_force,
                                 real_number &obstacle_force_x,
                                 real_number &obstacle_force_y,
                                 const Parameters &params);

void interact_fluid_fluid(particles &vd,
                          const vect_dist_key_dx &a_key,
                          const real_number &massa,
                          const real_number &rhoa,
                          const real_number &Pa,
                          const Point<DIM, real_number> &xa,
                          const Point<DIM, real_number> &va,
                          const Point<DIM, real_number> &xb,
                          const Point<DIM, real_number> &r_ab,
                          const real_number &r2,
                          const unsigned long &b_key,
                          const Parameters &params);

void interact_fluid_boundary_new_regularize(particles &vd,
                                            vect_dist_key_dx fluid_key,
                                            const real_number &rhof,
                                            const Point<DIM, real_number> &r_wall_to_fluid,
                                            unsigned long &boundary_key,
                                            const Parameters &params);

void interact_fluid_fluid_regularize(particles &vd,
                                     const vect_dist_key_dx &a_key,
                                     const real_number &massa,
                                     const real_number &rhoa,
                                     const Point<DIM, real_number> &r_ab,
                                     const real_number &r2,
                                     const unsigned long &b_key,
                                     const Parameters &params);

#endif // INTERACTIONS_H
