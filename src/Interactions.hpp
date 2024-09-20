#ifndef INTERACTIONS_H
#define INTERACTIONS_H

#include "Definitions.hpp"
#include "VectorUtilities.hpp"
#include "Physics.hpp"
#include "Kernel.hpp"

void interact_fluid_boundary_new(particles &vd,
                                 vect_dist_key_dx fluid_key,
                                 const double &massf,
                                 const double &rhof,
                                 const double &Pf,
                                 const Point<DIM, double> &xf,
                                 const Point<DIM, double> &vf,
                                 const Point<DIM, double> &xw,
                                 const Point<DIM, double> &r_wall_to_fluid,
                                 unsigned long &boundary_key,
                                 bool accumulate_force,
                                 double &obstacle_force_x,
                                 double &obstacle_force_y,
                                 const Parameters &params);

void interact_fluid_boundary_old(particles &vd,
                                 vect_dist_key_dx &fluid_key,
                                 const double &massf,
                                 const double &rhof,
                                 const double &Pf,
                                 const Point<DIM, double> &xf,
                                 const Point<DIM, double> &vf,
                                 const Point<DIM, double> &xb,
                                 const Point<DIM, double> &r_ab,
                                 const double &r2,
                                 unsigned long &boundary_key,
                                 bool accumulate_force,
                                 double &obstacle_force_x,
                                 double &obstacle_force_y,
                                 const Parameters &params);

void interact_fluid_fluid(particles &vd,
                          const vect_dist_key_dx &a_key,
                          const double &massa,
                          const double &rhoa,
                          const double &Pa,
                          const Point<DIM, double> &xa,
                          const Point<DIM, double> &va,
                          const Point<DIM, double> &xb,
                          const Point<DIM, double> &r_ab,
                          const double &r2,
                          const unsigned long &b_key,
                          const Parameters &params);

void interact_fluid_boundary_new_regularize(particles &vd,
                                            vect_dist_key_dx fluid_key,
                                            const double &rhof,
                                            const Point<DIM, double> &r_wall_to_fluid,
                                            unsigned long &boundary_key,
                                            const Parameters &params);

void interact_fluid_fluid_regularize(particles &vd,
                                     const vect_dist_key_dx &a_key,
                                     const double &massa,
                                     const double &rhoa,
                                     const Point<DIM, double> &r_ab,
                                     const double &r2,
                                     const unsigned long &b_key,
                                     const Parameters &params);

#endif // INTERACTIONS_H
