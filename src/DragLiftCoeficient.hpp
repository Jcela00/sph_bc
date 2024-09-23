#ifndef DRAG_LIFT_COEFICIENT_HPP
#define DRAG_LIFT_COEFICIENT_HPP

#include "Definitions.hpp"
#include "VectorUtilities.hpp"

void CalcDragLift(particles &vd,
                  Vcluster<> &v_cl,
                  real_number t,
                  std::ofstream &avgvelstream,
                  real_number obstacle_force_x,
                  real_number obstacle_force_y,
                  const Parameters &params,
                  size_t write);

#endif // DRAG_LIFT_COEFICIENT_HPP
