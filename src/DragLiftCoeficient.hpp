#ifndef DRAG_LIFT_COEFICIENT_HPP
#define DRAG_LIFT_COEFICIENT_HPP

#include "Definitions.hpp"
#include "VectorUtilities.hpp"

void CalcDragLift(particles &vd,
                  Vcluster<> &v_cl,
                  double t,
                  std::ofstream &avgvelstream,
                  double obstacle_force_x,
                  double obstacle_force_y,
                  const Parameters &params,
                  size_t write);

#endif // DRAG_LIFT_COEFICIENT_HPP
