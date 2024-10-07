#ifndef DRAG_LIFT_COEFICIENT_HPP
#define DRAG_LIFT_COEFICIENT_HPP

#include "Definitions.hpp"
#include "VectorUtilities.hpp"

void CalcDragLift(real_number t,
                  std::ofstream &avgvelstream,
                  Point<3, real_number> &VxDragLift,
                  const Parameters &params,
                  size_t write);

#endif // DRAG_LIFT_COEFICIENT_HPP
