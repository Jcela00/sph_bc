#include "DragLiftCoeficient.hpp"

void CalcDragLift(real_number t,
                  std::ofstream &avgvelstream,
                  Point<3, real_number> &VxDragLift,
                  const Parameters &params,
                  size_t write)
{
    // VxDragLif contains vx, dragForce, liftForce
    Vcluster<> &v_cl = create_vcluster();

    // normalization factor, 1/normalization is equal to dx^2/Channel_area
    size_t normalization = params.Nfluid[0] * params.Nfluid[1];

    VxDragLift[0] = VxDragLift[0] / (real_number)normalization;

    // compute drag and lift
    VxDragLift[1] = VxDragLift[1] / (params.eta * VxDragLift[0]);
    VxDragLift[2] = VxDragLift[2] / (params.eta * VxDragLift[0]);
    real_number Reynolds = VxDragLift[0] * params.ObstacleBase / params.nu;
    if (v_cl.getProcessUnitID() == 0)
    {
        if (write == 0)
        {
            avgvelstream << "t, <vx>, drag, lift, Re" << std::endl;
        }

        avgvelstream << t << ", " << VxDragLift[0] << ", " << VxDragLift[1] << ", " << VxDragLift[2] << ", " << Reynolds << std::endl;
    }
}
