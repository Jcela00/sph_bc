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
    size_t normalization = 0;

    if constexpr (DIM == 2)
        normalization = params.Nfluid[0] * params.Nfluid[1];
    else if constexpr (DIM == 3)
        normalization = params.Nfluid[0] * params.Nfluid[1] * params.Nfluid[2];

    VxDragLift[0] = VxDragLift[0] / static_cast<real_number>(normalization);

    // compute drag and lift
    // VxDragLift[1] = VxDragLift[1] / (params.eta * VxDragLift[0]);
    // VxDragLift[2] = VxDragLift[2] / (params.eta * VxDragLift[0]);

    if (v_cl.getProcessUnitID() == 0)
    {
        if (write == 0)
        {
            avgvelstream << "t, <vx>, dragforce, liftforce" << std::endl;
        }

        avgvelstream << t << ", " << VxDragLift[0] << ", " << VxDragLift[1] << ", " << VxDragLift[2] << std::endl;
    }
}
