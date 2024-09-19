#include "DragLiftCoeficient.hpp"

void CalcDragLift(particles &vd,
                  Vcluster<> &v_cl,
                  double t,
                  std::ofstream &avgvelstream,
                  double obstacle_force_x,
                  double obstacle_force_y,
                  const Parameters &params,
                  size_t write)
{

    double vx = 0.0;

    auto it = vd.getDomainIterator();
    while (it.isNext())
    {
        auto key = it.get();

        if (vd.template getProp<type>(key) == FLUID)
        {
            vx += vd.template getProp<velocity>(key)[0];
        }

        ++it;
    }

    // normalization factor, 1/normalization is equal to dx^2/Channel area
    size_t normalization = params.Nfluid[0] * params.Nfluid[1];

    // add across processors
    v_cl.sum(vx);
    v_cl.execute();

    vx = vx / (double)normalization;

    // compute drag and lift
    double drag = obstacle_force_x / (params.eta * vx);
    double lift = obstacle_force_y / (params.eta * vx);

    if (v_cl.getProcessUnitID() == 0)
    {
        if (write == 0)
        {
            avgvelstream << "t, <vx>, drag, lift" << std::endl;
        }

        avgvelstream << t << ", " << vx << ", " << drag << ", " << lift << std::endl;
    }
}
