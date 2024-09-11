#include "DragLiftCoeficient.hpp"

void CalcDragLift(particles &vd,
                  Vcluster<> &v_cl,
                  double t,
                  std::ofstream &avgvelstream,
                  double obstacle_force_x,
                  double obstacle_force_y,
                  const Parameters &params)
{
    auto it = vd.getDomainIterator();
    double vx{0.0}, vy{0.0}, v{0.0};
    int counter = 0;
    while (it.isNext())
    {
        auto key = it.get();

        if (vd.template getProp<type>(key) == FLUID)
        {
            vx += vd.template getProp<velocity>(key)[0];
            vy += vd.template getProp<velocity>(key)[1];
            v += getVectorNorm(vd.template getProp<velocity>(key));
            counter++;
        }

        ++it;
    }

    v_cl.sum(vx);
    v_cl.sum(vy);
    v_cl.sum(v);
    v_cl.sum(counter);
    v_cl.execute();
    vx = vx / counter;
    vy = vy / counter;
    v = v / counter;
    double drag_alt = obstacle_force_x / (params.eta * vx);
    obstacle_force_x = obstacle_force_x / (params.eta * v);
    obstacle_force_y = obstacle_force_y / (params.eta * v);

    if (v_cl.getProcessUnitID() == 0)
    {
        avgvelstream << "t = " << t << ", vx = " << vx << ", vy = " << vy << ", v = " << v << ", drag = " << obstacle_force_x << ", drag_alt = " << drag_alt << ", lift = " << obstacle_force_y << std::endl;
    }
}
