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

    double vx = 0.0, v = 0.0;

    size_t counter = 0;
    auto it = vd.getDomainIterator();
    while (it.isNext())
    {
        auto key = it.get();

        if (vd.template getProp<type>(key) == FLUID)
        {
            vx += vd.template getProp<velocity>(key)[0];
            // vy += vd.template getProp<velocity>(key)[1];
            // vtx += vd.template getProp<v_transport>(key)[0];
            // vty += vd.template getProp<v_transport>(key)[1];
            v += getVectorNorm(vd.template getProp<velocity>(key));
            // vt += getVectorNorm(vd.template getProp<v_transport>(key));

            counter++;
        }

        ++it;
    }

    v_cl.sum(vx);
    // v_cl.sum(vy);
    v_cl.sum(v);
    // v_cl.sum(vtx);
    // v_cl.sum(vty);
    // v_cl.sum(vt);
    v_cl.sum(counter);
    v_cl.execute();
    vx = vx / (double)counter;
    // vy = vy / (double)counter;
    v = v / (double)counter;
    // vtx = vtx / (double)counter;
    // vty = vty / (double)counter;
    // vt = vt / (double)counter;
    // double drag_alt = obstacle_force_x / (params.eta * vx);

    double drag = obstacle_force_x / (params.eta * v);
    double drag_vx = obstacle_force_x / (params.eta * vx);
    double drag_12 = obstacle_force_x / (params.eta * 1.2e-4);

    // double lift_alt = obstacle_force_y / (params.eta * v);
    // obstacle_force_y = obstacle_force_y / (params.eta * vx);

    if (v_cl.getProcessUnitID() == 0)
    {
        if (write == 0)
        {
            avgvelstream << "t, <v>,<vx>, drag, drag_vx, drag_12" << std::endl;
        }

        avgvelstream << t << ", " << v << ", " << vx << ", " << drag << ", " << drag_vx << ", " << drag_12 << std::endl;
    }
}
