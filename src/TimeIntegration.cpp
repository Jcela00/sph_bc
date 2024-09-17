#include "TimeIntegration.hpp"

void max_velocity(particles &vd, Vcluster<> &v_cl, double &max_vel)
{
    // Calculate the maximum acceleration
    auto part = vd.getDomainIterator();

    while (part.isNext())
    {
        auto a = part.get();

        Point<DIM, double> vel(vd.getProp<velocity>(a));
        double vel2 = norm2(vel);

        if (vel2 >= max_vel)
            max_vel = vel2;

        ++part;
    }
    max_vel = std::sqrt(max_vel);

    // Vcluster<> &v_cl = create_vcluster();
    v_cl.max(max_vel);
    v_cl.execute();
}

double calc_deltaT(particles &vd, Vcluster<> &v_cl, const Parameters &params)
{
    double Maxvel = 0.0;
    max_velocity(vd, v_cl, Maxvel);

    double dt_u = 0.25 * params.H / (params.cbar + abs(Maxvel));
    double dt_visc = 0.125 * params.H * params.H / (params.nu);
    double dt_g = 0.25 * sqrt(params.H / getVectorNorm(params.gravity_vector));
    double dt = params.CFLnumber * std::min({dt_u, dt_visc, dt_g});
    // if (dt < DtMin)
    // 	dt = DtMin;

    return dt;
}

Point<DIM, double> SolidBodyAcceleration(double t, const Parameters &params)
{
    Point<DIM, double> a;
    if (params.SCENARIO == MOVING_OBSTACLE)
    {
        double period = 10.0;
        double amplitude = 0.25;
        a.get(0) = amplitude * (4.0 * M_PI * M_PI / (period * period)) * sin(2.0 * M_PI * t / period);
        a.get(1) = 0.0;
    }
    else
    {
        a = {0.0, 0.0};
    }

    return a;
}
