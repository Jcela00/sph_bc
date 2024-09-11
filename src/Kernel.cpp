#include "Kernel.hpp"

double Wab(double r, const double H, const double Kquintic)
{
    const double q = r / H;
    const double tmp3 = (3.0 - q) * (3.0 - q) * (3.0 - q) * (3.0 - q) * (3.0 - q);
    const double tmp2 = (2.0 - q) * (2.0 - q) * (2.0 - q) * (2.0 - q) * (2.0 - q);
    const double tmp1 = (1.0 - q) * (1.0 - q) * (1.0 - q) * (1.0 - q) * (1.0 - q);

    if (q < 0.0)
        return 0.0;
    else if (q < 1.0)
        return Kquintic * (tmp3 - 6.0 * tmp2 + 15.0 * tmp1);
    else if (q < 2.0)
        return Kquintic * (tmp3 - 6.0 * tmp2);
    else if (q < 3.0)
        return Kquintic * tmp3;
    else
        return 0.0;
}

Point<DIM, double> DWab(const Point<DIM, double> &dx, const double r, const double H, const double Kquintic)
{
    Point<DIM, double> DW;
    const double q = r / H;

    const double tmp3 = (3.0 - q) * (3.0 - q) * (3.0 - q) * (3.0 - q);
    const double tmp2 = (2.0 - q) * (2.0 - q) * (2.0 - q) * (2.0 - q);
    const double tmp1 = (1.0 - q) * (1.0 - q) * (1.0 - q) * (1.0 - q);

    const double c1 = -5.0 * Kquintic / H;
    double factor;

    if (q < 1.0)
    {
        if (q > 0.0)
            factor = c1 * (tmp3 - 6.0 * tmp2 + 15.0 * tmp1) / r;
        else
            factor = 0.0;
    }
    else if (q < 2.0)
        factor = c1 * (tmp3 - 6.0 * tmp2) / r;
    else if (q < 3.0)
        factor = c1 * tmp3 / r;
    else
        factor = 0.0;

    DW = factor * dx;

    return DW;
}