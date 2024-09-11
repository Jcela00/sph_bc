#include "VectorUtilities.hpp"

// Function definitions

double dotProduct(const Point<DIM, double> &v, const Point<DIM, double> &w)
{
    double result = 0.0;
    if constexpr (DIM == 2)
        result = v.get(0) * w.get(0) + v.get(1) * w.get(1);
    else if constexpr (DIM == 3)
        result = v.get(0) * w.get(0) + v.get(1) * w.get(1) + v.get(2) * w.get(2);

    return result;
}

std::array<Point<DIM, double>, DIM> dyadicProduct(const Point<DIM, double> &v, const Point<DIM, double> &w)
{
    std::array<Point<DIM, double>, DIM> dyad;
    for (size_t i = 0; i < DIM; i++)
    {
        dyad[i] = v.get(i) * w;
    }
    return dyad;
}

Point<DIM, double> matVec(const std::array<Point<DIM, double>, DIM> &m, const Point<DIM, double> &v)
{
    Point<DIM, double> res;
    for (size_t i = 0; i < DIM; i++)
    {
        res.get(i) = dotProduct(m[i], v);
    }
    return res;
}

double getVectorNorm(const Point<DIM, double> &v)
{
    return sqrt(norm2(v));
}

void normalizeVector(Point<DIM, double> &v)
{
    const double norm = getVectorNorm(v);
    if (norm > 0.0)
        v = v / norm;
}

Point<DIM, double> getPerpendicularUnit2D(const Point<DIM, double> &v)
{
    Point<DIM, double> perp;
    perp.get(0) = -v.get(1);
    perp.get(1) = v.get(0);
    return perp;
}

std::array<Point<DIM, double>, 3> getBoundaryPositions(const Point<DIM, double> &r, const Point<DIM, double> &normal, double dp)
{
    // aplies offset to a vector and returns the vectors pointing at the virtual wall particles
    // r needs to be pointing to the wall particle

    Point<DIM, double> offset_1 = -0.5 * normal * dp; // offset from wall to first particle
    Point<DIM, double> offset_2 = -1.0 * normal * dp; // offset from first to second particle, and from second to third
    std::array<Point<DIM, double>, 3> r_virtual;
    r_virtual[0] = r + offset_1;
    r_virtual[1] = r_virtual[0] + offset_2;
    r_virtual[2] = r_virtual[1] + offset_2;

    return r_virtual;
}

Point<DIM, double> crossProduct(const Point<DIM, double> &v, const Point<DIM, double> &w)
{
    Point<DIM, double> res;
    if constexpr (DIM == 2)
    {
        res.get(0) = v.get(0) * w.get(1) - v.get(1) * w.get(0);
    }
    else if constexpr (DIM == 3)
    {
        res.get(0) = v.get(1) * w.get(2) - v.get(2) * w.get(1);
        res.get(1) = v.get(2) * w.get(0) - v.get(0) * w.get(2);
        res.get(2) = v.get(0) * w.get(1) - v.get(1) * w.get(0);
    }
    return res;
}


void ApplyRotation(Point<DIM, double> &x, const double theta, const Point<DIM, double> centre)
{

    Point<DIM, double> x_rotated;
    x_rotated.get(0) = cos(theta) * (x.get(0) - centre.get(0)) - sin(theta) * (x.get(1) - centre.get(1)) + centre.get(0);
    x_rotated.get(1) = sin(theta) * (x.get(0) - centre.get(0)) + cos(theta) * (x.get(1) - centre.get(1)) + centre.get(1);

    x = x_rotated;
}