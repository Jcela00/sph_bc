#ifndef EQSTATE_CU
#define EQSTATE_CU

#include "Definitions.hpp"

template <typename vd_type>
__global__ void EqStateGPU(vd_type vd, real_number rho0, real_number B, real_number gamma, real_number xi)
{
    auto a = GET_PARTICLE(vd);

    if (vd.template getProp<vd0_type>(a) == FLUID)
    {
        real_number rho_a = vd.template getProp<vd1_rho>(a);
        real_number rho_frac = rho_a / rho0;

        vd.template getProp<vd2_pressure>(a) = B * (std::pow(rho_frac, gamma) - 1.0) + xi;
    }
}

inline void EqState(particles &vd, const real_number rho0, const real_number B, const real_number gamma, const real_number xi)
{
    auto it = vd.getDomainIteratorGPU();
    CUDA_LAUNCH(EqStateGPU, it, vd.toKernel(), rho0, B, gamma, xi);
}

#endif