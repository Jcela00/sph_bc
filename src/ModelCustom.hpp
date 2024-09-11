#ifndef MODELCUSTOM_H
#define MODELCUSTOM_H

#include "Definitions.hpp"

struct ModelCustom
{
    template <typename Decomposition, typename vector>
    inline void addComputation(Decomposition &dec, vector &vd, size_t v, size_t p)
    {
        if (vd.template getProp<type>(p) == FLUID)
            dec.addComputationCost(v, 4);
        else
            dec.addComputationCost(v, 4);
    }

    template <typename Decomposition>
    inline void applyModel(Decomposition &dec, size_t v)
    {
        dec.setSubSubDomainComputationCost(v, dec.getSubSubDomainComputationCost(v) * dec.getSubSubDomainComputationCost(v));
    }

    double distributionTol()
    {
        return 1.01;
    }
};

#endif // MODELCUSTOM_H