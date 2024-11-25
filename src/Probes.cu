#include "Probes.hpp"

void PlaceProbes(probe_particles &probes,
                 const int k0,
                 const int kmax,
                 const Point<DIM, real_number> Corner,
                 const Point<DIM, real_number> UnitOffset,
                 Obstacle *obstacle_ptr,
                 const std::vector<int> FixedProbeIndices,
                 const std::vector<real_number> FixedProbeValues)
{
    for (unsigned int k = k0; k < kmax; k++)
    {
        Point<DIM, real_number> probe_position = Corner + k * UnitOffset;
        probes.add();
        for (int xyz = 0; xyz < DIM; xyz++)
        {
            probes.getLastPos()[xyz] = probe_position.get(xyz);
        }

        auto it = std::find(FixedProbeIndices.begin(), FixedProbeIndices.end(), k);
        if (it != FixedProbeIndices.end()) // if k is in FixedProbeIndices
        {
            // k is a fixed probe
            probes.template getLastProp<probe0_type>() = FIXED_PROBE;

            // Get the index of k in FixedProbeIndices
            size_t index = std::distance(FixedProbeIndices.begin(), it);

            // Assign the corresponding FixedProbeValues[index]
            probes.template getLastProp<probe1_quantity>() = FixedProbeValues[index];
        }
        else if (obstacle_ptr->isInsidePlusTol(probe_position)) // if probe is inside the obstacle it gets a 0 value
        {
            probes.template getLastProp<probe0_type>() = FIXED_PROBE;
            probes.template getLastProp<probe1_quantity>() = 0.0;
        }
        else
        {
            probes.template getLastProp<probe0_type>() = VARIABLE_PROBE;
            probes.template getLastProp<probe1_quantity>() = 0.0;
        }
    }
}