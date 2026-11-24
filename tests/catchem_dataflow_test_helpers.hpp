#pragma once
#include <cstddef>

inline double dataflow_pattern(int field, int timestep, int column, int level = 0, int species = 0) {
    return field * 100000.0 + timestep * 10000.0 + species * 1000.0 + level * 100.0 + column;
}
