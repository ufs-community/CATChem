// src/core/catchem_grid_manager.hpp
#pragma once
#include "catchem_interop_field.hpp"
#include <memory>

namespace catchem {

    struct GridGeometry {
        int nx = 1;
        int ny = 1;
        int nz = 1;

        std::shared_ptr<InteropField<double, 2>> lat;
        std::shared_ptr<InteropField<double, 2>> lon;
        std::shared_ptr<InteropField<double, 2>> grid_area;
        // dz and z_levels can be added as 1D fields if needed globally
    };

    class GridManager {
    public:
        GridGeometry geometry;
        bool is_initialized = false;

        GridManager(int nx, int ny, int nz);

        // Bindings to support Fortran Interop arrays if allocated externally
        void bind_lat(double* ptr);
        void bind_lon(double* ptr);
        void bind_area(double* ptr);
    };

} // namespace catchem
