// src/core/catchem_diagnostic.hpp
#pragma once
#include <Kokkos_Core.hpp>
#include <string>
#include <vector>

namespace catchem {

    enum class DiagType { SCALAR, FIELD_1D, FIELD_2D, FIELD_3D };

    class DiagnosticField {
    public:
        using HostSpace = Kokkos::HostSpace;
        using DeviceSpace = Kokkos::DefaultExecutionSpace::memory_space;

        using View2D = Kokkos::View<double**, Kokkos::LayoutLeft, DeviceSpace>;
        using HostView2D = Kokkos::View<double**, Kokkos::LayoutLeft, HostSpace>;

        using View3D = Kokkos::View<double***, Kokkos::LayoutLeft, DeviceSpace>;
        using HostView3D = Kokkos::View<double***, Kokkos::LayoutLeft, HostSpace>;

        std::string name;
        std::string description;
        std::string units;
        DiagType type;
        std::vector<int> dimensions;

        View2D device_view_2d;
        HostView2D host_view_2d;

        View3D device_view_3d;
        HostView3D host_view_3d;

        bool is_gpu_target;

        DiagnosticField(const std::string& name_val, const std::string& desc_val, const std::string& units_val,
                        DiagType type_val, const std::vector<int>& dims);

        void sync_to_host();
        void sync_to_device();
        void reset();
    };

} // namespace catchem
