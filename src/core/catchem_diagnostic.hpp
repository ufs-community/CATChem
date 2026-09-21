// src/core/catchem_diagnostic.hpp
#pragma once
#include "catchem_field_contract.hpp"
#include "catchem_interop_field.hpp"
#include "catchem_kokkos_compat.hpp"
#include <string>
#include <vector>

namespace catchem {

    enum class DiagType { SCALAR, FIELD_1D, FIELD_2D, FIELD_3D };
    enum class DiagnosticPolicy { Instantaneous, TimestepAccumulated, Persistent };

    class DiagnosticField {
    public:
        using Mdspan2D = typename MdspanTypeHelper<double, 2>::type;
        using Mdspan3D = typename MdspanTypeHelper<double, 3>::type;

        std::string name;
        std::string description;
        std::string units;
        DiagType type;
        std::vector<int> dimensions;
        std::vector<SemanticAxis> axes;
        /// Human-readable label for each position of the field's single packed
        /// (Species or Category) dimension; empty when no dimension is packed.
        /// Index i is the 0-based slot the scheme writes, NOT a global species index.
        std::vector<std::string> unpack_labels;
        std::size_t generation = 0;
        std::size_t registration_generation = 0;
        AvailabilityState availability = AvailabilityState::Unavailable;
        bool generation_failed = false;
        LatestWriter latest_writer = LatestWriter::Uninitialized;
        DiagnosticPolicy reset_policy = DiagnosticPolicy::Instantaneous;
        double reset_value = 0.0;

#ifdef CATCHEM_ENABLE_KOKKOS
        using HostSpace = Kokkos::HostSpace;
        using DeviceSpace = Kokkos::DefaultExecutionSpace::memory_space;

        using View2D = Kokkos::View<double**, Kokkos::LayoutLeft, DeviceSpace>;
        using HostView2D = Kokkos::View<double**, Kokkos::LayoutLeft, HostSpace>;

        using View3D = Kokkos::View<double***, Kokkos::LayoutLeft, DeviceSpace>;
        using HostView3D = Kokkos::View<double***, Kokkos::LayoutLeft, HostSpace>;

        View2D device_view_2d;
        HostView2D host_view_2d;

        View3D device_view_3d;
        HostView3D host_view_3d;
#else
        // Host-only builds own their storage directly; mdspans over `storage`
        // serve as both the "host" and "device" views.
        std::vector<double> storage;
#endif

        bool is_gpu_target;

        DiagnosticField(const std::string& name_val, const std::string& desc_val, const std::string& units_val,
                        DiagType type_val, const std::vector<int>& dims,
                        DiagnosticPolicy policy = DiagnosticPolicy::Instantaneous, double reset = 0.0,
                        std::vector<SemanticAxis> semantic_axes = {}, std::vector<std::string> labels = {});

        void sync_to_host();
        void sync_to_device();
        void reset();
        void advance_generation(std::size_t value);
        void mark_host_modified() { latest_writer = LatestWriter::HostCurrent; }
        void mark_device_modified() { latest_writer = LatestWriter::DeviceCurrent; }
    };

} // namespace catchem
