// src/core/catchem_diagnostic_manager.hpp
#pragma once
#include "catchem_diagnostic.hpp"
#include <memory>
#include <string>
#include <unordered_map>

namespace catchem {

    class DiagnosticManager {
    private:
        std::unordered_map<std::string, std::shared_ptr<DiagnosticField>> fields;

    public:
        DiagnosticManager() = default;

        void register_field(const std::string& name, const std::string& desc, const std::string& units, DiagType type,
                            const std::vector<int>& dims);

        bool has_field(const std::string& name) const;
        std::shared_ptr<DiagnosticField> get_field(const std::string& name);

        Kokkos::View<double**, Kokkos::LayoutLeft, Kokkos::DefaultExecutionSpace::memory_space>
        get_device_view_2d(const std::string& name);

        Kokkos::View<double***, Kokkos::LayoutLeft, Kokkos::DefaultExecutionSpace::memory_space>
        get_device_view_3d(const std::string& name);

        void* get_host_pointer(const std::string& name);
        void sync_to_host();
        void sync_to_device();
        void reset_all();
        std::vector<std::string> get_registered_names() const;
    };

} // namespace catchem
