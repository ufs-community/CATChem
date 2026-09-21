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
        void register_field_contract(const std::string& name, const std::string& desc, const std::string& units,
                                     DiagType type, const std::vector<int>& dims, DiagnosticPolicy policy,
                                     double reset_value, const std::vector<SemanticAxis>& axes,
                                     const std::vector<std::string>& unpack_labels = {}, bool strict_labels = true);

        bool has_field(const std::string& name) const;
        std::shared_ptr<DiagnosticField> get_field(const std::string& name);
        /// Semantic axis per dimension (length == rank). Throws on unknown name.
        const std::vector<SemanticAxis>& get_axes(const std::string& name) const;
        /// Labels for the single packed (Species/Category) dimension; empty when the
        /// field has none. Throws on unknown name.
        const std::vector<std::string>& get_unpack_labels(const std::string& name) const;

#ifdef CATCHEM_ENABLE_KOKKOS
        Kokkos::View<double**, Kokkos::LayoutLeft, Kokkos::DefaultExecutionSpace::memory_space>
        get_device_view_2d(const std::string& name);

        Kokkos::View<double***, Kokkos::LayoutLeft, Kokkos::DefaultExecutionSpace::memory_space>
        get_device_view_3d(const std::string& name);
#else
        DiagnosticField::Mdspan2D get_device_view_2d(const std::string& name);
        DiagnosticField::Mdspan3D get_device_view_3d(const std::string& name);
#endif

        const void* get_host_read_pointer(const std::string& name);
        void* get_host_write_pointer(const std::string& name);
        void* get_host_pointer(const std::string& name) { return const_cast<void*>(get_host_read_pointer(name)); }
        void mark_host_modified(const std::string& name);
        void mark_device_modified(const std::string& name);
        void sync_to_host();
        void sync_to_device();
        void reset_all();
        void begin_timestep();
        void mark_generation_failed();
        std::size_t generation() const noexcept { return generation_; }
        std::vector<std::string> get_registered_names() const;

    private:
        std::size_t generation_ = 0;
        /// Insertion order of field names. `fields` is an unordered_map, so this
        /// side-list is what makes get_registered_names() (and therefore the
        /// axes-driven diagnostic writer's variable order) deterministic.
        std::vector<std::string> registration_order_;
    };

} // namespace catchem
