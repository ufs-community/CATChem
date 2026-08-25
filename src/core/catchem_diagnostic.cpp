// src/core/catchem_diagnostic.cpp
#include "catchem_diagnostic_manager.hpp"
#include <stdexcept>

namespace catchem {

    DiagnosticField::DiagnosticField(const std::string& name_val, const std::string& desc_val,
                                     const std::string& units_val, DiagType type_val, const std::vector<int>& dims)
        : name(name_val), description(desc_val), units(units_val), type(type_val), dimensions(dims) {
        is_gpu_target = !std::is_same_v<HostSpace, DeviceSpace>;

        if (type == DiagType::FIELD_2D) {
            if (dims.size() != 2)
                throw std::invalid_argument("2D field requires 2 dimensions");
            if (is_gpu_target)
                device_view_2d = View2D("dev_" + name, dims[0], dims[1]);
            host_view_2d = HostView2D("host_" + name, dims[0], dims[1]);
            if (!is_gpu_target)
                device_view_2d = host_view_2d;
        } else if (type == DiagType::FIELD_3D) {
            if (dims.size() != 3)
                throw std::invalid_argument("3D field requires 3 dimensions");
            if (is_gpu_target)
                device_view_3d = View3D("dev_" + name, dims[0], dims[1], dims[2]);
            host_view_3d = HostView3D("host_" + name, dims[0], dims[1], dims[2]);
            if (!is_gpu_target)
                device_view_3d = host_view_3d;
        } else {
            throw std::invalid_argument("Unsupported DiagType");
        }
    }

    void DiagnosticField::sync_to_host() {
        if constexpr (!std::is_same_v<HostSpace, DeviceSpace>) {
            if (type == DiagType::FIELD_2D) {
                Kokkos::deep_copy(host_view_2d, device_view_2d);
            } else if (type == DiagType::FIELD_3D) {
                Kokkos::deep_copy(host_view_3d, device_view_3d);
            }
        }
    }

    void DiagnosticField::sync_to_device() {
        if constexpr (!std::is_same_v<HostSpace, DeviceSpace>) {
            if (type == DiagType::FIELD_2D) {
                Kokkos::deep_copy(device_view_2d, host_view_2d);
            } else if (type == DiagType::FIELD_3D) {
                Kokkos::deep_copy(device_view_3d, host_view_3d);
            }
        }
    }

    void DiagnosticField::reset() {
        if (type == DiagType::FIELD_2D) {
            Kokkos::deep_copy(device_view_2d, 0.0);
            Kokkos::deep_copy(host_view_2d, 0.0);
        } else if (type == DiagType::FIELD_3D) {
            Kokkos::deep_copy(device_view_3d, 0.0);
            Kokkos::deep_copy(host_view_3d, 0.0);
        }
    }

    void DiagnosticManager::register_field(const std::string& name, const std::string& desc, const std::string& units,
                                           DiagType type, const std::vector<int>& dims) {
        fields[name] = std::make_shared<DiagnosticField>(name, desc, units, type, dims);
    }

    bool DiagnosticManager::has_field(const std::string& name) const {
        return fields.find(name) != fields.end();
    }

    std::shared_ptr<DiagnosticField> DiagnosticManager::get_field(const std::string& name) {
        if (!has_field(name))
            throw std::invalid_argument("Field not found: " + name);
        return fields.at(name);
    }

    Kokkos::View<double**, Kokkos::LayoutLeft, Kokkos::DefaultExecutionSpace::memory_space>
    DiagnosticManager::get_device_view_2d(const std::string& name) {
        auto field = get_field(name);
        if (field->type != DiagType::FIELD_2D)
            throw std::invalid_argument("Field is not 2D: " + name);
        return field->device_view_2d;
    }

    Kokkos::View<double***, Kokkos::LayoutLeft, Kokkos::DefaultExecutionSpace::memory_space>
    DiagnosticManager::get_device_view_3d(const std::string& name) {
        auto field = get_field(name);
        if (field->type != DiagType::FIELD_3D)
            throw std::invalid_argument("Field is not 3D: " + name);
        return field->device_view_3d;
    }

    void* DiagnosticManager::get_host_pointer(const std::string& name) {
        auto field = get_field(name);
        if (field->type == DiagType::FIELD_2D) {
            return static_cast<void*>(field->host_view_2d.data());
        } else if (field->type == DiagType::FIELD_3D) {
            return static_cast<void*>(field->host_view_3d.data());
        }
        return nullptr;
    }

    void DiagnosticManager::sync_to_host() {
        for (auto& [key, field] : fields) {
            field->sync_to_host();
        }
    }

    void DiagnosticManager::sync_to_device() {
        for (auto& [key, field] : fields) {
            field->sync_to_device();
        }
    }

    void DiagnosticManager::reset_all() {
        for (auto& [key, field] : fields) {
            field->reset();
        }
    }

    std::vector<std::string> DiagnosticManager::get_registered_names() const {
        std::vector<std::string> names;
        for (const auto& [name, field] : fields) {
            names.push_back(name);
        }
        return names;
    }

} // namespace catchem
