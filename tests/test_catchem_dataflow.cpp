#include "catchem_api.hpp"
#include "catchem_core.hpp"
#include "catchem_dataflow_test_helpers.hpp"
#include "catchem_diagnostic_manager.hpp"
#include "catchem_species_metadata.hpp"
#include "catchem_state_manager.hpp"
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <vector>

// Unlike assert(), CHECK always evaluates its expression: several calls below
// populate pointers consumed on later lines, and NDEBUG builds must not skip
// the side effects (or the checks themselves).
#define CHECK(expr)                                                                                                    \
    do {                                                                                                               \
        if (!(expr)) {                                                                                                 \
            std::cerr << "CHECK failed: " << #expr << " at " << __FILE__ << ':' << __LINE__ << std::endl;              \
            std::abort();                                                                                              \
        }                                                                                                              \
    } while (0)

int main() {
#ifdef CATCHEM_ENABLE_KOKKOS
    Kokkos::initialize();
    std::cout << "Kokkos execution space: " << Kokkos::DefaultExecutionSpace::name()
              << ", concurrency=" << Kokkos::DefaultExecutionSpace().concurrency() << '\n';
#endif
    constexpr int nc = 3, nl = 2, ns = 3;
    {
        catchem::MechanismDefinition first;
        first.identity = "first-order";
        first.species = {{}, {}, {}};
        first.species[0].short_name = "O3";
        first.species[1].short_name = "NO2";
        first.species[2].short_name = "CO";
        first.rebuild_index();
        CHECK(first.index_of("o3") == 0);
        CHECK(first.index_of("CO") == 2);

        catchem::MechanismDefinition reordered = first;
        std::swap(reordered.species[0], reordered.species[2]);
        reordered.identity = "reordered";
        reordered.rebuild_index();
        CHECK(reordered.index_of("O3") == 2);
        CHECK(reordered.index_of("co") == 0);
    }
    {
        catchem::ConfigManager duplicate_config;
        catchem::SpeciesConfig first;
        first.name = "O3";
        catchem::SpeciesConfig second;
        second.name = "o3";
        duplicate_config.data.species = {first, second};
        catchem::ChemState chemistry_state;
        bool duplicate_rejected = false;
        try {
            chemistry_state.load_from_config_manager(duplicate_config);
        } catch (const std::invalid_argument&) {
            duplicate_rejected = true;
        }
        CHECK(duplicate_rejected);
    }

    void* core = catchem_core_create(nc, nl, ns);
    CHECK(core);
    void* state = catchem_core_get_state_manager(core);
    std::vector<double> field(nc * nl, 1.0);

    CHECK(catchem_state_bind_met_3d_checked(nullptr, "T", field.data(), nc, nl, 1) == CATCHEM_NULL_ARGUMENT);
    CHECK(catchem_state_bind_met_3d_checked(state, "T", field.data(), nc, nl + 1, 1) == CATCHEM_EXTENT_MISMATCH);
    CHECK(catchem_state_bind_met_3d_checked(state, "T", field.data(), nc, nl, 1) == CATCHEM_SUCCESS);
    auto* state_object = static_cast<catchem::StateManager*>(state);
    std::vector<double> latitude(nc, 45.0);
    state_object->bind_met_field_2d("LAT", latitude.data());
    CHECK(state_object->meteorology().LAT->contract.persistence == catchem::PersistencePolicy::Persistent);
    CHECK(state_object->meteorology().LAT->contract.canonical_name == "LAT");
    CHECK(state_object->meteorology().T->contract.units == "K");
    state_object->begin_import_generation();
    CHECK(state_object->find_2d_ptr({"LAT"}) == latitude.data());
    CHECK(state_object->find_3d_ptr({"T"}) == nullptr);
    CHECK(catchem_state_bind_met_3d_checked(state, "T", field.data(), nc, nl, 1) == CATCHEM_SUCCESS);

    std::vector<double> chemistry(nc * nl * ns);
    for (int s = 0; s < ns; ++s)
        for (int l = 0; l < nl; ++l)
            for (int c = 0; c < nc; ++c)
                chemistry[c + nc * (l + nl * s)] = dataflow_pattern(1, 1, c, l, s);
    CHECK(catchem_state_bind_unified_chemistry_checked(nullptr, chemistry.data(), nc, nl, ns) == CATCHEM_NULL_ARGUMENT);
    CHECK(catchem_state_bind_unified_chemistry_checked(state, chemistry.data(), nc + 1, nl, ns) ==
          CATCHEM_EXTENT_MISMATCH);
    CHECK(catchem_state_bind_unified_chemistry_checked(state, chemistry.data(), nc, nl, ns) == CATCHEM_SUCCESS);
    state_object->chemistry().species_name_to_index["O3"] = 1;
    CHECK(catchem_state_get_species_index(state, "o3") == 2);
    CHECK(catchem_state_get_species_index(state, "O3") == 2);
    state_object->chemistry().conc->sync_to_device();
    CHECK(state_object->chemistry().conc->latest_writer == catchem::LatestWriter::Synchronized);
    chemistry[0] = 99.0;
    CHECK(catchem_state_mark_chem_host_modified(state) == CATCHEM_SUCCESS);
    CHECK(state_object->chemistry().conc->latest_writer == catchem::LatestWriter::HostCurrent);
    state_object->chemistry().conc->sync_to_device();
    CHECK(state_object->chemistry().conc->view()(0, 0, 0) == 99.0);
    double* slab = nullptr;
#ifdef CATCHEM_ENABLE_KOKKOS
    // The checked concentration pointer is the NUOPC/read boundary.  A
    // device-side writer must be visible there without requiring the caller
    // to know CATChem's internal mirror state.
    auto device_concentration = state_object->chemistry().conc->device_write();
    Kokkos::deep_copy(device_concentration, 314.0);
#endif
    CHECK(catchem_state_get_species_conc_pointer_checked(state, 1, nc, nl, &slab) == CATCHEM_SUCCESS);
    CHECK(slab == chemistry.data());
#ifdef CATCHEM_ENABLE_KOKKOS
    CHECK(slab[0] == 314.0);
#endif
    CHECK(catchem_state_get_species_conc_pointer_checked(state, ns, nc, nl, &slab) == CATCHEM_SUCCESS);
    CHECK(slab == chemistry.data() + nc * nl * (ns - 1));
    CHECK(catchem_state_get_species_conc_pointer_checked(state, 0, nc, nl, &slab) == CATCHEM_INVALID_INDEX);
    CHECK(catchem_state_get_species_conc_pointer_checked(state, 1, nc + 1, nl, &slab) == CATCHEM_EXTENT_MISMATCH);

    catchem_diag_register(core, "shape", "shape test", "1", 3, nc, nl, ns);
    int dims[3] = {nc, nl, ns};
    void* diag = nullptr;
    CHECK(catchem_diag_get_pointer_checked(core, "shape", 3, dims, &diag) == CATCHEM_SUCCESS);
    CHECK(diag != nullptr);
    int wrong[3] = {nc, nl, ns + 1};
    CHECK(catchem_diag_get_pointer_checked(core, "shape", 3, wrong, &diag) == CATCHEM_EXTENT_MISMATCH);
    CHECK(catchem_diag_get_pointer_checked(core, "missing", 3, dims, &diag) == CATCHEM_MISSING_FIELD);
    CHECK(catchem_diag_register_checked(core, "shape", "shape test", "1", 3, nc, nl, ns + 1) ==
          CATCHEM_EXTENT_MISMATCH);

    auto* values = static_cast<double*>(catchem_diag_get_pointer(core, "shape"));
    values[0] = 42.0;
    CHECK(catchem_diag_mark_host_modified(core, "shape") == CATCHEM_SUCCESS);
    catchem_diag_sync_to_host(core); // must not replace the current host writer
    CHECK(values[0] == 42.0);

    catchem_diag_register(core, "transition", "writer transition", "1", 2, nc, 1, 0);
    {
        auto manager = static_cast<catchem::Core*>(core)->get_diagnostic_manager();
        auto* transition_host = static_cast<double*>(manager->get_host_pointer("transition"));
        transition_host[0] = 17.0;
        manager->mark_host_modified("transition");
        auto transition_device = manager->get_device_view_2d("transition");
        CHECK(transition_device(0, 0) == 17.0);
        CHECK(manager->get_field("transition")->latest_writer == catchem::LatestWriter::DeviceCurrent);
    }

    catchem_core_destroy(core);
#ifdef CATCHEM_ENABLE_KOKKOS
    Kokkos::finalize();
#endif
    return 0;
}
