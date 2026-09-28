#include "catchem_api.hpp"
#include "catchem_core.hpp"
#include "catchem_diagnostic_manager.hpp"
#include "catchem_kokkos_compat.hpp"
#include "catchem_process_registry.hpp"
#include "catchem_state_manager.hpp"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

extern "C" {
void catchem_register_carbchem_cpp();
}

int main(int argc, char* argv[]) {
    Kokkos::initialize(argc, argv);
    {
        std::cout << "\n==========================================" << std::endl;
        std::cout << "RUNNING TEST: CarbChem Process Unit Test" << std::endl;
        std::cout << "==========================================\n" << std::endl;

        catchem_register_carbchem_cpp();
        assert(catchem::ProcessRegistry::get_instance().has_process("carbchem"));

        int n_cols = 4;
        int n_levels = 5;
        int n_species = 22;

        auto core = std::make_shared<catchem::Core>(n_cols, n_levels, n_species);
        auto state = core->get_state_manager();
        auto runtime_config = std::make_shared<catchem::ConfigManager>();
        runtime_config->load_from_file("CATChem_new_config.yml");
        state->attach_config_manager(runtime_config);

        std::string species_path = "CATChem_species.yml";
        std::vector<std::string> candidates = {species_path, "tests/" + species_path, "../tests/" + species_path,
                                               "../../tests/" + species_path};
        for (const auto& candidate : candidates) {
            std::ifstream f(candidate);
            if (f.good()) {
                species_path = candidate;
                break;
            }
        }
        state->load_species_config(species_path);

        std::vector<double> temperature(n_cols * n_levels, 280.0);
        std::vector<double> airden_dry(n_cols * n_levels, 1.2);
        std::vector<double> delp(n_cols * n_levels, 1000.0);
        std::vector<double> pmid(n_cols * n_levels, 90000.0);
        std::vector<double> chem_conc(n_cols * n_levels * n_species, 0.0);

        state->bind_met_field_3d("T", temperature.data());
        state->bind_met_field_3d("AIRDEN_DRY", airden_dry.data());
        state->bind_met_field_3d("DELP", delp.data());
        state->bind_met_field_3d("PMID", pmid.data());
        state->bind_unified_chemistry(chem_conc.data());
        state->clock().timestep = 3600.0;

        // Populate every mechanism species through the semantic concentration
        // view.  Distinct values make species-slot transposition visible while
        // keeping this test independent of any particular mechanism names.
        auto concentration = state->chemistry().conc->mdspan();
        for (int species = 0; species < n_species; ++species)
            for (int level = 0; level < n_levels; ++level)
                for (int column = 0; column < n_cols; ++column)
                    concentration(column, level, species) = 1.0 + 0.1 * static_cast<double>(species) +
                                                            0.01 * static_cast<double>(level) +
                                                            0.001 * static_cast<double>(column);
        const std::vector<double> initial_concentration = chem_conc;

        auto carbchem = catchem::ProcessRegistry::get_instance().create("carbchem");
        assert(carbchem != nullptr);
        carbchem->init(state);

        // The configured diagnostic set controls the packed species extent.
        auto configured_diag_species = runtime_config->data.processes.at("carbchem").diag_species;
        if (configured_diag_species.empty()) {
            configured_diag_species.reserve(state->chemistry().species_list.size());
            for (const auto& species : state->chemistry().species_list)
                configured_diag_species.push_back(species.short_name);
        }
        const int n_diag = static_cast<int>(configured_diag_species.size());
        {
            const auto manager = core->get_diagnostic_manager();
            assert(manager->has_field("carbchem_prod_mass"));
            assert(manager->get_field("carbchem_prod_mass")->dimensions ==
                   std::vector<int>({n_cols, n_levels, n_diag}));
            assert(manager->get_field("carbchem_loss_flux")->dimensions == std::vector<int>({n_cols, n_diag}));
            assert(manager->get_field("carbchem_phobic_mass")->dimensions ==
                   std::vector<int>({n_cols, n_levels, n_diag}));
            assert(manager->get_field("carbchem_phobic_flux")->dimensions == std::vector<int>({n_cols, n_diag}));
            assert(manager->get_unpack_labels("carbchem_phobic_flux") == configured_diag_species);
            std::cout << "  PASS configured diagnostic set: fields and labels follow diag_species" << std::endl;
        }

        // A diag_species subset must shrink the species dimension to 1 and map
        // to the GLOBAL 1-based species index (the space the scheme matches on).
        // register_field throws on different-dims re-registration, so this runs
        // on a fresh Core.  init() reads only config + mechanism.
        {
            auto core2 = std::make_shared<catchem::Core>(n_cols, n_levels, n_species);
            auto state2 = core2->get_state_manager();
            auto cfg2 = std::make_shared<catchem::ConfigManager>();
            cfg2->load_from_file("CATChem_new_config.yml");
            cfg2->data.processes["carbchem"].diagnostics = true;
            state2->attach_config_manager(cfg2);
            state2->load_species_config(species_path);
            const std::string selected_species = state2->chemistry().species_list.front().short_name;
            cfg2->data.processes["carbchem"].diag_species = {selected_species};
            auto carbchem2 = catchem::ProcessRegistry::get_instance().create("carbchem");
            carbchem2->init(state2);
            const auto mgr2 = core2->get_diagnostic_manager();
            assert(mgr2->get_field("carbchem_prod_mass")->dimensions == std::vector<int>({n_cols, n_levels, 1}));
            assert(mgr2->get_field("carbchem_loss_flux")->dimensions == std::vector<int>({n_cols, 1}));
            std::cout << "  PASS diag_species subset: fields register with species dim = 1" << std::endl;
        }

        carbchem->run(state);
        state->sync_to_host();

        // Verify diagnostic values, not only registration.  For every packed
        // slot the column flux must be the pressure-weighted integral of the
        // reported per-level transfer mass.  For species handled by the active
        // scheme, that mass must also equal the absolute concentration change.
        {
            constexpr double gravity = 9.80665;
            const auto manager = core->get_diagnostic_manager();
            const auto labels = manager->get_unpack_labels("carbchem_phobic_mass");
            const auto* mass = static_cast<const double*>(manager->get_host_read_pointer("carbchem_phobic_mass"));
            const auto* flux = static_cast<const double*>(manager->get_host_read_pointer("carbchem_phobic_flux"));
            assert(mass != nullptr && flux != nullptr);
            assert(labels.size() == static_cast<std::size_t>(n_diag));

            bool observed_transfer = false;
            for (int slot = 0; slot < n_diag; ++slot) {
                std::string canonical = labels[static_cast<std::size_t>(slot)];
                std::transform(canonical.begin(), canonical.end(), canonical.begin(),
                               [](unsigned char value) { return std::toupper(value); });
                const auto species_it = state->chemistry().species_name_to_index.find(canonical);
                assert(species_it != state->chemistry().species_name_to_index.end());
                const int species = species_it->second;

                for (int column = 0; column < n_cols; ++column) {
                    double expected_flux = 0.0;
                    bool slot_transferred = false;
                    for (int level = 0; level < n_levels; ++level) {
                        const std::size_t mass_index =
                            static_cast<std::size_t>(column) +
                            static_cast<std::size_t>(n_cols) *
                                (static_cast<std::size_t>(level) + static_cast<std::size_t>(n_levels) * slot);
                        const double transferred_mass = mass[mass_index];
                        assert(std::isfinite(transferred_mass));
                        assert(transferred_mass >= 0.0);
                        expected_flux +=
                            transferred_mass *
                            delp[static_cast<std::size_t>(column) + static_cast<std::size_t>(n_cols) * level] /
                            gravity / state->clock().timestep;

                        const std::size_t concentration_index =
                            static_cast<std::size_t>(column) +
                            static_cast<std::size_t>(n_cols) *
                                (static_cast<std::size_t>(level) + static_cast<std::size_t>(n_levels) * species);
                        const double concentration_change =
                            std::abs(chem_conc[concentration_index] - initial_concentration[concentration_index]) *
                            1.0e-9;
                        if (transferred_mass > 1.0e-30) {
                            slot_transferred = true;
                            const double scale = std::max(transferred_mass, 1.0e-30);
                            assert(std::abs(concentration_change - transferred_mass) / scale < 2.0e-5);
                        }
                    }
                    const std::size_t flux_index =
                        static_cast<std::size_t>(column) + static_cast<std::size_t>(n_cols) * slot;
                    assert(std::isfinite(flux[flux_index]));
                    const double flux_scale = std::max(std::abs(expected_flux), 1.0e-30);
                    assert(std::abs(flux[flux_index] - expected_flux) / flux_scale < 2.0e-5);
                    observed_transfer = observed_transfer || slot_transferred;
                }
            }
            assert(observed_transfer && "configured diagnostics must include at least one species handled by CarbChem");
            std::cout << "  PASS packed diagnostic values follow concentration transfer and column flux" << std::endl;
        }

        std::cout << "SUCCESS: CarbChem process executed successfully." << std::endl;
    }
    Kokkos::finalize();
    return 0;
}
