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
void catchem_register_dust_cpp();
}

int main(int argc, char* argv[]) {
    Kokkos::initialize(argc, argv);
    {
        std::cout << "\n==========================================" << std::endl;
        std::cout << "RUNNING TEST: Dust Process Unit Test" << std::endl;
        std::cout << "==========================================\n" << std::endl;

        catchem_register_dust_cpp();
        assert(catchem::ProcessRegistry::get_instance().has_process("dust"));

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

        // The science bridge consumes dust bins in species-list order.  Keep
        // this mechanism-independent: the active dust metadata must define a
        // physically ordered bin sequence, not a set of tracer-name rules.
        std::vector<size_t> dust_indices;
        for (size_t index = 0; index < state->chemistry().species_list.size(); ++index) {
            if (state->chemistry().species_list[index].is_dust)
                dust_indices.push_back(index);
        }
        assert(!dust_indices.empty());
        for (size_t bin = 0; bin < dust_indices.size(); ++bin) {
            const auto& species = state->chemistry().species_list[dust_indices[bin]];
            assert(species.lower_radius < species.radius && species.radius < species.upper_radius);
            if (bin > 0)
                assert(state->chemistry().species_list[dust_indices[bin - 1]].radius < species.radius);
        }

        std::vector<double> u10m(n_cols, 15.0);
        std::vector<double> v10m(n_cols, 0.0);
        std::vector<double> ustar(n_cols, 0.8);
        std::vector<double> ustar_threshold(n_cols, 0.4);
        std::vector<double> airden(n_cols * n_levels, 1.2);
        std::vector<double> delp(n_cols * n_levels, 1000.0);
        std::vector<double> clay_fraction(n_cols, 0.2);
        std::vector<double> lake_fraction(n_cols, 0.0);
        std::vector<double> snow_fraction(n_cols, 0.0);
        std::vector<double> vegetation_fraction(n_cols, 0.3);
        std::vector<double> leaf_area_index(n_cols, 1.0);
        std::vector<double> bxheight(n_cols * n_levels, 100.0);
        std::vector<double> rdrag(n_cols, 0.01);
        std::vector<double> sand_fraction(n_cols, 0.4);
        std::vector<double> surface_soil_moisture(n_cols, 0.1);
        std::vector<double> skin_temperature(n_cols, 290.0);
        std::vector<double> roughness_length(n_cols, 0.05);
        std::vector<double> soil_moisture(n_cols * n_levels, 0.1);
        std::vector<double> lwi(n_cols, 1.0), ssm(n_cols, 0.1), z0(n_cols, 0.05);
        std::vector<double> chem_conc(n_cols * n_levels * n_species, 0.0);

        state->bind_met_field_3d("air_density", airden.data());
        state->bind_met_field_3d("DELP", delp.data());
        state->bind_met_field_3d("box_height", bxheight.data());
        state->bind_met_field_2d("clay_fraction", clay_fraction.data());
        state->bind_met_field_2d("lake_fraction", lake_fraction.data());
        state->bind_met_field_2d("snow_fraction", snow_fraction.data());
        state->bind_met_field_2d("vegetation_fraction", vegetation_fraction.data());
        const auto gvf_field = state->find_field<2>("GVF");
        assert(gvf_field != nullptr);
        assert(gvf_field->contract.units == "frac");
        state->bind_met_field_2d("leaf_area_index", leaf_area_index.data());
        state->bind_met_field_2d("RDRAG", rdrag.data());
        state->bind_met_field_2d("sand_fraction", sand_fraction.data());
        state->bind_met_field_3d("soil_moisture", soil_moisture.data());
        state->bind_met_field_2d("surface_soil_moisture", surface_soil_moisture.data());
        state->bind_met_field_2d("LWI", lwi.data());
        state->bind_met_field_2d("SSM", ssm.data());
        state->bind_met_field_2d("Z0", z0.data());
        state->bind_met_field_2d("skin_temperature", skin_temperature.data());
        state->bind_met_field_2d("u_10m", u10m.data());
        state->bind_met_field_2d("v_10m", v10m.data());
        state->bind_met_field_2d("friction_velocity", ustar.data());
        state->bind_met_field_2d("threshold_friction_velocity", ustar_threshold.data());
        state->bind_met_field_2d("roughness_length", roughness_length.data());
        state->bind_unified_chemistry(chem_conc.data());

        auto dust = catchem::ProcessRegistry::get_instance().create("dust");
        assert(dust != nullptr);
        // The shared test config leaves processes.dust.diagnostics off; force
        // it on so init() registers the diagnostic fields this test checks.
        runtime_config->data.processes["dust"].diagnostics = true;
        dust->init(state);

        // Per-bin diagnostics must register as compact [ncols, n_dust] fields
        // (not full-width species arrays) so the NUOPC driver can emit one 3D
        // (nx, ny, nbin) variable; per-column diagnostics stay [ncols, 1].
        {
            const auto manager = core->get_diagnostic_manager();
            const int n_dust = static_cast<int>(dust_indices.size());
            assert(manager->has_field("dust_emission_bin"));
            assert(manager->get_field("dust_emission_bin")->dimensions == std::vector<int>({n_cols, n_dust}));
            assert(manager->has_field("dust_utar_threshold"));
            assert(manager->get_field("dust_utar_threshold")->dimensions == std::vector<int>({n_cols, n_dust}));
            assert(manager->get_field("dust_emission_total")->dimensions == std::vector<int>({n_cols, 1}));
        }

        // A diag_species subset must shrink the per-bin fields to [ncols, 1].
        // register_field throws "Incompatible diagnostic re-registration" on
        // different-dims re-registration, so this MUST run on a fresh
        // Core/state.  init() reads only the config and the chemistry
        // mechanism, so no met fields need binding here.
        {
            auto core2 = std::make_shared<catchem::Core>(n_cols, n_levels, n_species);
            auto state2 = core2->get_state_manager();
            auto cfg2 = std::make_shared<catchem::ConfigManager>();
            cfg2->load_from_file("CATChem_new_config.yml");
            cfg2->data.processes["dust"].diagnostics = true;
            cfg2->data.processes["dust"].diag_species = {
                state->chemistry().species_list[dust_indices.back()].short_name};
            state2->attach_config_manager(cfg2);
            state2->load_species_config(species_path);
            auto dust2 = catchem::ProcessRegistry::get_instance().create("dust");
            dust2->init(state2);
            const auto mgr2 = core2->get_diagnostic_manager();
            assert(mgr2->get_field("dust_emission_bin")->dimensions == std::vector<int>({n_cols, 1}));
            assert(mgr2->get_field("dust_utar_threshold")->dimensions == std::vector<int>({n_cols, 1}));
            assert(mgr2->get_field("dust_emission_total")->dimensions == std::vector<int>({n_cols, 1}));
            std::cout << "  PASS diag_species subset: per-bin fields register as [ncols, 1]" << std::endl;
        }

        // FENGSHA gates emission on saltation: the White horizontal flux is
        // max(0, R*ustar - ustar_threshold*H) * (...)^2, so emission only
        // occurs when the drag-scaled friction velocity exceeds the moisture-
        // adjusted threshold.  Exercise BOTH regimes through the real process
        // to prove the gate works in each direction, not just that output is
        // finite/non-negative.  The met arrays are bound by pointer, so we
        // mutate the gating inputs in place and re-run.
        const auto surface_index = [&](size_t species_index) {
            // Surface layer (level 0) of the given species, column 0.
            return species_index * static_cast<size_t>(n_cols * n_levels);
        };
        const auto total_dust = [&]() {
            double sum = 0.0;
            for (const auto species_index : dust_indices)
                for (int col = 0; col < n_cols; ++col)
                    sum += chem_conc[species_index * static_cast<size_t>(n_cols * n_levels) + col];
            return sum;
        };
        const auto reset_chem = [&]() { std::fill(chem_conc.begin(), chem_conc.end(), 0.0); };

        // --- Scenario A: saltation-favorable -> dust MUST emit --------------
        // Strong wind, low threshold, near-unity drag partition, dry soil.
        std::fill(ustar.begin(), ustar.end(), 0.8);
        std::fill(ustar_threshold.begin(), ustar_threshold.end(), 0.15);
        std::fill(rdrag.begin(), rdrag.end(), 1.0);
        std::fill(surface_soil_moisture.begin(), surface_soil_moisture.end(), 0.02);
        std::fill(soil_moisture.begin(), soil_moisture.end(), 0.02);
        reset_chem();
        dust->run(state);
        state->sync_to_host();

        for (const auto species_index : dust_indices)
            assert(chem_conc[surface_index(species_index)] >= 0.0); // physical floor
        const double emitting_total = total_dust();
        assert(emitting_total > 0.0 && "FENGSHA must emit dust under saltation-favorable inputs");
        std::cout << "  Scenario A (emitting): total surface dust = " << emitting_total << std::endl;

        // Regression guard for the diagnostic index space: the scheme fills
        // per-bin slots by LOCAL bin position (diagnostic_species_id(diag_idx)
        // == species_idx over the 1..n_dust slice).  With the old GLOBAL ids
        // the compact field stayed ~all zeros; under the default (all-bin)
        // diag set the emitting column must show positive per-bin emission.
        {
            const auto manager = core->get_diagnostic_manager();
            const int n_dust = static_cast<int>(dust_indices.size());
            const double* emission_bin =
                static_cast<const double*>(manager->get_host_pointer("dust_emission_bin"));
            assert(emission_bin != nullptr);
            bool any_positive = false;
            for (int i = 0; i < n_cols * n_dust; ++i) {
                assert(std::isfinite(emission_bin[i]));
                assert(emission_bin[i] >= 0.0);
                if (emission_bin[i] > 0.0)
                    any_positive = true;
            }
            assert(any_positive && "dust_emission_bin must be populated in LOCAL bin index space");
            std::cout << "  Scenario A (emitting): dust_emission_bin populated" << std::endl;
        }

        // --- Scenario B: saltation-suppressed -> dust MUST be zero ----------
        // Same wind, but the drag-scaled friction velocity R*ustar is far
        // below the threshold, so no saltation and therefore no emission.
        std::fill(ustar.begin(), ustar.end(), 0.8);
        std::fill(ustar_threshold.begin(), ustar_threshold.end(), 2.0);
        std::fill(rdrag.begin(), rdrag.end(), 0.01);
        reset_chem();
        dust->run(state);
        state->sync_to_host();

        const double suppressed_total = total_dust();
        assert(suppressed_total == 0.0 && "FENGSHA must not emit dust when saltation is suppressed");
        std::cout << "  Scenario B (suppressed): total surface dust = " << suppressed_total << std::endl;

        std::cout << "SUCCESS: Dust process emits under favorable inputs and is silent when suppressed." << std::endl;
    }
    Kokkos::finalize();
    return 0;
}
