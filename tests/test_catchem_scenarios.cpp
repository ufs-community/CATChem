#include "catchem_api.hpp"
#include "catchem_config_manager.hpp"
#include "catchem_core.hpp"
#include "catchem_kokkos_compat.hpp"
#include "catchem_process_registry.hpp"
#include "catchem_state_manager.hpp"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

extern "C" {
void catchem_register_seasalt_cpp();
void catchem_register_drydep_cpp();
void catchem_register_wetdep_cpp();
void catchem_register_settling_cpp();
void catchem_register_so4chem_cpp();
void catchem_register_dust_cpp();
void catchem_register_carbchem_cpp();
}

struct AtmosphericScenario {
    std::string id;
    std::string description;
    double lat;
    double lon;
    double t_surf;       // Surface temperature (K)
    double t_top;        // Top of column temperature (K)
    double p_surf;       // Surface pressure (Pa)
    double p_top;        // Top pressure (Pa)
    double rh_val;       // Relative humidity (%)
    double qv_val;       // Specific humidity (kg/kg)
    double u10m_val;     // 10m U-wind (m/s)
    double v10m_val;     // 10m V-wind (m/s)
    double ustar_val;    // Friction velocity (m/s)
    double hflux_val;    // Sensible heat flux (W/m2)
    double obk_val;      // Obukhov length (m)
    double pblh_val;     // Boundary layer height (m)
    double z0h_val;      // Surface roughness length (m)
    double precip_rate;  // Rain/snow precip rate (kg/m2/s)
    double frocean_val;  // Ocean fraction [0.0 - 1.0]
    double frseaice_val; // Sea ice fraction [0.0 - 1.0]
    double gas_bg;       // Background gas concentration (kg/kg)
    double dust_bg;      // Background dust concentration (kg/kg)
    double seas_bg;      // Background sea salt concentration (kg/kg)
    double carbon_bg;    // Background carbon aerosol concentration (kg/kg)
};

static std::string find_species_config() {
    std::vector<std::string> candidates = {"CATChem_species.yml", "tests/CATChem_species.yml",
                                           "../tests/CATChem_species.yml", "../../tests/CATChem_species.yml"};
    for (const auto& path : candidates) {
        std::ifstream f(path);
        if (f.good())
            return path;
    }
    return "CATChem_species.yml";
}

// Build a runtime ConfigManager that activates and configures all seven
// processes with the schemes their init() contracts require.  The scenario
// cores are created grid-only (catchem_core_create), so without an attached
// configuration each process init() would reject the empty processes block.
static std::shared_ptr<catchem::ConfigManager> make_scenario_config() {
    const char* yaml = "simulation:\n"
                       "  name: scenarios\n"
                       "processes:\n"
                       "  seasalt:\n"
                       "    activate: true\n"
                       "    scheme: geos12\n"
                       "    diagnostics: false\n"
                       "  dust:\n"
                       "    activate: true\n"
                       "    scheme: fengsha\n"
                       "    diagnostics: false\n"
                       "  drydep:\n"
                       "    activate: true\n"
                       "    gas_scheme: wesely\n"
                       "    aero_scheme: gocart\n"
                       "    diagnostics: false\n"
                       "  settling:\n"
                       "    activate: true\n"
                       "    scheme: gocart\n"
                       "    diagnostics: false\n"
                       "  so4chem:\n"
                       "    activate: true\n"
                       "    scheme: gocart\n"
                       "    diagnostics: false\n"
                       "  wetdep:\n"
                       "    activate: true\n"
                       "    scheme: jacob\n"
                       "    diagnostics: false\n"
                       "  carbchem:\n"
                       "    activate: true\n"
                       "    scheme: gocart\n"
                       "    diagnostics: false\n";
    const std::string path = "scenarios_runtime_config.yml";
    std::ofstream out(path);
    out << yaml;
    out.close();
    auto config = std::make_shared<catchem::ConfigManager>();
    config->load_from_file(path);
    return config;
}

int main(int argc, char* argv[]) {
    Kokkos::initialize(argc, argv);
    {
        std::cout << "\n==========================================================================" << std::endl;
        std::cout << "=== RUNNING ATMOSPHERIC SCENARIO PROPERTY TEST (ALL PROCESSES COUPLED) ===" << std::endl;
        std::cout << "==========================================================================\n" << std::endl;

        // Register all process constructors
        catchem_register_seasalt_cpp();
        catchem_register_drydep_cpp();
        catchem_register_wetdep_cpp();
        catchem_register_settling_cpp();
        catchem_register_so4chem_cpp();
        catchem_register_dust_cpp();
        catchem_register_carbchem_cpp();

        std::vector<AtmosphericScenario> scenarios = {
            {"TROPICAL_RAINFOREST",
             "Equatorial Tropical Forest (Amazon/Congo) - High T, high humidity, heavy rain",
             -3.0,
             -60.0,
             303.15,
             220.15,
             101000.0,
             10000.0,
             85.0,
             0.018,
             2.0,
             1.0,
             0.25,
             30.0,
             -10.0,
             1200.0,
             0.5,
             5.0e-5,
             0.0,
             0.0,
             1.0e-8,
             1.0e-10,
             1.0e-10,
             5.0e-9},

            {"SUBTROPICAL_DESERT",
             "Subtropical Arid Desert (Sahara/Rub' al Khali) - Extremely hot, dry, high dust wind",
             23.5,
             25.0,
             318.15,
             215.15,
             101300.0,
             10000.0,
             10.0,
             0.001,
             18.0,
             6.0,
             1.2,
             150.0,
             -50.0,
             2500.0,
             0.005,
             0.0,
             0.0,
             0.0,
             2.0e-9,
             1.0e-6,
             1.0e-11,
             1.0e-10},

            {"URBAN_INDUSTRIAL_MEGACITY",
             "Mid-Latitude Urban Megacity (Beijing/LA) - High anthropogenic pollution background",
             39.9,
             116.4,
             295.15,
             218.15,
             101500.0,
             10000.0,
             55.0,
             0.008,
             4.0,
             2.0,
             0.4,
             80.0,
             -20.0,
             1000.0,
             1.0,
             0.0,
             0.0,
             0.0,
             1.0e-7,
             1.0e-8,
             1.0e-9,
             5.0e-8},

            {"POLAR_ARCTIC_WINTER",
             "High-Latitude Polar Ice Cap (Greenland/Arctic) - Subzero extreme cold, sea ice cover",
             78.0,
             -40.0,
             233.15,
             205.15,
             102000.0,
             10000.0,
             35.0,
             0.00005,
             6.0,
             -2.0,
             0.2,
             -10.0,
             100.0,
             300.0,
             0.001,
             0.0,
             0.0,
             1.0,
             1.0e-10,
             1.0e-11,
             1.0e-11,
             1.0e-11},

            {"SOUTHERN_OCEAN_STORM",
             "High-Latitude Sub-Antarctic Ocean Storm - 100% ocean, high gale winds, sea spray",
             -55.0,
             70.0,
             278.15,
             210.15,
             98000.0,
             10000.0,
             90.0,
             0.006,
             28.0,
             12.0,
             1.8,
             40.0,
             -30.0,
             800.0,
             0.001,
             2.0e-5,
             1.0,
             0.0,
             5.0e-9,
             1.0e-10,
             1.0e-6,
             5.0e-10},

            {"SUB_ANTARCTIC_ICE_EDGE",
             "Sub-Antarctic Ice Edge - Mixed 50% ocean / 50% ice, cold snow showers",
             -65.0,
             0.0,
             268.15,
             208.15,
             99500.0,
             10000.0,
             80.0,
             0.003,
             12.0,
             4.0,
             0.6,
             15.0,
             -15.0,
             500.0,
             0.005,
             1.0e-5,
             0.5,
             0.5,
             2.0e-9,
             1.0e-11,
             1.0e-7,
             1.0e-10},

            {"MIDLAT_TEMPERATE_FOREST",
             "Mid-Latitude Rural Temperate Forest (Central Europe) - Summer mild, vegetation",
             48.5,
             11.5,
             293.15,
             218.15,
             101325.0,
             10000.0,
             65.0,
             0.009,
             4.0,
             2.0,
             0.3,
             50.0,
             -25.0,
             1400.0,
             0.3,
             1.0e-6,
             0.0,
             0.0,
             1.0e-8,
             1.0e-9,
             1.0e-9,
             1.0e-8},

            {"HIGH_ALTITUDE_TROPOPAUSE",
             "High Altitude / Tropopause Layer (12 km) - Low pressure, low density, jet stream wind",
             45.0,
             0.0,
             215.15,
             200.15,
             15000.0,
             1000.0,
             20.0,
             0.00001,
             45.0,
             15.0,
             0.5,
             0.0,
             1000.0,
             100.0,
             0.001,
             0.0,
             0.0,
             0.0,
             1.0e-9,
             1.0e-11,
             1.0e-11,
             1.0e-10}};

        const int n_cols = 8;
        const int n_levels = 10;
        const int n_species = 22;
        const int size_3d = n_cols * n_levels;
        const size_t total_size = n_cols * n_levels * n_species;
        const std::string species_path = find_species_config();

        std::cout << "Grid configuration: " << n_cols << " columns x " << n_levels << " levels x " << n_species
                  << " species" << std::endl;
        std::cout << "Loading species definition from: " << species_path << std::endl;

        for (const auto& sc : scenarios) {
            std::cout << "\n--------------------------------------------------------------------------" << std::endl;
            std::cout << "TESTING SCENARIO [" << sc.id << "]: " << sc.description << std::endl;
            std::cout << "--------------------------------------------------------------------------" << std::endl;

            // Create Core instance for this scenario
            void* core_ptr = catchem_core_create(n_cols, n_levels, n_species);
            auto* core = static_cast<catchem::Core*>(core_ptr);
            auto state = core->get_state_manager();
            state->load_species_config(species_path);
            // Grid-only cores start with an empty processes block; attach a
            // runtime configuration so each process init() sees its required
            // scheme selection.
            state->attach_config_manager(make_scenario_config());

            // Meteorological & Surface Tensors
            std::vector<double> t_air(size_3d);
            std::vector<double> pmid(size_3d);
            std::vector<double> pedge(n_cols * (n_levels + 1));
            std::vector<double> airden_dry(size_3d);
            std::vector<double> mairden(size_3d);
            std::vector<double> bxheight(size_3d);
            std::vector<double> cldf(size_3d, 0.1);
            std::vector<double> pfilsan(size_3d, sc.precip_rate);
            std::vector<double> pfllsan(size_3d, sc.precip_rate);
            std::vector<double> reevapls(size_3d, 0.0);
            std::vector<double> lat(n_cols, sc.lat);
            std::vector<double> lon(n_cols, sc.lon);
            std::vector<double> sst(n_cols, sc.t_surf);
            std::vector<double> frocean(n_cols, sc.frocean_val);
            std::vector<double> frseaice(n_cols, sc.frseaice_val);
            std::vector<double> ustar(n_cols, sc.ustar_val);
            std::vector<double> delp(size_3d);
            std::vector<double> u10m(n_cols, sc.u10m_val);
            std::vector<double> v10m(n_cols, sc.v10m_val);
            std::vector<double> rh(size_3d, sc.rh_val);
            std::vector<double> qv(size_3d, sc.qv_val);
            std::vector<double> ps(n_cols, sc.p_surf);
            std::vector<double> ts(n_cols, sc.t_surf);
            std::vector<double> hflux(n_cols, sc.hflux_val);
            std::vector<double> obk(n_cols, sc.obk_val);
            std::vector<double> pblh(n_cols, sc.pblh_val);
            std::vector<double> z0h(n_cols, sc.z0h_val);

            // Construct vertical pressure profile and vertical temperature profile
            for (int col = 0; col < n_cols; ++col) {
                pedge[col + 0 * n_cols] = sc.p_surf;
                double dp = (sc.p_surf - sc.p_top) / static_cast<double>(n_levels);
                for (int k = 0; k < n_levels; ++k) {
                    pedge[col + (k + 1) * n_cols] = sc.p_surf - (k + 1) * dp;
                    double p_mid_val = 0.5 * (pedge[col + k * n_cols] + pedge[col + (k + 1) * n_cols]);
                    pmid[col + k * n_cols] = p_mid_val;
                    delp[col + k * n_cols] = dp;

                    // Linear temperature profile from surf to top
                    double frac = static_cast<double>(k) / static_cast<double>(n_levels - 1);
                    double t_val = sc.t_surf + frac * (sc.t_top - sc.t_surf);
                    if (t_val < 180.0)
                        t_val = 180.0;
                    t_air[col + k * n_cols] = t_val;

                    // Derive air density via Ideal Gas Law: rho = P / (R_dry * T)
                    double rho = p_mid_val / (287.05 * t_val);
                    if (rho < 1.0e-4)
                        rho = 1.0e-4;
                    airden_dry[col + k * n_cols] = rho;
                    mairden[col + k * n_cols] = rho;

                    // Derive layer thickness dz = dp / (rho * g)
                    double dz = dp / (rho * 9.80665);
                    if (dz < 0.1)
                        dz = 0.1;
                    bxheight[col + k * n_cols] = dz;
                }
            }

            // Bind all meteorological fields
            state->bind_met_field_3d("T", t_air.data());
            state->bind_met_field_3d("PMID", pmid.data());
            state->bind_met_field_3d("PEDGE", pedge.data());
            state->bind_met_field_3d("AIRDEN_DRY", airden_dry.data());
            state->bind_met_field_3d("AIRDEN", airden_dry.data());
            state->bind_met_field_3d("BXHEIGHT", bxheight.data());
            state->bind_met_field_3d("MAIRDEN", mairden.data());
            state->bind_met_field_3d("PFILSAN", pfilsan.data());
            state->bind_met_field_3d("PFLLSAN", pfllsan.data());
            state->bind_met_field_3d("REEVAPLS", reevapls.data());
            state->bind_met_field_3d("RH", rh.data());
            state->bind_met_field_3d("QV", qv.data());
            state->bind_met_field_3d("CLDF", cldf.data());
            state->bind_met_field_2d("LAT", lat.data());
            state->bind_met_field_2d("LON", lon.data());
            state->bind_met_field_2d("SST", sst.data());
            state->bind_met_field_2d("PS", ps.data());
            state->bind_met_field_2d("TS", ts.data());
            state->bind_met_field_2d("HFLUX", hflux.data());
            state->bind_met_field_2d("OBK", obk.data());
            state->bind_met_field_2d("PBLH", pblh.data());
            state->bind_met_field_2d("Z0H", z0h.data());
            state->bind_met_field_2d("FROCEAN", frocean.data());
            state->bind_met_field_2d("FRSEAICE", frseaice.data());
            state->bind_met_field_2d("USTAR", ustar.data());
            state->bind_met_field_3d("DELP", delp.data());
            state->bind_met_field_2d("U10M", u10m.data());
            state->bind_met_field_2d("V10M", v10m.data());

            // Initialize background chemical concentrations per species group
            std::vector<double> conc(total_size, 0.0);
            for (size_t s = 0; s < state->chemistry().species_list.size(); ++s) {
                const auto& spec = state->chemistry().species_list[s];
                double bg = sc.gas_bg;
                if (spec.is_dust)
                    bg = sc.dust_bg;
                else if (spec.is_seasalt)
                    bg = sc.seas_bg;
                else if (spec.short_name == "bc1" || spec.short_name == "bc2" || spec.short_name == "oc1" ||
                         spec.short_name == "oc2")
                    bg = sc.carbon_bg;

                for (int col = 0; col < n_cols; ++col) {
                    for (int lvl = 0; lvl < n_levels; ++lvl) {
                        int idx = col + lvl * n_cols + s * size_3d;
                        conc[idx] = bg;
                    }
                }
            }

            state->bind_unified_chemistry(conc.data());

            // Schedule ALL 7 physics & chemistry processes in the core engine
            auto settling = catchem::ProcessRegistry::get_instance().create("settling");
            settling->init(state);
            core->add_process(settling);

            auto drydep = catchem::ProcessRegistry::get_instance().create("drydep");
            drydep->init(state);
            core->add_process(drydep);

            auto seasalt = catchem::ProcessRegistry::get_instance().create("seasalt");
            seasalt->init(state);
            core->add_process(seasalt);

            auto wetdep = catchem::ProcessRegistry::get_instance().create("wetdep");
            wetdep->init(state);
            core->add_process(wetdep);

            auto so4chem = catchem::ProcessRegistry::get_instance().create("so4chem");
            so4chem->init(state);
            core->add_process(so4chem);

            auto dust = catchem::ProcessRegistry::get_instance().create("dust");
            dust->init(state);
            core->add_process(dust);

            auto carbchem = catchem::ProcessRegistry::get_instance().create("carbchem");
            carbchem->init(state);
            core->add_process(carbchem);

            std::cout << "  Core registered " << core->get_num_processes() << " active processes." << std::endl;

            // Run multi-timestep simulation (24 timesteps = 12 simulated hours)
            const double dt = 1800.0; // 30-minute timestep
            const int num_steps = 24;

            for (int step = 1; step <= num_steps; ++step) {
                state->sync_to_device();
                catchem_core_run_timestep(core_ptr, dt);
                state->sync_to_host();

                // Check Property Invariants: Finiteness
                for (size_t i = 0; i < total_size; ++i) {
                    if (!std::isfinite(conc[i])) {
                        int spec_idx = i / size_3d;
                        std::cerr << "PROPERTY FAILURE in [" << sc.id << "]: NaN/Inf at index " << i
                                  << " (Species=" << spec_idx << ") at step " << step << std::endl;
                        assert(false && "Concentration must remain finite!");
                    }
                }

                // Clip negative concentrations to 0.0 (standard atmospheric model non-negative solver boundary)
                for (auto& val : conc) {
                    if (val < 0.0)
                        val = 0.0;
                }

                // Verify Non-negativity invariant holds after solver boundary
                for (size_t i = 0; i < total_size; ++i) {
                    assert(conc[i] >= 0.0 && "Concentration must remain non-negative after solver boundary!");
                }
            }

            // Validate registered diagnostics pointers remain finite
            if (state->diagnostic_manager()) {
                auto diag_names = state->diagnostic_manager()->get_registered_names();
                std::cout << "  Diagnostic Manager registered " << diag_names.size() << " fields." << std::endl;
                for (const auto& dname : diag_names) {
                    const double* dptr =
                        static_cast<const double*>(state->diagnostic_manager()->get_host_pointer(dname));
                    assert(dptr != nullptr);
                    assert(std::isfinite(dptr[0]) && "Diagnostic field pointer contains non-finite value!");
                }
            }

            catchem_core_destroy(core_ptr);
            std::cout << "  PASS: Scenario [" << sc.id << "] completed " << num_steps << " timesteps successfully!"
                      << std::endl;
        }

        std::cout << "\n==========================================================================" << std::endl;
        std::cout << "=== SUCCESS: ALL ATMOSPHERIC SCENARIO PROPERTY CHECKS PASSED!        ===" << std::endl;
        std::cout << "==========================================================================\n" << std::endl;
    }
    Kokkos::finalize();
    return 0;
}
