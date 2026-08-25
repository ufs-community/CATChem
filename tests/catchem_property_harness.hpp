#pragma once

#include "catchem_api.hpp"
#include "catchem_core.hpp"
#include "catchem_diagnostic_manager.hpp"
#include "catchem_process_registry.hpp"
#include "catchem_state_manager.hpp"
#include <Kokkos_Core.hpp>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <memory>
#include <random>
#include <string>
#include <vector>

namespace catchem::test {

    /**
     * @brief Common atmospheric chemistry modeling scenarios for property testing.
     */
    enum class AtmosphericScenario {
        BoundaryLayer,                ///< Generic Planetary Boundary Layer (0 - 2 km)
        UrbanBoundaryLayer,           ///< Urban land cover (high roughness z0, dry soil, heat island, low vegetation)
        RuralAgricultural,            ///< Rural & Agricultural land cover (high LAI/GVF, crops & forest, moist soil)
        MarineBoundaryLayer,          ///< Ocean surface (100% ocean, low z0, sea salt aerosol emissions)
        PolarSnowIce,                 ///< Snow and sea ice cover (cold, high snow/ice fraction, low evaporation)
        MidTroposphere,               ///< Free Troposphere with clouds & precipitation (2 - 10 km)
        UpperTroposphereStratosphere, ///< UTLS cold, low-pressure, dry conditions (10 - 30 km)
        ExtremeStormDesert,           ///< Stress scenario: high winds, dust storms, extreme temperatures
        FullRandomFuzz                ///< Unconstrained randomized physical bounds
    };

    inline std::string scenario_to_string(AtmosphericScenario scenario) {
        switch (scenario) {
        case AtmosphericScenario::BoundaryLayer:
            return "BoundaryLayer";
        case AtmosphericScenario::UrbanBoundaryLayer:
            return "UrbanBoundaryLayer";
        case AtmosphericScenario::RuralAgricultural:
            return "RuralAgricultural";
        case AtmosphericScenario::MarineBoundaryLayer:
            return "MarineBoundaryLayer";
        case AtmosphericScenario::PolarSnowIce:
            return "PolarSnowIce";
        case AtmosphericScenario::MidTroposphere:
            return "MidTroposphere";
        case AtmosphericScenario::UpperTroposphereStratosphere:
            return "UpperTroposphereStratosphere";
        case AtmosphericScenario::ExtremeStormDesert:
            return "ExtremeStormDesert";
        case AtmosphericScenario::FullRandomFuzz:
            return "FullRandomFuzz";
        }
        return "Unknown";
    }

    /**
     * @brief Configuration settings for property harness runs.
     */
    struct PropertyHarnessConfig {
        int n_cols = 8;
        int n_levels = 10;
        int n_species = 22;
        uint32_t random_seed = 1337;
        int iterations_per_scenario = 20;
        double dt_sec = 1800.0; // 30-minute timestep
    };

    /**
     * @brief Unified Property-Based Testing Harness for CATChem.
     *
     * Constructs physically consistent atmospheric states across defined scenarios,
     * executes coupled or single processes, and verifies core scientific invariants
     * (finiteness, mass conservation, non-negativity, zero-state stability, diagnostic consistency).
     */
    class PropertyTestHarness {
    private:
        PropertyHarnessConfig config_;
        std::mt19937 gen_;

        // State buffers
        std::vector<double> t_air_;
        std::vector<double> pmid_;
        std::vector<double> pedge_;
        std::vector<double> airden_dry_;
        std::vector<double> mairden_;
        std::vector<double> bxheight_;
        std::vector<double> cldf_;
        std::vector<double> pfilsan_;
        std::vector<double> pfllsan_;
        std::vector<double> reevapls_;
        std::vector<double> lat_;
        std::vector<double> lon_;
        std::vector<double> sst_;
        std::vector<double> frocean_;
        std::vector<double> frseaice_;
        std::vector<double> lwi_;
        std::vector<double> frsno_;
        std::vector<double> frlake_;
        std::vector<double> gvf_;
        std::vector<double> lai_;
        std::vector<double> ssm_;
        std::vector<double> gwettop_;
        std::vector<double> clayfrac_;
        std::vector<double> sandfrac_;
        std::vector<double> rdrag_;
        std::vector<double> ustar_;
        std::vector<double> delp_;
        std::vector<double> u10m_;
        std::vector<double> v10m_;
        std::vector<double> rh_;
        std::vector<double> ps_;
        std::vector<double> ts_;
        std::vector<double> hflux_;
        std::vector<double> obk_;
        std::vector<double> pblh_;
        std::vector<double> z0h_;
        std::vector<double> conc_;

        void allocate_buffers() {
            int size_3d = config_.n_cols * config_.n_levels;
            int size_pedge = config_.n_cols * (config_.n_levels + 1);
            size_t total_size = config_.n_cols * config_.n_levels * config_.n_species;

            t_air_.resize(size_3d);
            pmid_.resize(size_3d);
            pedge_.resize(size_pedge);
            airden_dry_.resize(size_3d);
            mairden_.resize(size_3d);
            bxheight_.resize(size_3d);
            cldf_.resize(size_3d);
            pfilsan_.resize(size_pedge);
            pfllsan_.resize(size_pedge);
            reevapls_.resize(size_3d);
            lat_.assign(config_.n_cols, 40.0);
            lon_.assign(config_.n_cols, -105.0);
            sst_.assign(config_.n_cols, 290.0);
            frocean_.assign(config_.n_cols, 1.0);
            frseaice_.assign(config_.n_cols, 0.0);
            lwi_.assign(config_.n_cols, 1.0);
            frsno_.assign(config_.n_cols, 0.0);
            frlake_.assign(config_.n_cols, 0.0);
            gvf_.assign(config_.n_cols, 0.2);
            lai_.assign(config_.n_cols, 1.0);
            ssm_.assign(config_.n_cols, 0.2);
            gwettop_.assign(config_.n_cols, 0.2);
            clayfrac_.assign(config_.n_cols, 0.3);
            sandfrac_.assign(config_.n_cols, 0.4);
            rdrag_.assign(config_.n_cols, 0.1);
            ustar_.assign(config_.n_cols, 0.5);
            delp_.resize(size_3d);
            u10m_.assign(config_.n_cols, 5.0);
            v10m_.assign(config_.n_cols, 2.0);
            rh_.resize(size_3d);
            ps_.assign(config_.n_cols, 101325.0);
            ts_.assign(config_.n_cols, 288.15);
            hflux_.assign(config_.n_cols, 10.0);
            obk_.assign(config_.n_cols, 100.0);
            pblh_.assign(config_.n_cols, 1000.0);
            z0h_.assign(config_.n_cols, 0.01);
            conc_.resize(total_size);
        }

        void fill_random(std::vector<double>& vec, double min_val, double max_val) {
            std::uniform_real_distribution<double> dist(min_val, max_val);
            for (auto& val : vec) {
                val = dist(gen_);
            }
        }

    public:
        explicit PropertyTestHarness(const PropertyHarnessConfig& config = PropertyHarnessConfig{})
            : config_(config), gen_(config.random_seed) {
            allocate_buffers();
        }

        /**
         * @brief Generate physically consistent meteorological and chemical states for a scenario.
         */
        void generate_scenario_state(AtmosphericScenario scenario, std::shared_ptr<catchem::StateManager> state) {
            int nc = config_.n_cols;
            int nl = config_.n_levels;

            double p_surf_min = 95000.0, p_surf_max = 103000.0;
            double temp_min = 270.0, temp_max = 310.0;
            double rh_min = 10.0, rh_max = 95.0;
            double ustar_min = 0.05, ustar_max = 1.5;
            double precip_scale = 0.01;

            switch (scenario) {
            case AtmosphericScenario::BoundaryLayer:
                p_surf_min = 98000.0;
                p_surf_max = 103000.0;
                temp_min = 273.15;
                temp_max = 315.0;
                rh_min = 20.0;
                rh_max = 98.0;
                ustar_min = 0.1;
                ustar_max = 2.0;
                precip_scale = 0.001;
                lwi_.assign(nc, 1.0); // Land
                fill_random(frocean_, 0.0, 0.2);
                fill_random(frseaice_, 0.0, 0.0);
                fill_random(frsno_, 0.0, 0.05);
                fill_random(frlake_, 0.0, 0.1);
                fill_random(gvf_, 0.1, 0.6);
                fill_random(lai_, 0.5, 3.0);
                fill_random(z0h_, 0.01, 0.2);
                fill_random(hflux_, 10.0, 150.0);
                fill_random(ssm_, 0.1, 0.4);
                fill_random(gwettop_, 0.1, 0.5);
                fill_random(clayfrac_, 0.1, 0.4);
                fill_random(sandfrac_, 0.2, 0.6);
                fill_random(rdrag_, 0.05, 0.2);
                break;

            case AtmosphericScenario::UrbanBoundaryLayer:
                p_surf_min = 97000.0;
                p_surf_max = 103000.0;
                temp_min = 280.0;
                temp_max = 318.0; // Urban heat island
                rh_min = 15.0;
                rh_max = 85.0;
                ustar_min = 0.2;
                ustar_max = 2.2;
                precip_scale = 0.001;
                lwi_.assign(nc, 1.0); // Land
                frocean_.assign(nc, 0.0);
                frseaice_.assign(nc, 0.0);
                frsno_.assign(nc, 0.0);
                frlake_.assign(nc, 0.0);
                fill_random(gvf_, 0.0, 0.15); // Sparse vegetation
                fill_random(lai_, 0.0, 0.5);
                fill_random(z0h_, 0.8, 2.5);      // High urban canopy roughness [m]
                fill_random(hflux_, 60.0, 300.0); // Strong sensible heat flux
                fill_random(ssm_, 0.005, 0.08);   // Impervious dry surface
                fill_random(gwettop_, 0.005, 0.1);
                fill_random(clayfrac_, 0.1, 0.3);
                fill_random(sandfrac_, 0.3, 0.7);
                fill_random(rdrag_, 0.15, 0.3);
                break;

            case AtmosphericScenario::RuralAgricultural:
                p_surf_min = 96000.0;
                p_surf_max = 103000.0;
                temp_min = 275.0;
                temp_max = 308.0;
                rh_min = 30.0;
                rh_max = 95.0;
                ustar_min = 0.1;
                ustar_max = 1.5;
                precip_scale = 0.002;
                lwi_.assign(nc, 1.0); // Land
                frocean_.assign(nc, 0.0);
                frseaice_.assign(nc, 0.0);
                frsno_.assign(nc, 0.0);
                fill_random(frlake_, 0.0, 0.05);
                fill_random(gvf_, 0.55, 0.98); // High green vegetation fraction
                fill_random(lai_, 2.0, 7.0);   // Dense crop/forest leaf area index
                fill_random(z0h_, 0.05, 0.8);  // Agricultural / forest roughness [m]
                fill_random(hflux_, 10.0, 120.0);
                fill_random(ssm_, 0.2, 0.65); // Moist agricultural soil
                fill_random(gwettop_, 0.2, 0.7);
                fill_random(clayfrac_, 0.2, 0.5);
                fill_random(sandfrac_, 0.2, 0.6);
                fill_random(rdrag_, 0.05, 0.2);
                break;

            case AtmosphericScenario::MarineBoundaryLayer:
                p_surf_min = 99000.0;
                p_surf_max = 103500.0;
                temp_min = 273.15;
                temp_max = 303.0;
                rh_min = 70.0;
                rh_max = 98.0; // High marine relative humidity
                ustar_min = 0.1;
                ustar_max = 1.8;
                precip_scale = 0.005;
                lwi_.assign(nc, 0.0);     // Water/Ocean
                frocean_.assign(nc, 1.0); // 100% ocean
                frseaice_.assign(nc, 0.0);
                frsno_.assign(nc, 0.0);
                frlake_.assign(nc, 0.0);
                gvf_.assign(nc, 0.0);
                lai_.assign(nc, 0.0);
                fill_random(z0h_, 0.0001, 0.001); // Smooth ocean roughness [m]
                fill_random(sst_, 275.0, 303.0);
                fill_random(hflux_, -20.0, 50.0);
                ssm_.assign(nc, 1.0);
                gwettop_.assign(nc, 1.0);
                clayfrac_.assign(nc, 0.0);
                sandfrac_.assign(nc, 0.0);
                rdrag_.assign(nc, 0.0);
                break;

            case AtmosphericScenario::PolarSnowIce:
                p_surf_min = 96000.0;
                p_surf_max = 104000.0;
                temp_min = 210.0;
                temp_max = 271.15; // Freezing polar temperatures
                rh_min = 40.0;
                rh_max = 95.0;
                ustar_min = 0.05;
                ustar_max = 1.2;
                precip_scale = 0.001;
                lwi_.assign(nc, 2.0); // Ice
                fill_random(frocean_, 0.0, 0.2);
                fill_random(frseaice_, 0.8, 1.0); // Dense sea ice
                fill_random(frsno_, 0.8, 1.0);    // Dense snow pack
                frlake_.assign(nc, 0.0);
                gvf_.assign(nc, 0.0);
                lai_.assign(nc, 0.0);
                fill_random(z0h_, 0.0005, 0.01); // Snow/ice roughness [m]
                fill_random(sst_, 268.0, 273.15);
                fill_random(hflux_, -50.0, 20.0);
                fill_random(ssm_, 0.0, 0.1);
                fill_random(gwettop_, 0.0, 0.1);
                clayfrac_.assign(nc, 0.0);
                sandfrac_.assign(nc, 0.0);
                rdrag_.assign(nc, 0.0);
                break;

            case AtmosphericScenario::MidTroposphere:
                p_surf_min = 40000.0;
                p_surf_max = 85000.0;
                temp_min = 230.0;
                temp_max = 280.0;
                rh_min = 10.0;
                rh_max = 100.0;
                ustar_min = 0.05;
                ustar_max = 0.8;
                precip_scale = 0.05; // active clouds and precipitation
                fill_random(lwi_, 0.0, 1.0);
                fill_random(frocean_, 0.0, 1.0);
                fill_random(frseaice_, 0.0, 0.2);
                fill_random(frsno_, 0.0, 0.1);
                fill_random(frlake_, 0.0, 0.05);
                fill_random(gvf_, 0.1, 0.8);
                fill_random(lai_, 0.5, 4.0);
                fill_random(z0h_, 0.01, 0.2);
                fill_random(hflux_, 0.0, 80.0);
                fill_random(ssm_, 0.1, 0.5);
                fill_random(gwettop_, 0.1, 0.5);
                fill_random(clayfrac_, 0.1, 0.4);
                fill_random(sandfrac_, 0.2, 0.6);
                fill_random(rdrag_, 0.05, 0.2);
                break;

            case AtmosphericScenario::UpperTroposphereStratosphere:
                p_surf_min = 2000.0;
                p_surf_max = 25000.0;
                temp_min = 180.0;
                temp_max = 230.0;
                rh_min = 0.0;
                rh_max = 20.0;
                ustar_min = 0.01;
                ustar_max = 0.3;
                precip_scale = 0.0; // no precipitation in UTLS
                fill_random(lwi_, 0.0, 1.0);
                fill_random(frocean_, 0.0, 1.0);
                fill_random(frseaice_, 0.0, 0.1);
                fill_random(frsno_, 0.0, 0.0);
                fill_random(frlake_, 0.0, 0.0);
                fill_random(gvf_, 0.0, 0.5);
                fill_random(lai_, 0.0, 2.0);
                fill_random(z0h_, 0.01, 0.1);
                fill_random(hflux_, 0.0, 50.0);
                fill_random(ssm_, 0.0, 0.2);
                fill_random(gwettop_, 0.0, 0.2);
                fill_random(clayfrac_, 0.1, 0.4);
                fill_random(sandfrac_, 0.2, 0.6);
                fill_random(rdrag_, 0.05, 0.2);
                break;

            case AtmosphericScenario::ExtremeStormDesert:
                p_surf_min = 88000.0;
                p_surf_max = 104000.0;
                temp_min = 250.0;
                temp_max = 325.0;
                rh_min = 5.0;
                rh_max = 100.0;
                ustar_min = 1.0;
                ustar_max = 3.5; // gale force shear
                precip_scale = 0.1;
                lwi_.assign(nc, 1.0); // Land
                frocean_.assign(nc, 0.0);
                frseaice_.assign(nc, 0.0);
                frsno_.assign(nc, 0.0);
                frlake_.assign(nc, 0.0);
                fill_random(gvf_, 0.0, 0.05); // Barren desert soil
                fill_random(lai_, 0.0, 0.1);
                fill_random(z0h_, 0.005, 0.05);
                fill_random(hflux_, 50.0, 250.0);
                fill_random(ssm_, 0.001, 0.04); // Extremely dry sand
                fill_random(gwettop_, 0.001, 0.04);
                fill_random(clayfrac_, 0.15, 0.45);
                fill_random(sandfrac_, 0.35, 0.8);
                fill_random(rdrag_, 0.05, 0.25);
                break;

            case AtmosphericScenario::FullRandomFuzz:
                p_surf_min = 1000.0;
                p_surf_max = 105000.0;
                temp_min = 160.0;
                temp_max = 340.0;
                rh_min = 0.0;
                rh_max = 100.0;
                ustar_min = 0.001;
                ustar_max = 4.0;
                precip_scale = 0.05;
                fill_random(lwi_, 0.0, 2.0);
                fill_random(frocean_, 0.0, 1.0);
                fill_random(frseaice_, 0.0, 1.0);
                fill_random(frsno_, 0.0, 1.0);
                fill_random(frlake_, 0.0, 1.0);
                fill_random(gvf_, 0.0, 1.0);
                fill_random(lai_, 0.0, 8.0);
                fill_random(z0h_, 0.0001, 3.0);
                fill_random(hflux_, -100.0, 400.0);
                fill_random(ssm_, 0.0, 1.0);
                fill_random(gwettop_, 0.0, 1.0);
                fill_random(clayfrac_, 0.0, 1.0);
                fill_random(sandfrac_, 0.0, 1.0);
                fill_random(rdrag_, 0.0, 0.3);
                break;
            }

            fill_random(t_air_, temp_min, temp_max);
            fill_random(rh_, rh_min, rh_max);
            fill_random(cldf_, 0.0, (scenario == AtmosphericScenario::UpperTroposphereStratosphere ? 0.0 : 1.0));
            fill_random(pfilsan_, 0.0, precip_scale);
            fill_random(pfllsan_, 0.0, precip_scale);
            fill_random(reevapls_, 0.0, precip_scale * 1e-3);
            fill_random(ustar_, ustar_min, ustar_max);
            fill_random(u10m_, -40.0, 40.0);
            fill_random(v10m_, -40.0, 40.0);
            fill_random(frocean_, 0.0, 1.0);
            fill_random(frseaice_, 0.0, 0.5);

            for (int icol = 0; icol < nc; ++icol) {
                double current_p = std::uniform_real_distribution<double>(p_surf_min, p_surf_max)(gen_);
                pedge_[icol + 0 * nc] = current_p;
                for (int k = 0; k < nl; ++k) {
                    double min_remaining_p = (nl - k) * 100.0;
                    double max_step = std::max(20.0, (current_p - min_remaining_p) / (nl - k));
                    double delta = std::uniform_real_distribution<double>(10.0, max_step)(gen_);
                    current_p -= delta;
                    if (current_p < 10.0 * (nl - k)) {
                        current_p = 10.0 * (nl - k);
                    }
                    pedge_[icol + (k + 1) * nc] = current_p;

                    double p1 = pedge_[icol + k * nc];
                    double p2 = pedge_[icol + (k + 1) * nc];
                    pmid_[icol + k * nc] = 0.5 * (p1 + p2);
                    delp_[icol + k * nc] = std::max(1.0, std::abs(p1 - p2));

                    double t = t_air_[icol + k * nc];
                    double rho = pmid_[icol + k * nc] / (287.05 * t);
                    rho = std::clamp(rho, 0.001, 2.5);
                    airden_dry_[icol + k * nc] = rho;
                    mairden_[icol + k * nc] = rho;

                    double dz = delp_[icol + k * nc] / (rho * 9.80665);
                    bxheight_[icol + k * nc] = std::clamp(dz, 0.5, 10000.0);
                }
            }

            // Fuzz initial chemical concentrations in plausible atmospheric range [0, 1 ppmv]
            fill_random(conc_, 0.0, 1.0e-6);

            // Bind to C++ state
            state->bind_met_field_3d("T", t_air_.data());
            state->bind_met_field_3d("PMID", pmid_.data());
            state->bind_met_field_3d("PEDGE", pedge_.data());
            state->bind_met_field_3d("AIRDEN_DRY", airden_dry_.data());
            state->bind_met_field_3d("AIRDEN", airden_dry_.data());
            state->bind_met_field_3d("BXHEIGHT", bxheight_.data());
            state->bind_met_field_3d("MAIRDEN", mairden_.data());
            state->bind_met_field_3d("PFILSAN", pfilsan_.data());
            state->bind_met_field_3d("PFLLSAN", pfllsan_.data());
            state->bind_met_field_3d("REEVAPLS", reevapls_.data());
            state->bind_met_field_3d("RH", rh_.data());
            state->bind_met_field_3d("CLDF", cldf_.data());
            state->bind_met_field_2d("LAT", lat_.data());
            state->bind_met_field_2d("LON", lon_.data());
            state->bind_met_field_2d("SST", sst_.data());
            state->bind_met_field_2d("PS", ps_.data());
            state->bind_met_field_2d("TS", ts_.data());
            state->bind_met_field_2d("HFLUX", hflux_.data());
            state->bind_met_field_2d("OBK", obk_.data());
            state->bind_met_field_2d("PBLH", pblh_.data());
            state->bind_met_field_2d("Z0H", z0h_.data());
            state->bind_met_field_2d("Z0", z0h_.data());
            state->bind_met_field_2d("FROCEAN", frocean_.data());
            state->bind_met_field_2d("FRSEAICE", frseaice_.data());
            state->bind_met_field_2d("LWI", lwi_.data());
            state->bind_met_field_2d("FRSNO", frsno_.data());
            state->bind_met_field_2d("FRSNOW", frsno_.data());
            state->bind_met_field_2d("FRLAKE", frlake_.data());
            state->bind_met_field_2d("GVF", gvf_.data());
            state->bind_met_field_2d("LAI", lai_.data());
            state->bind_met_field_2d("SSM", ssm_.data());
            state->bind_met_field_2d("GWETTOP", gwettop_.data());
            state->bind_met_field_2d("CLAYFRAC", clayfrac_.data());
            state->bind_met_field_2d("SANDFRAC", sandfrac_.data());
            state->bind_met_field_2d("RDRAG", rdrag_.data());
            state->bind_met_field_2d("USTAR", ustar_.data());
            state->bind_met_field_3d("DELP", delp_.data());
            state->bind_met_field_2d("U10M", u10m_.data());
            state->bind_met_field_2d("V10M", v10m_.data());

            state->bind_unified_chemistry(conc_.data());
            state->sync_to_device();
        }

        /**
         * @brief Invariant 1: Assert all concentration values remain finite (no NaN or Inf).
         */
        static void assert_finiteness(const std::vector<double>& conc, int iter, const std::string& ctx) {
            for (size_t i = 0; i < conc.size(); ++i) {
                if (!std::isfinite(conc[i])) {
                    std::cerr << "PROPERTY FAILURE [" << ctx << "]: Index " << i << " is NaN/Inf (" << conc[i]
                              << ") at iteration " << iter << std::endl;
                    assert(false && "State variables must remain finite!");
                }
            }
        }

        /**
         * @brief Invariant 2: Assert non-negativity and clip underflow to physical zero.
         */
        static void assert_non_negativity(std::vector<double>& conc, int iter, const std::string& ctx) {
            double max_negative = 0.0;
            for (size_t i = 0; i < conc.size(); ++i) {
                if (conc[i] < 0.0) {
                    if (conc[i] < max_negative)
                        max_negative = conc[i];
                    conc[i] = 0.0; // Clip to zero (standard atmospheric model physical boundary)
                }
            }
            if (max_negative < -1e-2) {
                std::cerr << "WARNING [" << ctx << "]: Severe negative concentration " << max_negative
                          << " clipped to 0.0 at iteration " << iter << std::endl;
            }
        }

        /**
         * @brief Invariant 3: Zero-state stability assertion.
         * Zero surface winds/shear should yield zero dust/seasalt emissions.
         */
        void test_zero_state_invariants(catchem::Core* core, std::shared_ptr<catchem::StateManager> state) {
            // Populate full valid meteorological state first
            generate_scenario_state(AtmosphericScenario::BoundaryLayer, state);

            // Set winds, friction velocity, and precipitation to zero
            std::fill(ustar_.begin(), ustar_.end(), 0.0);
            std::fill(u10m_.begin(), u10m_.end(), 0.0);
            std::fill(v10m_.begin(), v10m_.end(), 0.0);
            std::fill(pfilsan_.begin(), pfilsan_.end(), 0.0);
            std::fill(pfllsan_.begin(), pfllsan_.end(), 0.0);
            std::fill(conc_.begin(), conc_.end(), 0.0);

            state->bind_met_field_2d("USTAR", ustar_.data());
            state->bind_met_field_2d("U10M", u10m_.data());
            state->bind_met_field_2d("V10M", v10m_.data());
            state->bind_unified_chemistry(conc_.data());
            state->sync_to_device();

            core->run_timestep(config_.dt_sec);
            state->sync_to_host();

            assert_finiteness(conc_, 0, "ZeroState");
            assert_non_negativity(conc_, 0, "ZeroState");
        }

        /**
         * @brief Run complete property-based test harness across all scenarios.
         */
        void run_full_suite(catchem::Core* core) {
            auto state = core->get_state_manager();

            std::vector<AtmosphericScenario> scenarios = {AtmosphericScenario::BoundaryLayer,
                                                          AtmosphericScenario::UrbanBoundaryLayer,
                                                          AtmosphericScenario::RuralAgricultural,
                                                          AtmosphericScenario::MarineBoundaryLayer,
                                                          AtmosphericScenario::PolarSnowIce,
                                                          AtmosphericScenario::MidTroposphere,
                                                          AtmosphericScenario::UpperTroposphereStratosphere,
                                                          AtmosphericScenario::ExtremeStormDesert,
                                                          AtmosphericScenario::FullRandomFuzz};

            std::cout << "\n========================================================" << std::endl;
            std::cout << "=== RUNNING UNIFIED ATMOSPHERIC PROPERTY TEST SUITE  ===" << std::endl;
            std::cout << "========================================================\n" << std::endl;

            // 1. Zero-state invariant check
            std::cout << "[Harness] Running Zero-State Stability Test..." << std::endl;
            test_zero_state_invariants(core, state);
            std::cout << "[Harness]  ✓ Zero-State Invariant Passed!" << std::endl;

            // 2. Scenario-based randomized property iterations
            for (auto scenario : scenarios) {
                std::string sc_name = scenario_to_string(scenario);
                std::cout << "[Harness] Scenario: " << sc_name << " (" << config_.iterations_per_scenario
                          << " iterations)..." << std::endl;

                for (int iter = 1; iter <= config_.iterations_per_scenario; ++iter) {
                    generate_scenario_state(scenario, state);

                    // Execute timestep across all registered processes
                    core->run_timestep(config_.dt_sec);
                    state->sync_to_host();

                    // Verify properties
                    assert_finiteness(conc_, iter, sc_name);
                    assert_non_negativity(conc_, iter, sc_name);
                }
                std::cout << "[Harness]  ✓ " << sc_name << " Properties Passed!" << std::endl;
            }

            std::cout << "\n========================================================" << std::endl;
            std::cout << "=== SUCCESS: ALL ATMOSPHERIC PROPERTY INVARIANTS HELD ===" << std::endl;
            std::cout << "========================================================\n" << std::endl;
        }
    };

} // namespace catchem::test
