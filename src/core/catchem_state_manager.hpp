#pragma once
#include "catchem_chem_state.hpp"
#include "catchem_config_manager.hpp"
#include "catchem_constants.hpp"
#include "catchem_interop_field.hpp"
#include "catchem_met_state.hpp"
#include "catchem_met_utilities.hpp"
#include "catchem_time_state.hpp"
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>
#include <yaml-cpp/yaml.h>

namespace catchem {

    class DiagnosticManager;

    class StateManager {
    public:
        int n_cols;
        int n_levels;
        int n_species;

        std::string config_file_path;
        std::string trace_id;

        std::shared_ptr<ConfigManager> config_mgr;
        std::shared_ptr<DiagnosticManager> diag_mgr;

        std::unordered_map<std::string, std::shared_ptr<InteropField<double, 1>>> fields_1d;
        std::unordered_map<std::string, std::shared_ptr<InteropField<double, 2>>> fields_2d;
        std::unordered_map<std::string, std::shared_ptr<InteropField<double, 3>>> fields_3d;

        // Structured sub-states
        MetState met;
        ChemState chem;
        TimeState time;

        StateManager(int nc, int nl, int ns);

        void load_species_config(const std::string& filename) { chem.load_species_config(filename); }

        void bind_met_field_2d(const std::string& name, double* ptr) {
            auto field = std::make_shared<InteropField<double, 2>>(ptr, std::vector<int>{n_cols, n_levels});
            if (name == "PS")
                met.PS = field;
            else if (name == "TS")
                met.TS = field;
            else if (name == "PBLH")
                met.PBLH = field;
            else if (name == "USTAR")
                met.USTAR = field;
            else if (name == "HFLUX")
                met.HFLUX = field;
            else if (name == "OBK")
                met.OBK = field;
            else if (name == "LAT")
                met.LAT = field;
            else if (name == "LON")
                met.LON = field;
            met.fields_2d[name] = field;
        }

        void bind_met_field_3d(const std::string& name, double* ptr) {
            int nl = (name == "PEDGE" || name == "PFILSAN" || name == "PFLLSAN") ? n_levels + 1 : n_levels;
            auto field = std::make_shared<InteropField<double, 3>>(
                ptr, std::vector<int>{n_cols, nl, 1}); // Using 1 for single-field layout
            if (name == "T")
                met.T = field;
            else if (name == "QV")
                met.QV = field;
            else if (name == "RH")
                met.RH = field;
            else if (name == "PMID")
                met.PMID = field;
            else if (name == "PEDGE")
                met.PEDGE = field;
            else if (name == "AIRDEN")
                met.AIRDEN = field;
            else if (name == "AIRDEN_DRY")
                met.AIRDEN_DRY = field;
            else if (name == "BXHEIGHT")
                met.BXHEIGHT = field;
            met.fields_3d[name] = field;
        }

        void bind_unified_chemistry(double* ptr) {
            chem.conc = std::make_shared<InteropField<double, 3>>(ptr, std::vector<int>{n_cols, n_levels, n_species});
        }

        /**
         * @brief Derives layer thicknesses (BXHEIGHT / dz) hydrostatically.
         *
         * Computes the vertical distance of each grid cell based on pressure edges, temperature, and moisture content.
         * @throws std::runtime_error If required input fields (PEDGE, T, QV, BXHEIGHT) are not bound.
         */
        void derive_bxheight() {
            if (!met.PEDGE)
                throw std::runtime_error("derive_bxheight failed: PEDGE field is not bound.");
            if (!met.T)
                throw std::runtime_error("derive_bxheight failed: T field is not bound.");
            if (!met.QV)
                throw std::runtime_error("derive_bxheight failed: QV field is not bound.");
            if (!met.BXHEIGHT)
                throw std::runtime_error("derive_bxheight failed: BXHEIGHT field is not bound.");

            int nc = n_cols;
            int nl = n_levels;

            auto pedge = met.PEDGE->view();
            auto temp = met.T->view();
            auto qv = met.QV->view();
            auto bxheight = met.BXHEIGHT->view();

            Kokkos::parallel_for(
                "derive_bxheight_kernel",
                Kokkos::MDRangePolicy<Kokkos::DefaultExecutionSpace, Kokkos::Rank<2>>({0, 0}, {nc, nl}),
                KOKKOS_LAMBDA(int icol, int ilev) {
                    double p_lower = pedge(icol, ilev, 0);
                    double p_upper = pedge(icol, ilev + 1, 0);

                    if (p_upper > 0.0) {
                        double virtual_t = met_utilities::virtual_temperature(temp(icol, ilev, 0), qv(icol, ilev, 0));
                        bxheight(icol, ilev, 0) =
                            (constants::RD / constants::G0) * virtual_t * std::log(p_lower / p_upper);
                    } else {
                        bxheight(icol, ilev, 0) = 0.0;
                    }
                });
        }

        /**
         * @brief Derives dry air density (AIRDEN_DRY) using the Ideal Gas Law.
         *
         * Calculates dry air mass densities based on mid-point pressures, temperature, and specific humidity.
         * @throws std::runtime_error If required input fields (PMID, T, QV, AIRDEN_DRY) are not bound.
         */
        void derive_airden_dry() {
            if (!met.PMID)
                throw std::runtime_error("derive_airden_dry failed: PMID field is not bound.");
            if (!met.T)
                throw std::runtime_error("derive_airden_dry failed: T field is not bound.");
            if (!met.QV)
                throw std::runtime_error("derive_airden_dry failed: QV field is not bound.");
            if (!met.AIRDEN_DRY)
                throw std::runtime_error("derive_airden_dry failed: AIRDEN_DRY field is not bound.");

            int nc = n_cols;
            int nl = n_levels;

            auto pmid = met.PMID->view();
            auto temp = met.T->view();
            auto qv = met.QV->view();
            auto airden_dry = met.AIRDEN_DRY->view();

            Kokkos::parallel_for(
                "derive_airden_dry_kernel",
                Kokkos::MDRangePolicy<Kokkos::DefaultExecutionSpace, Kokkos::Rank<2>>({0, 0}, {nc, nl}),
                KOKKOS_LAMBDA(int icol, int ilev) {
                    double q = qv(icol, ilev, 0);
                    double avgw = (constants::AIR_MW / constants::H2O_MW) * q / (1.0 - q);
                    double xh2o = avgw / (1.0 + avgw);

                    double p_dry = pmid(icol, ilev, 0) * (1.0 - xh2o);
                    airden_dry(icol, ilev, 0) = p_dry / (constants::RD * temp(icol, ilev, 0));
                });
        }

        void sync_to_device() {
            for (auto& [k, v] : fields_1d)
                v->sync_to_device();
            for (auto& [k, v] : fields_2d)
                v->sync_to_device();
            for (auto& [k, v] : fields_3d)
                v->sync_to_device();
            for (auto& [k, v] : met.fields_2d)
                v->sync_to_device();
            for (auto& [k, v] : met.fields_3d)
                v->sync_to_device();
            if (chem.conc)
                chem.conc->sync_to_device();
        }

        void sync_to_host() {
            for (auto& [k, v] : fields_1d)
                v->sync_to_host();
            for (auto& [k, v] : fields_2d)
                v->sync_to_host();
            for (auto& [k, v] : fields_3d)
                v->sync_to_host();
            for (auto& [k, v] : met.fields_2d)
                v->sync_to_host();
            for (auto& [k, v] : met.fields_3d)
                v->sync_to_host();
            if (chem.conc)
                chem.conc->sync_to_host();
        }

        void bind_field_1d(const std::string& name, double* ptr) {
            fields_1d[name] = std::make_shared<InteropField<double, 1>>(ptr, std::vector<int>{n_cols});
        }

        void bind_field_2d(const std::string& name, double* ptr) {
            fields_2d[name] = std::make_shared<InteropField<double, 2>>(ptr, std::vector<int>{n_cols, n_levels});
        }

        void bind_field_3d(const std::string& name, double* ptr) {
            fields_3d[name] =
                std::make_shared<InteropField<double, 3>>(ptr, std::vector<int>{n_cols, n_levels, n_species});
        }

        double* get_host_pointer_1d(const std::string& name) {
            if (fields_1d.find(name) == fields_1d.end())
                return nullptr;
            return fields_1d.at(name)->host_view.data();
        }

        double* get_host_pointer_2d(const std::string& name) {
            if (fields_2d.find(name) != fields_2d.end())
                return fields_2d.at(name)->host_view.data();
            if (met.fields_2d.find(name) != met.fields_2d.end())
                return met.fields_2d.at(name)->host_view.data();
            return nullptr;
        }

        double* get_host_pointer_3d(const std::string& name) {
            if (fields_3d.find(name) != fields_3d.end())
                return fields_3d.at(name)->host_view.data();
            if (met.fields_3d.find(name) != met.fields_3d.end())
                return met.fields_3d.at(name)->host_view.data();
            return nullptr;
        }
    };

} // namespace catchem
