#pragma once

#include "catchem_process_interface.hpp"
#include "catchem_state_manager.hpp"
#include <memory>
#include <string>
#include <vector>

namespace catchem {

    class DustProcess : public ProcessInterface {
    public:
        DustProcess();
        ~DustProcess() override = default;

        std::string get_name() const override { return "dust"; }
        ProcessContract get_contract() const override;
        void prepare_inputs(std::shared_ptr<StateManager> state) override;
        void init(std::shared_ptr<StateManager> state) override;
        void run(std::shared_ptr<StateManager> state) override;
        void finalize() override;

    private:
        std::string active_scheme_;
        bool diagnostics_enabled_;
        std::vector<int> diagnostic_species_id_;

        // Fengsha scheme tuning options, read from processes.dust.fengsha.*
        // during init() and forwarded to the science bridge on every run().
        // Defaults mirror DustSchemeFENGSHAConfig in DustCommon_Mod.F90.
        double fengsha_alpha_ = 0.2;
        double fengsha_gamma_ = 1.0;
        double fengsha_drylimit_factor_ = 1.0;
        double fengsha_moist_correction_factor_ = 1.0;
        double fengsha_kvhmax_ = 0.0002;
        int fengsha_drag_option_ = 1;
        int fengsha_horizflux_option_ = 1;
        int fengsha_moist_option_ = 1;
        int fengsha_distribution_option_ = 1;

        // Ginoux scheme tuning options, read from processes.dust.ginoux.*
        // during init().  Ch_DU holds one multiplier per dust size bin and
        // must stay the same length as DustSchemeGINOUXConfig%Ch_DU.
        std::vector<double> ginoux_ch_du_ = std::vector<double>(5, 1.0);
    };

} // namespace catchem
