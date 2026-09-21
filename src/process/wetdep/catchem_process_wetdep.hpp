#pragma once
#include "catchem_process_interface.hpp"
#include <string>
#include <vector>

namespace catchem {

    class WetDepProcess : public ProcessInterface {
    private:
        std::string active_scheme;
        bool diagnostics_enabled;
        std::vector<int> diagnostic_species_id;

        // Jacob scheme tuning options read from processes.wetdep.jacob.*
        // during init() and forwarded to the science bridge on every run().
        // Defaults mirror WetDepSchemeJACOBConfig in WetDepCommon_Mod.F90.
        double jacob_scale_factor = 1.0;
        double jacob_radius_threshold = 1.0;
        bool jacob_so4_gocart_resusp = true;
        double jacob_so4_washout_eff = 1.0;

    public:
        WetDepProcess();
        std::string get_name() const override { return "wetdep"; }
        ProcessContract get_contract() const override;
        void prepare_inputs(std::shared_ptr<StateManager> state) override;
        void init(std::shared_ptr<StateManager> state) override;
        void run(std::shared_ptr<StateManager> state) override;
        void finalize() override {}
    };

} // namespace catchem
