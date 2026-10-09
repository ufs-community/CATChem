#pragma once
#include "catchem_process_interface.hpp"
#include <string>
#include <vector>

namespace catchem {

    class SeaSaltProcess : public ProcessInterface {
    private:
        std::string active_scheme;
        bool diagnostics_enabled;
        std::vector<int> diagnostic_species_id;

        // Per-scheme tuning options read from processes.seasalt.<scheme>.*
        // during init() and forwarded to the science bridge on every run().
        // Defaults mirror the SeaSaltScheme*Config types in SeaSaltCommon_Mod.F90.
        double gong97_scale_factor = 1.0;
        bool gong97_weibull_flag = false;
        double gong03_scale_factor = 1.0;
        bool gong03_weibull_flag = false;
        double geos12_scale_factor = 1.0;
        bool geos12_weibull_flag = false;

    public:
        SeaSaltProcess();
        std::string get_name() const override { return "seasalt"; }
        ProcessContract get_contract() const override;
        void prepare_inputs(std::shared_ptr<StateManager> state) override;
        void init(std::shared_ptr<StateManager> state) override;
        void run(std::shared_ptr<StateManager> state) override;
        void finalize() override {}
    };

} // namespace catchem
