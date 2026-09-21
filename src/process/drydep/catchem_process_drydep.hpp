#pragma once
#include "catchem_process_interface.hpp"
#include <string>
#include <vector>

namespace catchem {

    class DryDepProcess : public ProcessInterface {
    private:
        std::string gas_scheme;
        std::string aero_scheme;
        bool diagnostics_enabled;
        std::vector<int> diagnostic_species_id;

        // Per-scheme tuning options read from processes.drydep.<scheme>.*
        // during init() and forwarded to the science bridge on every run().
        // Defaults mirror the DryDepScheme*Config types in DryDepCommon_Mod.F90.
        double wesely_scale_factor = 1.0;
        bool wesely_co2_effect = true;
        double wesely_co2_level = 600.0;
        double wesely_co2_reference = 380.0;
        double gocart_scale_factor = 1.0;
        bool gocart_resuspension = false;
        bool gocart_dust_resuspension_only = true;
        double zhang_scale_factor = 1.0;

    public:
        DryDepProcess();
        std::string get_name() const override { return "drydep"; }
        ProcessContract get_contract() const override;
        void prepare_inputs(std::shared_ptr<StateManager> state) override;
        void init(std::shared_ptr<StateManager> state) override;
        void run(std::shared_ptr<StateManager> state) override;
        void finalize() override {}
    };

} // namespace catchem
