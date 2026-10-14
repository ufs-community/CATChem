#pragma once

#include "catchem_process_interface.hpp"
#include "catchem_state_manager.hpp"
#include <memory>
#include <string>
#include <vector>

namespace catchem {

    class CarbChemProcess : public ProcessInterface {
    public:
        CarbChemProcess();
        ~CarbChemProcess() override = default;

        std::string get_name() const override { return "carbchem"; }
        ProcessContract get_contract() const override;
        void prepare_inputs(std::shared_ptr<StateManager> state) override;
        void init(std::shared_ptr<StateManager> state) override;
        void run(std::shared_ptr<StateManager> state) override;
        void finalize() override;

    private:
        std::string active_scheme;
        bool diagnostics_enabled;
        std::vector<int> diagnostic_species_id;

        // GOCART scheme options staged from the runtime configuration
        // (processes/carbchem/gocart/*) and forwarded to the science bridge.
        double gocart_time_days = 2.5;
    };

} // namespace catchem
