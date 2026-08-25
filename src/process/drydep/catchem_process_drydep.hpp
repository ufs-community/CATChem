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

    public:
        DryDepProcess();
        std::string get_name() const override { return "drydep"; }
        void init(std::shared_ptr<StateManager> state) override;
        void run(std::shared_ptr<StateManager> state) override;
        void finalize() override {}
    };

} // namespace catchem
