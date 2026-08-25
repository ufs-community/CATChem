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

    public:
        WetDepProcess();
        std::string get_name() const override { return "wetdep"; }
        void init(std::shared_ptr<StateManager> state) override;
        void run(std::shared_ptr<StateManager> state) override;
        void finalize() override {}
    };

} // namespace catchem
