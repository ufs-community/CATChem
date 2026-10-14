// src/process/gaschem/catchem_process_gaschem.hpp
#pragma once
#include "catchem_process_interface.hpp"
#include <memory>
#include <musica/micm/micm.hpp>
#include <musica/micm/state.hpp>
#include <string>

namespace catchem {

    class GasChemProcess : public ProcessInterface {
    private:
        std::string config_dir;
        // Optional explicit single-file MUSICA mechanism (v1 format). When set, it
        // takes precedence over config_dir and is passed directly to MICM.
        std::string config_file;
        std::unique_ptr<musica::MICM> micm_instance;
        std::unique_ptr<musica::State> micm_state;
        bool initialized = false;

    public:
        GasChemProcess();
        ~GasChemProcess() override;

        std::string get_name() const override { return std::string(ProcessNames::GasChem); }
        ProcessContract get_contract() const override;

        void init(std::shared_ptr<StateManager> state) override;
        void run(std::shared_ptr<StateManager> state) override;
        void finalize() override;
    };

} // namespace catchem

extern "C" {
void catchem_register_gaschem_cpp();
}
