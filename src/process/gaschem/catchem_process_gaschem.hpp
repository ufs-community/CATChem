// src/process/gaschem/catchem_process_gaschem.hpp
#pragma once
#include "catchem_process_interface.hpp"
#include <cstddef>
#include <memory>
#include <musica/micm/micm.hpp>
#include <musica/micm/state.hpp>
#include <string>
#include <vector>

namespace catchem {

    class GasChemProcess : public ProcessInterface {
    private:
        struct SpeciesMapping {
            int catchem_idx = 0;
            std::size_t micm_idx = 0;
            std::string name;
        };

        enum class RateParamKind { Photolysis, Loss };

        struct RateParamMapping {
            std::size_t micm_param_idx = 0;
            RateParamKind kind = RateParamKind::Photolysis;
            std::string diag_name;
        };

        // What run() does with a negative concentration, on input from the host
        // or on output from the solver. Non-finite values throw under both.
        enum class NegativeConcentrationPolicy { Reject, Clamp };

        std::string config_dir_;
        // Optional explicit single-file MUSICA mechanism (v1 format). When set, it
        // takes precedence over config_dir and is passed directly to MICM.
        std::string config_file_;
        std::size_t columns_per_batch_ = 1;
        NegativeConcentrationPolicy negative_concentration_policy_ = NegativeConcentrationPolicy::Reject;
        // Cells floored to zero since init(), and the most negative value seen.
        std::size_t negative_inputs_clamped_ = 0;
        std::size_t negative_outputs_clamped_ = 0;
        double min_negative_input_ = 0.0;
        double min_negative_output_ = 0.0;
        bool clamp_warning_issued_ = false;
        std::size_t vector_size_ = 1;
        std::size_t n_micm_species_ = 0;
        std::size_t n_rate_params_ = 0;
        std::vector<SpeciesMapping> species_mappings_;
        // MICM variables that no CATChem species or alias resolves to, sorted by name.
        std::vector<std::string> unmapped_micm_species_;
        std::vector<RateParamMapping> rate_param_mappings_;
        std::unique_ptr<musica::MICM> micm_instance_;
        std::unique_ptr<musica::State> micm_state_;
        bool initialized_ = false;

    public:
        GasChemProcess();
        ~GasChemProcess() override;

        std::string get_name() const override { return std::string(ProcessNames::GasChem); }
        ProcessContract get_contract() const override;

        std::size_t columns_per_batch() const noexcept { return columns_per_batch_; }
        std::size_t micm_grid_cell_count() const { return micm_state_ ? micm_state_->NumberOfGridCells() : 0; }
        std::size_t negative_inputs_clamped() const noexcept { return negative_inputs_clamped_; }
        std::size_t negative_outputs_clamped() const noexcept { return negative_outputs_clamped_; }
        double min_negative_input() const noexcept { return min_negative_input_; }

        void init(std::shared_ptr<StateManager> state) override;
        void run(std::shared_ptr<StateManager> state) override;
        void finalize() override;
    };

} // namespace catchem

extern "C" {
void catchem_register_gaschem_cpp();
}
