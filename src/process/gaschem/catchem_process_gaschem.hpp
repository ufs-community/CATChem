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
        // Whether this PET has already warned about a floored input and about a
        // floored solver output. Kept apart so that an early warning about host
        // noise cannot hide a later negative from the solver.
        bool input_clamp_warned_ = false;
        bool output_clamp_warned_ = false;
        std::size_t vector_size_ = 1;
        std::size_t n_micm_species_ = 0;
        std::size_t n_rate_params_ = 0;
        std::vector<SpeciesMapping> species_mappings_;
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
        double min_negative_output() const noexcept { return min_negative_output_; }
        bool input_clamp_warned() const noexcept { return input_clamp_warned_; }
        bool output_clamp_warned() const noexcept { return output_clamp_warned_; }

        void init(std::shared_ptr<StateManager> state) override;
        void run(std::shared_ptr<StateManager> state) override;
        void finalize() override;

    protected:
        // Test seam. run() calls this once per batch, after a solve that converged
        // and before the solver-output guard, with the MICM concentration array
        // (mol/m3, MICM ordering, padding cells included). It does nothing here.
        // A test overrides it to make the solver return a negative value, which
        // the mechanisms in the test suite never do on their own.
        virtual void adjust_solver_output_for_testing(std::vector<double>& /*concentrations*/, int /*batch*/) {}
    };

} // namespace catchem

extern "C" {
void catchem_register_gaschem_cpp();
}
