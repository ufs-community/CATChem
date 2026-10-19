#pragma once
#include "catchem_process_interface.hpp"
#include <string>
#include <vector>

namespace catchem {

    class SO4chemProcess : public ProcessInterface {
    private:
        std::string active_scheme_;
        bool diagnostics_enabled_;
        std::vector<int> diagnostic_species_id_;

        // GOCART scheme options staged from the runtime configuration
        // (processes/so4chem/gocart/*) and forwarded to the science bridge.
        bool gocart_update_so2_ = true;
        bool gocart_do_drydep_ = false;

        // Persistent column states
        // Stored as C int rather than char/bool so the Fortran bridge has an
        // exact, alias-safe integer(c_int) ABI.  std::vector<bool> is packed
        // and vector<char> cannot legally be reinterpreted as bool*.
        std::vector<int> firsttime_;
        std::vector<int> nymd_last_;
        std::vector<int> nhms_last_recycle_;
        std::vector<double> xh2o2_init_;
        std::vector<double> pso4_g_so2_;
        std::vector<double> pso4_aq_so2_;
        std::vector<double> pso2_dms_;
        std::vector<double> dms_flux_;

    public:
        SO4chemProcess();
        std::string get_name() const override { return "so4chem"; }
        ProcessContract get_contract() const override;
        void prepare_inputs(std::shared_ptr<StateManager> state) override;
        void init(std::shared_ptr<StateManager> state) override;
        void run(std::shared_ptr<StateManager> state) override;
        void finalize() override {}
    };

} // namespace catchem
