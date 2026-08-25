#pragma once
#include "catchem_process_interface.hpp"
#include <string>
#include <vector>

namespace catchem {

    class SO4chemProcess : public ProcessInterface {
    private:
        std::string active_scheme;
        bool diagnostics_enabled;
        std::vector<int> diagnostic_species_id;

        // Persistent column states
        std::vector<char> firsttime;
        std::vector<int> nymd_last;
        std::vector<int> nhms_last_recycle;
        std::vector<double> xh2o2_init;
        std::vector<double> pso4_so2;
        std::vector<double> pso4_g_so2;
        std::vector<double> pso4_aq_so2;
        std::vector<double> pso2_dms;
        std::vector<double> dms_flux;

    public:
        SO4chemProcess();
        std::string get_name() const override { return "so4chem"; }
        void init(std::shared_ptr<StateManager> state) override;
        void run(std::shared_ptr<StateManager> state) override;
        void finalize() override {}
    };

} // namespace catchem
