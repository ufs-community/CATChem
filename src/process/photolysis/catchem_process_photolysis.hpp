// src/process/photolysis/catchem_process_photolysis.hpp
#pragma once
#include "catchem_process_interface.hpp"
#include <memory>
#include <musica/tuvx/tuvx.hpp>
#include <string>
#include <unordered_set>

namespace catchem {

    class PhotolysisProcess : public ProcessInterface {
    private:
        std::string config_path;
        musica::TUVX* tuvx_instance = nullptr;
        musica::Mappings photo_mappings;

        // Keep raw map pointers alive
        musica::GridMap* grids = nullptr;
        musica::ProfileMap* profiles = nullptr;
        musica::RadiatorMap* radiators = nullptr;

        void register_profile_if_missing(const StateManager* state,
                                         const std::unordered_set<std::string>& config_defined_profiles,
                                         const char* name, const char* units, musica::Grid* grid, double default_val,
                                         std::size_t num_vals, musica::Error* err);

    public:
        PhotolysisProcess();
        ~PhotolysisProcess() override;

        std::string get_name() const override { return std::string(ProcessNames::Photolysis); }

        void init(std::shared_ptr<StateManager> state) override;
        void run(std::shared_ptr<StateManager> state) override;
        void finalize() override;
    };

} // namespace catchem
