#include "catchem_core.hpp"
#include "catchem_diagnostic_manager.hpp"
#include "catchem_process_registry.hpp"
#include <cassert>
#include <vector>

class ExtensionProcess final : public catchem::ProcessInterface {
public:
    explicit ExtensionProcess(bool& finalized) : finalized_(finalized) {}
    std::string get_name() const override { return "test-extension"; }
    catchem::ProcessContract get_contract() const override {
        return catchem::make_contract(get_name(), {catchem::host_field_3d("T", "K"), catchem::host_concentration()},
                                      {{"extension.marker", "", true}},
                                      {{"extension_value",
                                        "1",
                                        {catchem::SemanticAxis::Column, catchem::SemanticAxis::Singleton},
                                        catchem::PersistencePolicy::Timestep}});
    }
    void init(std::shared_ptr<catchem::StateManager> state) override {
        state->diagnostic_manager()->register_field("extension_value", "extension output", "1",
                                                    catchem::DiagType::FIELD_2D, {state->column_count(), 1});
    }
    void run(std::shared_ptr<catchem::StateManager> state) override {
        auto* concentration = state->chemistry().conc->host_write();
        concentration[0] += state->write_field<3>("T")[0] * 0.001;
        static_cast<double*>(state->diagnostic_manager()->get_host_write_pointer("extension_value"))[0] =
            concentration[0];
    }
    void finalize() override { finalized_ = true; }

private:
    bool& finalized_;
};

int main() {
    bool finalized = false;
    auto& registry = catchem::ProcessRegistry::get_instance();
    registry.register_process(
        "test-extension", [&] { return std::make_shared<ExtensionProcess>(finalized); },
        [&] { return ExtensionProcess(finalized).get_contract(); },
        [](const catchem::ProcessConfig& config) {
            if (config.scheme != "custom")
                throw std::invalid_argument("custom scheme required");
        });
    catchem::ProcessConfig settings;
    settings.scheme = "custom";
    registry.validate_settings("test-extension", settings);

    catchem::Core core(1, 1, 1);
    auto state = core.get_state_manager();
    auto mechanism = std::make_shared<catchem::MechanismDefinition>();
    catchem::SpeciesMetadata species;
    species.short_name = "unfamiliar-marker";
    species.roles = {"extension.marker"};
    mechanism->species.push_back(species);
    mechanism->rebuild_index();
    state->chemistry().mechanism = mechanism;
    std::vector<double> temperature{300.0}, concentration{1.0};
    state->bind_met_field_3d("T", temperature.data());
    state->bind_unified_chemistry(concentration.data());
    auto process = registry.create("test-extension");
    process->init(state);
    core.add_process(process);
    core.run_timestep(1.0);
    assert(concentration[0] == 1.3);
    core.shutdown();
    assert(finalized);
    return 0;
}
