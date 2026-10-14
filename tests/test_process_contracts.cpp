#include "catchem_config_manager.hpp"
#include "catchem_core_architecture_test_helpers.hpp"
#include "catchem_execution_plan.hpp"
#include "catchem_process_drydep.hpp"
#include "catchem_process_dust.hpp"
#include "catchem_process_registry.hpp"
#include "catchem_process_seasalt.hpp"
#include "catchem_process_so4chem.hpp"
#include "catchem_process_wetdep.hpp"
#include <algorithm>
#include <cassert>
#include <iostream>
#include <string>
#include <yaml-cpp/yaml.h>

extern "C" {
void catchem_register_drydep_cpp();
void catchem_register_so4chem_cpp();
}

class ContractProcess : public catchem::test::RecordingProcess {
public:
    using RecordingProcess::RecordingProcess;
    catchem::ProcessContract get_contract() const override {
        return catchem::make_contract(get_name(), {{"TEMPERATURE",
                                                    "K",
                                                    {catchem::SemanticAxis::Column, catchem::SemanticAxis::Level,
                                                     catchem::SemanticAxis::Singleton},
                                                    catchem::PersistencePolicy::Timestep,
                                                    catchem::FieldRequirement::Required,
                                                    catchem::AccessIntent::Read,
                                                    catchem::ExecutionSpaceIntent::Host}});
    }
};

class ProducerProcess : public catchem::test::RecordingProcess {
public:
    using RecordingProcess::RecordingProcess;
    catchem::ProcessContract get_contract() const override {
        auto output =
            catchem::host_field_3d("DERIVED", "1", catchem::FieldRequirement::Required, catchem::AccessIntent::Write);
        output.produced = true;
        return catchem::make_contract(get_name(), {output});
    }
};

class ConsumerProcess : public catchem::test::RecordingProcess {
public:
    using RecordingProcess::RecordingProcess;
    catchem::ProcessContract get_contract() const override {
        return catchem::make_contract(get_name(), {catchem::host_field_3d("DERIVED", "1")});
    }
};

int main() {
    std::vector<std::string> events;
    std::vector<std::shared_ptr<catchem::ProcessInterface>> processes;
    processes.push_back(std::make_shared<ContractProcess>("contract", events));
    catchem::ExecutionPlan plan;
    plan.compile(processes, nullptr);
    assert(!plan.validation().has_errors());
    assert(plan.contract(0).fields.front().canonical_name == "TEMPERATURE");

    const std::vector<catchem::ProcessContract> surface_builtins{
        catchem::DustProcess().get_contract(), catchem::SeaSaltProcess().get_contract(),
        catchem::DryDepProcess().get_contract(), catchem::SO4chemProcess().get_contract()};
    for (const auto& contract : surface_builtins) {
        assert(contract.structurally_valid());
        bool has_surface_input = false;
        for (const auto& field : contract.fields)
            has_surface_input = has_surface_input || field.axes.size() == 2;
        assert(has_surface_input);
    }
    const auto dust_contract = catchem::DustProcess().get_contract();
    const auto soil_field =
        std::find_if(dust_contract.fields.begin(), dust_contract.fields.end(),
                     [](const catchem::FieldAccessContract& field) { return field.canonical_name == "SOILM"; });
    assert(soil_field != dust_contract.fields.end());
    assert(soil_field->units == "m3/m3");
    const std::vector<catchem::SemanticAxis> soil_axes{catchem::SemanticAxis::Column, catchem::SemanticAxis::SoilLayer,
                                                       catchem::SemanticAxis::Singleton};
    assert(soil_field->axes == soil_axes);
    assert(catchem::WetDepProcess().get_contract().structurally_valid());

    std::vector<std::shared_ptr<catchem::ProcessInterface>> reversed;
    reversed.push_back(std::make_shared<ConsumerProcess>("consumer", events));
    reversed.push_back(std::make_shared<ProducerProcess>("producer", events));
    plan.compile(reversed, nullptr);
    assert(plan.validation().has_errors());
    assert(plan.validation().format().find("dependency-order") != std::string::npos);

    std::vector<std::shared_ptr<catchem::ProcessInterface>> ordered;
    ordered.push_back(std::make_shared<ProducerProcess>("producer", events));
    ordered.push_back(std::make_shared<ConsumerProcess>("consumer", events));
    plan.compile(ordered, nullptr);
    assert(!plan.validation().has_errors());

    // --- US2 settings-validator allowlists (FR-021, contract K-3) ------------
    // The three GOCART-faithful routing keys must be accepted by the
    // registered per-process validators, and a misspelled sibling in the same
    // scheme block must still be rejected, so a typo cannot silently leave a
    // scheme on its compiled default.
    catchem_register_drydep_cpp();
    catchem_register_so4chem_cpp();
    auto& registry = catchem::ProcessRegistry::get_instance();
    auto accept = [&](const char* process, const YAML::Node& settings) {
        catchem::ProcessConfig config;
        config.activate = true;
        config.set_settings_node(settings);
        try {
            registry.validate_settings(process, config);
        } catch (const std::exception& error) {
            std::cerr << "FAIL: validator rejected legal option for " << process << ": " << error.what() << '\n';
            return false;
        }
        return true;
    };
    auto reject = [&](const char* process, const YAML::Node& settings, const std::string& expected_key) {
        catchem::ProcessConfig config;
        config.activate = true;
        config.set_settings_node(settings);
        try {
            registry.validate_settings(process, config);
        } catch (const std::exception& error) {
            return std::string(error.what()).find(expected_key) != std::string::npos;
        }
        return false;
    };
    YAML::Node drydep_ok;
    drydep_ok["wesely"]["skip_so2"] = true;
    drydep_ok["gocart"]["skip_sulfate_aero"] = true;
    assert(accept("drydep", drydep_ok));
    YAML::Node drydep_bad;
    drydep_bad["wesely"]["skip_so3"] = true;
    assert(reject("drydep", drydep_bad, "wesely/skip_so3"));
    YAML::Node so4chem_ok;
    so4chem_ok["gocart"]["do_drydep"] = true;
    assert(accept("so4chem", so4chem_ok));
    YAML::Node so4chem_bad;
    so4chem_bad["gocart"]["do_wetdep"] = true;
    assert(reject("so4chem", so4chem_bad, "gocart/do_wetdep"));
    return 0;
}
