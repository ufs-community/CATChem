#include "catchem_core_architecture_test_helpers.hpp"
#include "catchem_execution_plan.hpp"
#include "catchem_process_drydep.hpp"
#include "catchem_process_dust.hpp"
#include "catchem_process_seasalt.hpp"
#include "catchem_process_so4chem.hpp"
#include "catchem_process_wetdep.hpp"
#include <algorithm>
#include <cassert>

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
    return 0;
}
