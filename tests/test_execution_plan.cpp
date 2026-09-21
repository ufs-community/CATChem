#include "catchem_core_architecture_test_helpers.hpp"
#include "catchem_execution_plan.hpp"
#include <cassert>

class OptionalProcess : public catchem::test::RecordingProcess {
public:
    using RecordingProcess::RecordingProcess;
    catchem::ProcessContract get_contract() const override {
        return catchem::make_contract(get_name(), {{"UNRELATED",
                                                    "1",
                                                    {catchem::SemanticAxis::Column},
                                                    catchem::PersistencePolicy::Timestep,
                                                    catchem::FieldRequirement::Optional,
                                                    catchem::AccessIntent::Read,
                                                    catchem::ExecutionSpaceIntent::Host}});
    }
};

class InactiveRequiredProcess : public catchem::test::RecordingProcess {
public:
    using RecordingProcess::RecordingProcess;
    catchem::ProcessContract get_contract() const override {
        return catchem::make_contract(get_name(), {catchem::host_field_3d("NOT_BOUND", "1")});
    }
};

int main() {
    catchem::StateManager state(2, 3, 1);
    std::vector<std::string> events;
    std::vector<std::shared_ptr<catchem::ProcessInterface>> processes{
        std::make_shared<OptionalProcess>("optional", events)};
    catchem::ExecutionPlan plan;
    plan.compile(processes, nullptr);
    // The registry may contain many factories, but only the processes selected
    // by this runtime instance's YAML are supplied to schedule compilation.
    auto inactive = std::make_shared<InactiveRequiredProcess>("inactive", events);
    (void)inactive;
    assert(plan.size() == 1);
    plan.prepare(0, state); // absent optional field must not fail
    return 0;
}
