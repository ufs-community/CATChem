#include "catchem_execution_plan.hpp"
#include <stdexcept>
#include <unordered_map>

namespace catchem {

    void ExecutionPlan::compile(const std::vector<std::shared_ptr<ProcessInterface>>& processes,
                                const MechanismDefinition* mechanism) {
        contracts_.clear();
        validation_ = {};
        std::unordered_map<std::string, std::pair<std::size_t, FieldAccessContract>> produced;
        for (std::size_t process_index = 0; process_index < processes.size(); ++process_index) {
            const auto& process = processes[process_index];
            auto contract = process->get_contract();
            if (contract.process_name.empty())
                contract.process_name = process->get_name();
            if (!contract.structurally_valid())
                validation_.add({IssueSeverity::Error, "process-contract", contract.process_name,
                                 "process contract is structurally incomplete", "declare identity, units, and axes"});
            if (mechanism) {
                for (const auto& requirement : contract.mechanism_requirements) {
                    if (!requirement.role.empty() && requirement.required && !mechanism->has_role(requirement.role))
                        validation_.add({IssueSeverity::Error, "mechanism-role", contract.process_name,
                                         "missing required role " + requirement.role, "map the role in the mechanism"});
                    if (!requirement.capability.empty() && requirement.required &&
                        !mechanism->has_capability(requirement.capability))
                        validation_.add({IssueSeverity::Error, "mechanism-capability", contract.process_name,
                                         "missing required capability " + requirement.capability,
                                         "declare capability in the mechanism"});
                }
            }
            for (const auto& field : contract.fields) {
                const auto key = canonicalize_field_identity(field.canonical_name);
                if (field.produced) {
                    const auto prior = produced.find(key);
                    if (prior != produced.end() &&
                        (prior->second.second.units != field.units || prior->second.second.axes != field.axes))
                        validation_.add({IssueSeverity::Error, "incompatible-producers", contract.process_name,
                                         "output " + key + " conflicts with an earlier producer",
                                         "use one units/rank/axis contract for the shared output"});
                    else
                        produced.emplace(key, std::make_pair(process_index, field));
                }
            }
            contracts_.push_back(std::move(contract));
        }

        for (std::size_t consumer_index = 0; consumer_index < contracts_.size(); ++consumer_index) {
            for (const auto& field : contracts_[consumer_index].fields) {
                if (!field.reads() || field.requirement != FieldRequirement::Required)
                    continue;
                const auto producer = produced.find(canonicalize_field_identity(field.canonical_name));
                if (producer != produced.end() && producer->second.first > consumer_index)
                    validation_.add({IssueSeverity::Error, "dependency-order", contracts_[consumer_index].process_name,
                                     "required field " + canonicalize_field_identity(field.canonical_name) +
                                         " is produced later in the schedule",
                                     "move its producer before this consumer in active_processes"});
            }
        }
    }

    void ExecutionPlan::prepare(std::size_t index, StateManager& state) const {
        const auto& process = contracts_.at(index);
        for (const auto& access : process.fields) {
            const bool found = state.prepare_field_access(access);
            if (!found && access.requirement == FieldRequirement::Required)
                throw std::runtime_error("Process " + process.process_name + " missing required field " +
                                         access.canonical_name);
        }
    }

    void ExecutionPlan::complete(std::size_t index, StateManager& state) const {
        for (const auto& access : contracts_.at(index).fields)
            state.complete_field_access(access);
    }

} // namespace catchem
