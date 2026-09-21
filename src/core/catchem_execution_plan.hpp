#pragma once

#include "catchem_process_interface.hpp"
#include "catchem_validation.hpp"
#include <memory>
#include <vector>

namespace catchem {

    class ExecutionPlan {
    public:
        void compile(const std::vector<std::shared_ptr<ProcessInterface>>& processes,
                     const MechanismDefinition* mechanism);
        const ValidationIssueReport& validation() const noexcept { return validation_; }
        const ProcessContract& contract(std::size_t index) const { return contracts_.at(index); }
        std::size_t size() const noexcept { return contracts_.size(); }
        void prepare(std::size_t index, StateManager& state) const;
        void complete(std::size_t index, StateManager& state) const;

    private:
        std::vector<ProcessContract> contracts_;
        ValidationIssueReport validation_;
    };

} // namespace catchem
