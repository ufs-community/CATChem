#pragma once

#include "catchem_process_interface.hpp"
#include "catchem_species_metadata.hpp"
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace catchem::test {

    inline MechanismDefinition synthetic_mechanism(std::vector<std::string> names, std::string identity = "synthetic") {
        MechanismDefinition mechanism;
        mechanism.identity = std::move(identity);
        mechanism.source = "test-fixture";
        for (auto& name : names) {
            SpeciesMetadata metadata;
            metadata.short_name = std::move(name);
            metadata.mw_g = 40.0;
            mechanism.species.push_back(std::move(metadata));
        }
        mechanism.rebuild_index();
        return mechanism;
    }

    class RecordingProcess : public ProcessInterface {
    public:
        RecordingProcess(std::string name, std::vector<std::string>& events, bool fail_init = false,
                         bool fail_run = false)
            : name_(std::move(name)), events_(events), fail_init_(fail_init), fail_run_(fail_run) {}
        std::string get_name() const override { return name_; }
        void init(std::shared_ptr<StateManager>) override {
            events_.push_back("init:" + name_);
            if (fail_init_)
                throw std::runtime_error("injected init failure");
        }
        void run(std::shared_ptr<StateManager>) override {
            events_.push_back("run:" + name_);
            if (fail_run_)
                throw std::runtime_error("injected run failure");
        }
        void finalize() override { events_.push_back("finalize:" + name_); }

    private:
        std::string name_;
        std::vector<std::string>& events_;
        bool fail_init_;
        bool fail_run_;
    };

} // namespace catchem::test
