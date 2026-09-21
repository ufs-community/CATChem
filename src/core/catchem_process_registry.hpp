#pragma once
#include "catchem_process_interface.hpp"
#include <functional>
#include <memory>
#include <mutex>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace catchem {

    using ProcessCreator = std::function<std::shared_ptr<ProcessInterface>()>;
    using ProcessContractFactory = std::function<ProcessContract()>;
    using ProcessSettingsValidator = std::function<void(const ProcessConfig&)>;

    struct ProcessRegistryEntry {
        ProcessCreator creator;
        ProcessContractFactory contract_factory;
        ProcessSettingsValidator settings_validator;
    };

    class ProcessRegistry {
    private:
        std::unordered_map<std::string, ProcessRegistryEntry> creators;
        mutable std::mutex mutex_;

        // Singleton private constructor
        ProcessRegistry() = default;

    public:
        static ProcessRegistry& get_instance() {
            static ProcessRegistry instance;
            return instance;
        }

        void register_process(const std::string& name, ProcessCreator creator,
                              ProcessContractFactory contract_factory = {},
                              ProcessSettingsValidator settings_validator = {}) {
            std::lock_guard<std::mutex> lock(mutex_);
            creators[name] = {std::move(creator), std::move(contract_factory), std::move(settings_validator)};
        }

        bool has_process(const std::string& name) const {
            std::lock_guard<std::mutex> lock(mutex_);
            return creators.find(name) != creators.end();
        }

        std::shared_ptr<ProcessInterface> create(const std::string& name) {
            ProcessCreator creator;
            {
                std::lock_guard<std::mutex> lock(mutex_);
                const auto found = creators.find(name);
                if (found == creators.end())
                    throw std::invalid_argument("Process not registered in C++: " + name);
                creator = found->second.creator;
            }
            return creator();
        }

        void validate_settings(const std::string& name, const ProcessConfig& settings) const {
            ProcessSettingsValidator validator;
            {
                std::lock_guard<std::mutex> lock(mutex_);
                const auto found = creators.find(name);
                if (found == creators.end())
                    throw std::invalid_argument("Process not registered in C++: " + name);
                validator = found->second.settings_validator;
            }
            if (validator)
                validator(settings);
        }

        void clear() {
            std::lock_guard<std::mutex> lock(mutex_);
            creators.clear();
        }
    };

    /**
     * @brief Build a settings validator that rejects unknown scheme options.
     *
     * Each process registers this alongside its creator so that every nested
     * `processes/<name>/<scheme>/<key>` entry in the runtime YAML must appear
     * in `accepted_paths`.  A misspelled or removed option therefore fails at
     * initialization with the accepted schema listed, instead of silently
     * leaving the scheme on its compiled default.  Paths use the same
     * slash-separated form as ProcessConfig::get_* (for example
     * "fengsha/alpha").
     */
    inline ProcessSettingsValidator make_settings_validator(std::string process_name,
                                                            std::vector<std::string> accepted_paths) {
        std::set<std::string> accepted(accepted_paths.begin(), accepted_paths.end());
        return [process_name = std::move(process_name), accepted = std::move(accepted)](const ProcessConfig& settings) {
            std::vector<std::string> unknown;
            for (const auto& path : settings.option_paths()) {
                if (!accepted.count(path))
                    unknown.push_back(path);
            }
            if (unknown.empty())
                return;
            std::ostringstream message;
            message << "Unknown option(s) for process '" << process_name << "':";
            for (const auto& key : unknown)
                message << " " << key << ",";
            message.seekp(-1, std::ios_base::end); // drop trailing comma
            message << ". Accepted options:";
            for (const auto& key : accepted)
                message << " " << key << ",";
            message.seekp(-1, std::ios_base::end);
            message << '.';
            throw std::invalid_argument(message.str());
        };
    }

} // namespace catchem
