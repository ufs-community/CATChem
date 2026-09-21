#include "catchem_config_manager.hpp"
#include "catchem_logger.hpp"
#include <algorithm>
#include <cctype>
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <stdexcept>

namespace catchem {

    // Default ON so standalone runs and tests keep echoing the effective
    // configuration; the NUOPC driver turns this off on non-root PETs.
    bool ConfigManager::echo_config_to_stdout = true;

    bool ValidationReport::has_errors() const {
        return std::any_of(issues.begin(), issues.end(),
                           [](const ValidationIssue& issue) { return issue.severity == ValidationSeverity::Error; });
    }

    std::string ValidationReport::format() const {
        std::ostringstream output;
        for (const auto& issue : issues) {
            output << (issue.severity == ValidationSeverity::Error ? "error" : "warning") << " [" << issue.category
                   << "] " << issue.path << ": " << issue.message;
            if (!issue.suggested_correction.empty())
                output << " (" << issue.suggested_correction << ")";
            output << '\n';
        }
        return output.str();
    }

    namespace {

        YAML::Node safe_get_map_child(const YAML::Node& node, std::string_view key) {
            if (!node.IsDefined() || !node.IsMap()) {
                return YAML::Node();
            }
            try {
                std::string k(key);
                YAML::Node child = node[k];
                if (child.IsDefined() && !child.IsNull()) {
                    return child;
                }
            } catch (...) {
            }
            return YAML::Node();
        }

        YAML::Node resolve_yaml_path(const YAML::Node& root, std::string_view path) {
            if (!root.IsDefined() || root.IsNull() || path.empty()) {
                return YAML::Node();
            }
            std::stringstream ss{std::string(path)};
            std::string segment;
            std::vector<YAML::Node> chain;
            chain.push_back(root);
            while (std::getline(ss, segment, '/')) {
                if (segment.empty())
                    continue;
                const YAML::Node& parent = chain.back();
                if (!parent.IsDefined() || parent.IsNull() || !parent.IsMap()) {
                    return YAML::Node();
                }
                YAML::Node child = safe_get_map_child(parent, segment);
                if (!child.IsDefined() || child.IsNull()) {
                    return YAML::Node();
                }
                chain.push_back(child);
            }
            return chain.back();
        }

        std::string resolve_relative_config_path(const std::string& config_file, const std::string& path) {
            if (path.empty() || path.front() == '/' || path == "none" || path == "NONE" || path == "null" ||
                path == "NULL")
                return path;
            if (config_file.empty())
                return path;
            auto slash = config_file.find_last_of('/');
            if (slash == std::string::npos)
                return path;
            return config_file.substr(0, slash + 1) + path;
        }

        std::string extract_file_setting(const YAML::Node& node) {
            if (!node)
                return "";
            if (node.IsScalar()) {
                try {
                    std::string path = node.as<std::string>();
                    if (!path.empty() && path != "none" && path != "NONE" && path != "null" && path != "NULL" &&
                        path != "~" && path != "true" && path != "false" && path != "TRUE" && path != "FALSE" &&
                        path != "1" && path != "0") {
                        return path;
                    }
                } catch (...) {
                }
                return "";
            }
            if (!node.IsMap())
                return "";
            for (const char* key : {"source_file", "input_file", "filename", "file", "path"}) {
                YAML::Node val = safe_get_map_child(node, key);
                if (val) {
                    std::string res = extract_file_setting(val);
                    if (!res.empty())
                        return res;
                }
            }
            return "";
        }

        template <typename T> T value_or(const YAML::Node& node, const T& default_value) {
            return node ? node.as<T>() : default_value;
        }

        std::vector<std::string> string_vector_or_empty(const YAML::Node& node) {
            std::vector<std::string> values;
            if (!node || !node.IsSequence()) {
                return values;
            }
            for (const auto& item : node) {
                values.push_back(item.as<std::string>());
            }
            return values;
        }

        std::vector<double> double_vector_or_empty(const YAML::Node& node) {
            std::vector<double> values;
            if (!node || !node.IsSequence()) {
                return values;
            }
            for (const auto& item : node) {
                values.push_back(item.as<double>());
            }
            return values;
        }

        void parse_processes(const YAML::Node& process_node, ConfigData& data) {
            if (!process_node || !process_node.IsMap()) {
                return;
            }
            for (const auto& process : process_node) {
                const std::string name = process.first.as<std::string>();
                const YAML::Node node = process.second;
                ProcessConfig config;
                config.activate = value_or<bool>(node["activate"], false);
                config.diagnostics = value_or<bool>(node["diagnostics"], false);
                if (node["scheme"]) {
                    config.scheme = node["scheme"].as<std::string>();
                }
                config.diag_species = string_vector_or_empty(node["diag_species"]);
                config.set_settings_node(YAML::Clone(node));
                data.processes[name] = config;
            }
        }

    } // namespace

    bool ProcessConfig::get_bool(std::string_view key, bool default_val) const {
        YAML::Node val = resolve_yaml_path(settings_node, key);
        if (val.IsDefined() && !val.IsNull()) {
            try {
                return val.as<bool>();
            } catch (...) {
                if (val.IsScalar()) {
                    std::string s = val.Scalar();
                    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return std::tolower(c); });
                    if (s == "true" || s == "1" || s == "yes" || s == "y" || s == "on")
                        return true;
                    if (s == "false" || s == "0" || s == "no" || s == "n" || s == "off")
                        return false;
                }
            }
        }
        return default_val;
    }

    double ProcessConfig::get_double(std::string_view key, double default_val) const {
        YAML::Node val = resolve_yaml_path(settings_node, key);
        if (val.IsDefined() && !val.IsNull()) {
            try {
                return val.as<double>();
            } catch (...) {
                if (val.IsScalar()) {
                    try {
                        return std::stod(val.Scalar());
                    } catch (...) {
                    }
                }
            }
        }
        return default_val;
    }

    int ProcessConfig::get_int(std::string_view key, int default_val) const {
        YAML::Node val = resolve_yaml_path(settings_node, key);
        if (val.IsDefined() && !val.IsNull()) {
            try {
                return val.as<int>();
            } catch (...) {
                if (val.IsScalar()) {
                    try {
                        return std::stoi(val.Scalar());
                    } catch (...) {
                    }
                }
            }
        }
        return default_val;
    }

    std::string ProcessConfig::get_string(std::string_view key, std::string_view default_val) const {
        YAML::Node val = resolve_yaml_path(settings_node, key);
        if (val.IsDefined() && !val.IsNull()) {
            try {
                std::string res = val.as<std::string>();
                if (res != "null" && res != "NULL" && res != "~") {
                    return res;
                }
            } catch (...) {
                if (val.IsScalar()) {
                    return val.Scalar();
                }
            }
        }
        return std::string(default_val);
    }

    std::vector<double> ProcessConfig::get_vector(std::string_view key) const {
        std::vector<double> values;
        YAML::Node val = resolve_yaml_path(settings_node, key);
        if (!val.IsDefined() || !val.IsSequence())
            return values;
        for (const auto& item : val) {
            try {
                values.push_back(item.as<double>());
            } catch (...) {
                throw std::invalid_argument("Process option '" + std::string(key) + "' must be a sequence of numbers.");
            }
        }
        return values;
    }

    std::vector<std::string> ProcessConfig::option_paths() const {
        // Keys the process layer consumes directly; they are not scheme
        // tuning options and must not appear in the validator schema.
        static const std::set<std::string> framework_keys = {"activate",   "diagnostics", "scheme",
                                                             "gas_scheme", "aero_scheme", "diag_species"};
        std::vector<std::string> paths;
        if (!settings_node.IsDefined() || !settings_node.IsMap())
            return paths;
        for (const auto& entry : settings_node) {
            const std::string block = entry.first.as<std::string>();
            if (framework_keys.count(block))
                continue;
            // Only nested scheme blocks are option containers.  Legacy flat
            // scalars (scheme_opt, dust_alpha, weibull, ...) are ignored by
            // the validator: they belong to the retired Fortran config path
            // and appear in several reference configurations.
            if (!entry.second.IsMap())
                continue;
            for (const auto& option : entry.second) {
                paths.push_back(block + "/" + option.first.as<std::string>());
            }
        }
        return paths;
    }

    bool ConfigManager::get_bool(std::string_view path, bool default_val) const {
        if (path.rfind("processes/", 0) == 0 || path.rfind("process/", 0) == 0) {
            std::size_t first_slash = path.find('/');
            std::size_t second_slash = path.find('/', first_slash + 1);
            if (second_slash != std::string_view::npos) {
                std::string proc_name(path.substr(first_slash + 1, second_slash - (first_slash + 1)));
                std::string_view key = path.substr(second_slash + 1);
                auto it = data.processes.find(proc_name);
                if (it != data.processes.end()) {
                    if (key == "activate")
                        return it->second.activate;
                    if (key == "diagnostics")
                        return it->second.diagnostics;
                    return it->second.get_bool(key, default_val);
                }
            }
        }
        YAML::Node node = resolve_yaml_path(root_node, path);
        if (node.IsDefined() && !node.IsNull()) {
            try {
                return node.as<bool>();
            } catch (const std::exception& e) {
                try {
                    std::string s = node.as<std::string>();
                    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return std::tolower(c); });
                    if (s == "true" || s == "1" || s == "yes" || s == "y" || s == "on")
                        return true;
                    if (s == "false" || s == "0" || s == "no" || s == "n" || s == "off")
                        return false;
                } catch (...) {
                }
            }
        }
        return default_val;
    }

    double ConfigManager::get_double(std::string_view path, double default_val) const {
        if (path.rfind("processes/", 0) == 0 || path.rfind("process/", 0) == 0) {
            std::size_t first_slash = path.find('/');
            std::size_t second_slash = path.find('/', first_slash + 1);
            if (second_slash != std::string_view::npos) {
                std::string proc_name(path.substr(first_slash + 1, second_slash - (first_slash + 1)));
                std::string_view key = path.substr(second_slash + 1);
                auto it = data.processes.find(proc_name);
                if (it != data.processes.end()) {
                    double val = it->second.get_double(key, default_val);
                    if (val != default_val)
                        return val;
                }
            }
        }
        YAML::Node node = resolve_yaml_path(root_node, path);
        if (node.IsDefined() && !node.IsNull()) {
            try {
                return node.as<double>();
            } catch (...) {
                if (node.IsScalar()) {
                    try {
                        return std::stod(node.Scalar());
                    } catch (...) {
                    }
                }
            }
        }
        return default_val;
    }

    int ConfigManager::get_int(std::string_view path, int default_val) const {
        if (path.rfind("processes/", 0) == 0 || path.rfind("process/", 0) == 0) {
            std::size_t first_slash = path.find('/');
            std::size_t second_slash = path.find('/', first_slash + 1);
            if (second_slash != std::string_view::npos) {
                std::string proc_name(path.substr(first_slash + 1, second_slash - (first_slash + 1)));
                std::string_view key = path.substr(second_slash + 1);
                auto it = data.processes.find(proc_name);
                if (it != data.processes.end()) {
                    int val = it->second.get_int(key, default_val);
                    if (val != default_val)
                        return val;
                }
            }
        }
        YAML::Node node = resolve_yaml_path(root_node, path);
        if (node.IsDefined() && !node.IsNull()) {
            try {
                return node.as<int>();
            } catch (...) {
                if (node.IsScalar()) {
                    try {
                        return std::stoi(node.Scalar());
                    } catch (...) {
                    }
                }
            }
        }
        return default_val;
    }

    std::string ConfigManager::get_string(std::string_view path, std::string_view default_val) const {
        if (path.rfind("processes/", 0) == 0 || path.rfind("process/", 0) == 0) {
            std::size_t first_slash = path.find('/');
            std::size_t second_slash = path.find('/', first_slash + 1);
            if (second_slash != std::string_view::npos) {
                std::string proc_name(path.substr(first_slash + 1, second_slash - (first_slash + 1)));
                std::string_view key = path.substr(second_slash + 1);
                auto it = data.processes.find(proc_name);
                if (it != data.processes.end()) {
                    if (key == "scheme") {
                        if (!it->second.scheme.empty())
                            return it->second.scheme;
                    }
                    std::string setting = it->second.get_string(key, "");
                    if (!setting.empty()) {
                        if (key.find("source_file") != std::string_view::npos ||
                            key.find("filename") != std::string_view::npos) {
                            return resolve_relative_config_path(config_file_path, setting);
                        }
                        return setting;
                    }
                }
            }
        }
        YAML::Node node = resolve_yaml_path(root_node, path);
        if (node.IsDefined() && !node.IsNull()) {
            try {
                std::string str_val = node.as<std::string>();
                if (str_val != "null" && str_val != "NULL" && str_val != "Null" && str_val != "~") {
                    if (path.find("source_file") != std::string_view::npos ||
                        path.find("filename") != std::string_view::npos) {
                        return resolve_relative_config_path(config_file_path, str_val);
                    }
                    return str_val;
                }
            } catch (...) {
            }
        }
        return std::string(default_val);
    }

    std::vector<std::string> ConfigManager::get_string_list(std::string_view path) const {
        std::vector<std::string> results;
        if (path.rfind("processes/", 0) == 0 || path.rfind("process/", 0) == 0) {
            std::size_t first_slash = path.find('/');
            std::size_t second_slash = path.find('/', first_slash + 1);
            if (second_slash != std::string_view::npos) {
                std::string proc_name(path.substr(first_slash + 1, second_slash - (first_slash + 1)));
                std::string_view key = path.substr(second_slash + 1);
                auto it = data.processes.find(proc_name);
                if (it != data.processes.end()) {
                    YAML::Node val = resolve_yaml_path(it->second.settings_node, key);
                    if (val && val.IsSequence()) {
                        for (const auto& item : val) {
                            if (item.IsScalar()) {
                                results.push_back(item.Scalar());
                            } else {
                                try {
                                    results.push_back(item.as<std::string>());
                                } catch (...) {
                                }
                            }
                        }
                        return results;
                    }
                }
            }
        }
        YAML::Node node = resolve_yaml_path(root_node, path);
        if (node && node.IsSequence()) {
            for (const auto& item : node) {
                if (item.IsScalar()) {
                    results.push_back(item.Scalar());
                } else {
                    try {
                        results.push_back(item.as<std::string>());
                    } catch (...) {
                    }
                }
            }
        }
        return results;
    }

    bool ConfigManager::is_process_active(std::string_view process_name) const {
        std::string p_path = "processes/" + std::string(process_name) + "/activate";
        if (get_bool(p_path, false))
            return true;
        std::string p_singular = "process/" + std::string(process_name) + "/activate";
        if (get_bool(p_singular, false))
            return true;
        auto it = data.processes.find(std::string(process_name));
        if (it != data.processes.end()) {
            return it->second.activate;
        }
        return false;
    }

    bool ConfigManager::is_category_active(std::string_view category_name) const {
        if (data.emission_mappings.find(std::string(category_name)) == data.emission_mappings.end()) {
            return false;
        }
        std::string path = "processes/extemis/" + std::string(category_name) + "/activate";
        if (get_bool(path, true))
            return true;
        std::string path_sing = "process/extemis/" + std::string(category_name) + "/activate";
        return get_bool(path_sing, true);
    }

    std::string ConfigManager::find_process_file_setting(std::string_view process_name) const {
        const auto direct = data.processes.find(std::string(process_name));
        if (direct != data.processes.end()) {
            std::string path = extract_file_setting(direct->second.settings_node);
            if (!path.empty()) {
                return resolve_relative_config_path(config_file_path, path);
            }
        }

        const auto dust = data.processes.find("dust");
        if (dust != data.processes.end()) {
            std::string path = extract_file_setting(safe_get_map_child(dust->second.settings_node, process_name));
            if (path.empty()) {
                path = extract_file_setting(dust->second.settings_node);
            }
            if (!path.empty()) {
                return resolve_relative_config_path(config_file_path, path);
            }
        }

        const auto extemis = data.processes.find("extemis");
        if (extemis != data.processes.end()) {
            std::string path = extract_file_setting(safe_get_map_child(extemis->second.settings_node, process_name));
            if (path.empty()) {
                path = extract_file_setting(safe_get_map_child(extemis->second.settings_node, "dust"));
            }
            if (!path.empty()) {
                return resolve_relative_config_path(config_file_path, path);
            }
        }

        return "";
    }

    void ConfigManager::load_from_file(const std::string& filename) {
        config_file_path = filename;
        validation_report.issues.clear();
        try {
            root_node = YAML::LoadFile(filename);
            // Echo the effective configuration to stdout so the run log (e.g.
            // from the NUOPC cap) records exactly what the core parsed.  Under
            // NUOPC the driver keeps this on the root PET only (every PET
            // parses the same file); see echo_config_to_stdout.
            if (echo_config_to_stdout) {
                std::cout << "--- CATChem configuration: " << filename << " ---" << std::endl;
                std::cout << root_node << std::endl;
                std::cout << "--- End CATChem configuration ---" << std::endl;
            }
            const YAML::Node& config = root_node;
            data.active_processes.clear();
            data.processes.clear();
            if (config["simulation"]) {
                auto sim = config["simulation"];
                data.simulation.name = value_or<std::string>(sim["name"], "");
                data.simulation.start_date = value_or<std::string>(sim["start_date"], "");
                data.simulation.end_date = value_or<std::string>(sim["end_date"], "");
                if (sim["species_filename"]) {
                    data.species_filename = sim["species_filename"].as<std::string>();
                    data.simulation.species_filename = data.species_filename;
                }
                if (sim["emission_filename"]) {
                    data.simulation.emission_filename = sim["emission_filename"].as<std::string>();
                }
                if (sim["verbose"]) {
                    data.simulation.verbose_enabled = value_or<bool>(sim["verbose"]["activate"], false);
                    // Optional explicit logger threshold.  When present it
                    // always wins over the CATCHEM_LOG_LEVEL environment
                    // variable; when absent the previous setting is kept.
                    if (sim["verbose"]["log_level"]) {
                        data.simulation.log_level = sim["verbose"]["log_level"].as<std::string>();
                        const auto level = Logger::level_from_string(data.simulation.log_level);
                        if (!level)
                            throw std::invalid_argument(
                                "simulation/verbose/log_level must be debug, info, warn, or error");
                        Logger::set_level(*level);
                    }
                }
                if (sim["nx"]) {
                    data.runtime.nx = sim["nx"].as<int>();
                }
                if (sim["ny"]) {
                    data.runtime.ny = sim["ny"].as<int>();
                }
                if (sim["nz"]) {
                    data.runtime.nz = sim["nz"].as<int>();
                }
                if (sim["timestep"]) {
                    data.runtime.dt = sim["timestep"].as<double>();
                } else if (sim["dt"]) {
                    data.runtime.dt = sim["dt"].as<double>();
                }
                if (sim["nsteps"]) {
                    data.runtime.nsteps = sim["nsteps"].as<int>();
                }
            }
            if (config["mechanism"] && config["mechanism"].IsMap()) {
                const auto mechanism = config["mechanism"];
                data.mechanism_identity = value_or<std::string>(mechanism["identity"], "");
                data.mechanism_capabilities = string_vector_or_empty(mechanism["capabilities"]);
            }
            if (config["physical_validation"] && config["physical_validation"].IsMap()) {
                const auto policy = value_or<std::string>(config["physical_validation"]["policy"], "reject");
                if (policy == "reject")
                    data.physical_validation_policy = PhysicalValidationPolicy::Reject;
                else if (policy == "warn_and_clamp")
                    data.physical_validation_policy = PhysicalValidationPolicy::WarnAndClamp;
                else if (policy == "count_and_continue")
                    data.physical_validation_policy = PhysicalValidationPolicy::CountAndContinue;
                else
                    throw std::invalid_argument(
                        "physical_validation/policy must be reject, warn_and_clamp, or count_and_continue");
            }
            if (config["grid"]) {
                const YAML::Node grid = config["grid"];
                data.grid.number_of_levels = value_or<int>(grid["number_of_levels"], data.grid.number_of_levels);
                data.grid.number_of_soil_layers =
                    value_or<int>(grid["number_of_soil_layers"], data.grid.number_of_soil_layers);
                data.runtime.nz = data.grid.number_of_levels;
            }
            if (config["timesteps"]) {
                const YAML::Node timesteps = config["timesteps"];
                data.timesteps.transport_timestep_in_s =
                    value_or<int>(timesteps["transport_timestep_in_s"], data.timesteps.transport_timestep_in_s);
                data.timesteps.chemistry_timestep_in_s =
                    value_or<int>(timesteps["chemistry_timestep_in_s"], data.timesteps.chemistry_timestep_in_s);
            }
            if (config["diagnostics"]) {
                const YAML::Node diagnostics = config["diagnostics"];
                if (diagnostics["output"]) {
                    const YAML::Node output = diagnostics["output"];
                    data.diagnostics.output.enabled = value_or<bool>(output["enabled"], false);
                    data.diagnostics.output.directory = value_or<std::string>(output["directory"], "");
                    data.diagnostics.output.prefix = value_or<std::string>(output["prefix"], "");
                    data.diagnostics.output.frequency = value_or<int>(output["frequency"], 0);
                    data.diagnostics.output.format = value_or<std::string>(output["format"], "");
                    data.diagnostics.output.compress_lev = value_or<int>(output["compress_lev"], 0);
                    data.diagnostics.output.process_diagnostics = value_or<bool>(output["process_diagnostics"], false);
                    data.diagnostics.output.diag_list = string_vector_or_empty(output["diag_list"]);
                    if (output["attributes"] && output["attributes"].IsMap()) {
                        for (const auto& attr : output["attributes"]) {
                            if (attr.second)
                                data.diagnostics.output.attributes[attr.first.as<std::string>()] =
                                    attr.second.as<std::string>();
                        }
                    }
                }
                if (diagnostics["collection"]) {
                    const YAML::Node collection = diagnostics["collection"];
                    data.diagnostics.collection.enabled = value_or<bool>(collection["enabled"], false);
                    data.diagnostics.collection.buffer_size = value_or<int>(collection["buffer_size"], 0);
                }
            }
            if (config["mie"] && config["mie"].IsMap()) {
                const YAML::Node mie = config["mie"];
                data.mie = MieConfig{};
                data.mie.directory = value_or<std::string>(mie["directory"], "./");
                if (mie["files"] && mie["files"].IsMap()) {
                    // yaml-cpp preserves the document order of map children, so the
                    // load sequence below is deterministic across ranks.
                    for (const auto& entry : mie["files"]) {
                        const std::string type = entry.first.as<std::string>();
                        const std::string file = entry.second ? entry.second.as<std::string>() : "";
                        if (type.empty())
                            throw std::invalid_argument("mie.files contains an empty aerosol type key");
                        if (file.empty())
                            throw std::invalid_argument("mie.files entry '" + type + "' has no optics file name");
                        data.mie.files.emplace_back(type, file);
                    }
                }
            }
            parse_processes(config["processes"], data);
            parse_processes(config["process"], data);
            if (config["run_phases"]) {
                for (const auto& phase : config["run_phases"]) {
                    const YAML::Node processes = phase.second["processes"];
                    if (!processes || !processes.IsSequence()) {
                        continue;
                    }
                    for (const auto& process : processes) {
                        data.active_processes.push_back(process.as<std::string>());
                    }
                }
            }
            is_loaded = true;
        } catch (const std::exception& e) {
            std::cerr << "Error loading configuration file " << filename << ": " << e.what() << std::endl;
            is_loaded = false;
            throw;
        }
    }

    void ConfigManager::load_species_file(const std::string& filename) {
        const YAML::Node species_root = YAML::LoadFile(filename);
        data.species.clear();
        if (!species_root || !species_root.IsSequence()) {
            return;
        }

        auto get_bool = [](const YAML::Node& n, const std::string& key, bool def) {
            if (n["__" + key]) {
                try {
                    return n["__" + key].as<bool>();
                } catch (...) {
                }
            }
            if (n[key]) {
                try {
                    return n[key].as<bool>();
                } catch (...) {
                }
            }
            return def;
        };
        auto get_double = [](const YAML::Node& n, const std::string& key, double def) {
            if (n["__" + key]) {
                try {
                    return n["__" + key].as<double>();
                } catch (...) {
                }
            }
            if (n[key]) {
                try {
                    return n[key].as<double>();
                } catch (...) {
                }
            }
            return def;
        };
        auto get_str = [](const YAML::Node& n, const std::string& key, const std::string& def) {
            if (n["__" + key]) {
                try {
                    return n["__" + key].as<std::string>();
                } catch (...) {
                }
            }
            if (n[key]) {
                try {
                    return n[key].as<std::string>();
                } catch (...) {
                }
            }
            return def;
        };

        for (const auto& species_node : species_root) {
            SpeciesConfig species;
            species.name = value_or<std::string>(species_node["name"], "");
            species.long_name = get_str(species_node, "long_name", species.name);
            species.description = get_str(species_node, "description", "");
            species.aliases = string_vector_or_empty(species_node["aliases"]);
            species.roles = string_vector_or_empty(species_node["roles"]);

            species.is_gas = get_bool(species_node, "is_gas", false);
            species.is_aerosol = get_bool(species_node, "is_aerosol", false);
            species.is_tracer = get_bool(species_node, "is_tracer", false);
            species.is_advected = get_bool(species_node, "is_advected", true);
            species.is_drydep = get_bool(species_node, "is_drydep", false);
            species.is_wetdep = get_bool(species_node, "is_wetdep", false);
            species.is_photolysis = get_bool(species_node, "is_photolysis", false);
            species.is_gocart_aero = get_bool(species_node, "is_gocart_aero", false);
            species.is_dust = get_bool(species_node, "is_dust", false);
            species.is_seasalt = get_bool(species_node, "is_seasalt", false);
            // Hydrophilic aerosols take on water at elevated RH and settle at
            // their wet size; hydrophobic ones (e.g. fresh carbon, dust) do
            // not.  Defaults to true so an omitted attribute keeps the common
            // (soluble aerosol) behavior; set false for dust and the fresh
            // (phobic) carbon bins.
            species.is_hydrophilic = get_bool(species_node, "is_hydrophilic", true);

            if (species_node["molecular weight [kg mol-1]"]) {
                species.molecular_weight_kg_mol = species_node["molecular weight [kg mol-1]"].as<double>();
            } else if (species_node["molecular_weight_kg_mol"]) {
                species.molecular_weight_kg_mol = species_node["molecular_weight_kg_mol"].as<double>();
            } else {
                species.molecular_weight_kg_mol = get_double(species_node, "mw_g", 0.0) / 1000.0;
            }
            species.mw_g = species.molecular_weight_kg_mol * 1000.0;

            species.density = get_double(species_node, "density", 0.0);
            species.radius = get_double(species_node, "radius", 0.0);
            species.lower_radius = get_double(species_node, "lower_radius", 0.0);
            species.upper_radius = get_double(species_node, "upper_radius", 0.0);
            species.viscosity = get_double(species_node, "viscosity", 0.0);

            species.dd_f0 = get_double(species_node, "dd_f0", 0.0);
            species.dd_hstar = get_double(species_node, "dd_hstar", 0.0);
            species.dd_DvzAerSnow = get_double(species_node, "dd_DvzAerSnow", 0.0);
            species.dd_DvzMinVal_snow = get_double(species_node, "dd_DvzMinVal_snow", 0.0);
            species.dd_DvzMinVal_land = get_double(species_node, "dd_DvzMinVal_land", 0.0);

            species.henry_k0 = get_double(species_node, "henry_k0", 0.0);
            species.henry_cr = get_double(species_node, "henry_cr", 0.0);
            species.henry_pKa = get_double(species_node, "henry_pKa", 0.0);
            species.wd_retfactor = get_double(species_node, "wd_retfactor", 0.0);
            species.wd_LiqAndGas = get_bool(species_node, "wd_LiqAndGas", false);
            species.wd_convfacI2G = get_double(species_node, "wd_convfacI2G", 0.0);

            if (species_node["__wd_rainouteff"]) {
                species.wd_rainouteff = species_node["__wd_rainouteff"].as<std::vector<double>>();
            } else if (species_node["wd_rainouteff"]) {
                species.wd_rainouteff = species_node["wd_rainouteff"].as<std::vector<double>>();
            }
            species.wd_reevap_frac = get_double(species_node, "wd_reevap_frac", 0.5);

            species.t_chem_loss = get_double(species_node, "t_chem_loss", -1.0);
            species.BackgroundVV = get_double(species_node, "BackgroundVV", 1.0e-20);
            species.mie_name = get_str(species_node, "mie_name", "");

            data.species.push_back(species);
        }
    }

    void ConfigManager::load_emission_mapping_file(const std::string& filename) {
        const YAML::Node emission_root = YAML::LoadFile(filename);
        data.emission_mappings.clear();
        if (!emission_root || !emission_root.IsMap()) {
            return;
        }
        for (const auto& category_node : emission_root) {
            const std::string category_name = category_node.first.as<std::string>();
            EmissionCategoryMapping category;
            for (const auto& field_node : category_node.second) {
                const std::string field_name = field_node.first.as<std::string>();
                const YAML::Node node = field_node.second;
                EmissionFieldMapping field;
                field.long_name = value_or<std::string>(node["long_name"], "");
                field.units = value_or<std::string>(node["units"], "");
                field.scale = double_vector_or_empty(node["scale"]);
                field.map = string_vector_or_empty(node["map"]);
                category.fields[field_name] = field;
            }
            data.emission_mappings[category_name] = category;
        }
    }

    const ValidationReport& ConfigManager::validate(bool strict) {
        validation_report.issues.clear();
        auto add = [this](std::string category, std::string path, std::string message, std::string correction,
                          ValidationSeverity severity = ValidationSeverity::Error) {
            validation_report.issues.push_back(
                {severity, std::move(category), std::move(path), std::move(message), std::move(correction)});
        };

        if (!root_node || !root_node.IsMap())
            add("schema", "/", "configuration root must be a mapping", "provide a YAML mapping");
        if (data.runtime.nx <= 0)
            add("range", "simulation/nx", "must be positive", "set nx > 0");
        if (data.runtime.ny <= 0)
            add("range", "simulation/ny", "must be positive", "set ny > 0");
        if (data.runtime.nz <= 0)
            add("range", "grid/number_of_levels", "must be positive", "set levels > 0");
        if (data.runtime.dt <= 0.0 || data.runtime.dt > 86400.0)
            add("range", "simulation/timestep", "must satisfy 0 < timestep <= 86400", "set a valid timestep");
        if (data.runtime.nsteps <= 0)
            add("range", "simulation/nsteps", "must be positive", "set nsteps > 0");
        if (data.grid.number_of_soil_layers < 0)
            add("range", "grid/number_of_soil_layers", "must not be negative", "set soil layers >= 0");
        if (data.timesteps.transport_timestep_in_s < 0)
            add("range", "timesteps/transport_timestep_in_s", "must not be negative", "set a nonnegative value");
        if (data.timesteps.chemistry_timestep_in_s < 0)
            add("range", "timesteps/chemistry_timestep_in_s", "must not be negative", "set a nonnegative value");
        if (data.diagnostics.output.enabled) {
            if (data.diagnostics.output.frequency <= 0)
                add("diagnostics", "diagnostics/output/frequency", "must be positive when output is enabled",
                    "set frequency > 0");
            if (data.diagnostics.output.compress_lev < 0 || data.diagnostics.output.compress_lev > 9)
                add("diagnostics", "diagnostics/output/compress_lev", "must be between 0 and 9",
                    "choose a compression level from 0 through 9");
            if (!data.diagnostics.output.format.empty() && data.diagnostics.output.format != "netcdf")
                add("diagnostics", "diagnostics/output/format", "unsupported output format", "use netcdf");
        }
        if (data.diagnostics.collection.enabled && data.diagnostics.collection.buffer_size <= 0)
            add("diagnostics", "diagnostics/collection/buffer_size", "must be positive when collection is enabled",
                "set buffer_size > 0");

        if (data.species_filename.empty())
            add("mechanism", "simulation/species_filename", "a configured run requires a mechanism",
                "provide a species mechanism file");
        if (data.species.empty())
            add("mechanism", "species", "mechanism must contain at least one species",
                "provide a nonempty species list");
        std::set<std::string> species_names;
        for (std::size_t index = 0; index < data.species.size(); ++index) {
            std::string canonical = data.species[index].name;
            std::transform(canonical.begin(), canonical.end(), canonical.begin(),
                           [](unsigned char value) { return static_cast<char>(std::toupper(value)); });
            const std::string path = "species[" + std::to_string(index) + "]/name";
            if (canonical.empty())
                add("mechanism", path, "species name must not be empty", "set a canonical name");
            else if (!species_names.insert(canonical).second)
                add("mechanism", path, "duplicate species name after case normalization",
                    "remove or rename the duplicate");
        }

        const std::set<std::string> known_processes = {"carbchem",   "drydep",  "dust",     "extemis", "gaschem",
                                                       "photolysis", "seasalt", "settling", "so4chem", "wetdep"};
        for (const auto& process : data.active_processes) {
            if (known_processes.find(process) == known_processes.end())
                add("process", "run_phases/processes", "unknown process name: " + process,
                    "use a registered CATChem process");
        }
        for (const auto& [category_name, category] : data.emission_mappings) {
            for (const auto& [field_name, field] : category.fields) {
                const std::string base = "emissions/" + category_name + "/" + field_name;
                if (field.map.size() != field.scale.size())
                    add("emissions", base, "map and scale arrays have different lengths",
                        "provide one scale for each mapped species");
                for (const auto& target : field.map) {
                    std::string canonical = target;
                    std::transform(canonical.begin(), canonical.end(), canonical.begin(),
                                   [](unsigned char value) { return static_cast<char>(std::toupper(value)); });
                    // Emission mappings may also carry meteorological inputs.
                    // These use the runtime-resolved MET_ namespace and are
                    // intentionally not mechanism species.  Keep this rule
                    // namespace-based so configurations remain mechanism
                    // agnostic (and do not require a list of field names).
                    if (canonical.rfind("MET_", 0) == 0)
                        continue;
                    if (species_names.find(canonical) == species_names.end())
                        add("emissions", base + "/map", "target is absent from active mechanism: " + target,
                            "map to a mechanism species");
                }
            }
        }

        if (strict && root_node && root_node.IsMap()) {
            const std::set<std::string> allowed = {
                "simulation",  "mechanism", "physical_validation", "grid",      "timesteps",
                "diagnostics", "mie",       "run_phases",          "processes", "process"};
            for (const auto& entry : root_node) {
                const std::string key = entry.first.as<std::string>();
                if (allowed.find(key) == allowed.end())
                    add("schema", key, "unknown top-level configuration key", "remove it or use compatibility mode");
            }
        }
        return validation_report;
    }

    void ConfigManager::validate_or_throw(bool strict) {
        const auto& report = validate(strict);
        if (report.has_errors())
            throw std::invalid_argument("Invalid CATChem configuration:\n" + report.format());
    }

} // namespace catchem
