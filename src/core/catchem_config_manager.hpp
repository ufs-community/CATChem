#pragma once
#include <string>
#include <string_view>
#include <yaml-cpp/yaml.h>

namespace catchem {

    struct RuntimeConfig {
        int nx = 1;
        int ny = 1;
        int nz = 1;
        double dt = 3600.0;
        int nsteps = 1;
    };

    struct ConfigData {
        RuntimeConfig runtime;
        // We can add FilePathConfig, etc. here later
    };

    class ConfigManager {
    public:
        ConfigData data;
        YAML::Node root_node;
        bool is_loaded = false;
        std::string config_file_path;

        ConfigManager() = default;
        void load_from_file(const std::string& filename);

        YAML::Node get_process_config(std::string_view process_name) const {
            if (is_loaded && root_node["process"]) {
                std::string key(process_name);
                if (root_node["process"][key]) {
                    return root_node["process"][key];
                }
            }
            return YAML::Node();
        }
    };

} // namespace catchem
