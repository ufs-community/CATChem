#include "catchem_api.hpp"
#include <yaml-cpp/yaml.h>

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

extern "C" {
void catchem_register_seasalt_cpp();
void catchem_register_dust_cpp();
void catchem_register_carbchem_cpp();
void catchem_register_settling_cpp();
void catchem_register_drydep_cpp();
void catchem_register_so4chem_cpp();
void catchem_register_wetdep_cpp();
}

void register_candidate_processes() {
    // Static process libraries are not pulled into every linker/toolchain
    // solely for their registration symbols.  Register them explicitly.
    catchem_register_seasalt_cpp();
    catchem_register_dust_cpp();
    catchem_register_carbchem_cpp();
    catchem_register_settling_cpp();
    catchem_register_drydep_cpp();
    catchem_register_so4chem_cpp();
    catchem_register_wetdep_cpp();
}

namespace {
    struct Arguments {
        std::string config;
        std::string profile;
        std::string snapshot;
        std::string initial_snapshot;
        int steps = 1;
        double dt = 10.0;
        bool zero_initial_chemistry = false;
    };

    Arguments parse_arguments(int argc, char** argv) {
        Arguments args;
        for (int i = 1; i < argc; i += 2) {
            if (i + 1 >= argc)
                throw std::runtime_error("expected a value after " + std::string(argv[i]));
            const std::string key = argv[i];
            if (key == "--config")
                args.config = argv[i + 1];
            else if (key == "--met-profile")
                args.profile = argv[i + 1];
            else if (key == "--snapshot")
                args.snapshot = argv[i + 1];
            else if (key == "--initial-snapshot")
                args.initial_snapshot = argv[i + 1];
            else if (key == "--steps")
                args.steps = std::stoi(argv[i + 1]);
            else if (key == "--dt")
                args.dt = std::stod(argv[i + 1]);
            else if (key == "--initial-chemistry")
                args.zero_initial_chemistry = std::string(argv[i + 1]) == "zero";
            else
                throw std::runtime_error("unknown argument: " + key);
        }
        if (args.config.empty() || args.profile.empty() || args.snapshot.empty() || args.steps < 1 || args.dt <= 0.0)
            throw std::runtime_error("--config, --met-profile, --snapshot, and positive --steps are required");
        return args;
    }

    std::vector<double> columnize(const YAML::Node& values, int columns, int vertical) {
        if (!values.IsSequence())
            throw std::runtime_error("profile field has an unexpected vertical extent");
        std::vector<double> result(static_cast<std::size_t>(columns) * vertical);
        // Surface fields are supplied one value per column.  Vertical profiles
        // remain a single shared profile unless a runner is later extended with
        // per-column nested profiles.
        if (vertical == 1 && static_cast<int>(values.size()) == columns) {
            for (int column = 0; column < columns; ++column)
                result[column] = values[column].as<double>();
            return result;
        }
        if (static_cast<int>(values.size()) != vertical)
            throw std::runtime_error("profile field has an unexpected vertical extent");
        for (int column = 0; column < columns; ++column)
            for (int level = 0; level < vertical; ++level)
                // InteropField uses layout_left: column is contiguous, then the
                // vertical axis.  This is the same packing used by the Fortran
                // core and its science bridges.
                result[static_cast<std::size_t>(column) + static_cast<std::size_t>(columns) * level] =
                    values[level].as<double>();
        return result;
    }

    void require_status(int status, const std::string& operation) {
        if (status != CATCHEM_SUCCESS)
            throw std::runtime_error(operation + " failed with status " + std::to_string(status));
    }

    std::size_t native_index(int column, int level, int species_index, int columns, int levels) {
        // CATChem's bound state follows the legacy Fortran layout: column is the
        // contiguous axis, then level, then species.
        return static_cast<std::size_t>(column) +
               static_cast<std::size_t>(columns) *
                   (static_cast<std::size_t>(level) + static_cast<std::size_t>(levels) * species_index);
    }

    void write_json(const std::string& path, const std::vector<double>& concentration, int columns, int levels,
                    int species, int steps) {
        std::ofstream output(path);
        if (!output)
            throw std::runtime_error("unable to create snapshot " + path);
        output << std::setprecision(17);
        output << "{\n  \"schema_version\": 1,\n  \"snapshots\": [{\n";
        output << "    \"checkpoint\": \"after_timestep_" << steps << "\",\n    \"fields\": {\n";
        output << "      \"concentration\": {\"units\": \"model-native\", \"shape\": [" << columns << ", " << levels
               << ", " << species << "], \"species\": [";
        static const std::vector<std::string> species_names = {
            "so2", "h2o2",  "oh",    "no3",   "so4",   "dms",   "dms_in", "msa",   "bc1",   "bc2",   "oc1",
            "oc2", "dust1", "dust2", "dust3", "dust4", "dust5", "seas1",  "seas2", "seas3", "seas4", "seas5"};
        if (species != static_cast<int>(species_names.size()))
            throw std::runtime_error("candidate parity runner species count does not match the Default fixture");
        for (int i = 0; i < species; ++i)
            output << (i ? ", " : "") << "\"" << species_names[i] << "\"";
        output << "], \"values\": [";
        bool first = true;
        for (int column = 0; column < columns; ++column)
            for (int level = 0; level < levels; ++level)
                for (int species_index = 0; species_index < species; ++species_index) {
                    output << (first ? "" : ", ")
                           << concentration[native_index(column, level, species_index, columns, levels)];
                    first = false;
                }
        output << "]}\n    }\n  }]\n}\n";
    }
} // namespace

int main(int argc, char** argv) {
    try {
        const auto args = parse_arguments(argc, argv);
        register_candidate_processes();
        const YAML::Node profile = YAML::LoadFile(args.profile);
        const int columns = profile["grid"]["columns"].as<int>();
        const int levels = profile["grid"]["levels"].as<int>();
        void* core = catchem_core_create_from_config_with_grid(args.config.c_str(), columns, levels);
        if (!core) {
            char detail[512] = {};
            catchem_get_last_error(detail, sizeof(detail));
            throw std::runtime_error(std::string("candidate core construction failed: ") + detail);
        }
        void* state = catchem_core_get_state_manager(core);
        if (!state)
            throw std::runtime_error("candidate state manager unavailable");

        std::vector<std::vector<double>> owned_fields;
        for (const auto& item : profile["met_2d"]) {
            owned_fields.push_back(columnize(item.second, columns, 1));
            require_status(catchem_state_bind_met_2d_checked(state, item.first.as<std::string>().c_str(),
                                                             owned_fields.back().data(), columns, 1),
                           "bind 2-D met");
        }
        for (const auto& item : profile["met_3d"]) {
            owned_fields.push_back(columnize(item.second, columns, levels));
            require_status(catchem_state_bind_met_3d_checked(state, item.first.as<std::string>().c_str(),
                                                             owned_fields.back().data(), columns, levels, 1),
                           "bind 3-D met");
        }
        for (const auto& item : profile["met_interface"]) {
            owned_fields.push_back(columnize(item.second, columns, levels + 1));
            require_status(catchem_state_bind_met_3d_axis_checked(state, item.first.as<std::string>().c_str(),
                                                                  owned_fields.back().data(), columns, levels + 1, 1,
                                                                  1),
                           "bind interface met");
        }
        for (const auto& item : profile["met_soil"]) {
            const int soil_levels = static_cast<int>(item.second.size());
            owned_fields.push_back(columnize(item.second, columns, soil_levels));
            require_status(catchem_state_bind_met_3d_axis_checked(state, item.first.as<std::string>().c_str(),
                                                                  owned_fields.back().data(), columns, soil_levels, 1,
                                                                  2),
                           "bind soil met");
        }

        // The Default fixture has 22 ordered species. The native legacy runner
        // uses the same checked-in species file and writes this identical order.
        const int species = 22;
        std::vector<double> concentration(static_cast<std::size_t>(columns) * levels * species, 0.0);
        if (!args.initial_snapshot.empty()) {
            const YAML::Node initial = YAML::LoadFile(args.initial_snapshot);
            const YAML::Node field = initial["snapshots"][0]["fields"]["concentration"];
            const auto shape = field["shape"];
            const auto values = field["values"];
            if (shape[0].as<int>() != columns || shape[1].as<int>() != levels || shape[2].as<int>() != species ||
                static_cast<int>(values.size()) != columns * levels * species)
                throw std::runtime_error("initial snapshot concentration extent does not match the candidate grid");
            for (int col = 0; col < columns; ++col)
                for (int lev = 0; lev < levels; ++lev)
                    for (int sp = 0; sp < species; ++sp) {
                        const std::size_t canonical =
                            static_cast<std::size_t>(sp) +
                            static_cast<std::size_t>(species) *
                                (static_cast<std::size_t>(lev) + static_cast<std::size_t>(levels) * col);
                        concentration[native_index(col, lev, sp, columns, levels)] = values[canonical].as<double>();
                    }
        } else if (!args.zero_initial_chemistry) {
            // Deterministic, positive tracer state makes the parity comparison
            // sensitive to process tendencies instead of passing trivially on
            // an all-zero chemistry field.
            for (int col = 0; col < columns; ++col)
                for (int lev = 0; lev < levels; ++lev)
                    for (int sp = 0; sp < species; ++sp)
                        concentration[native_index(col, lev, sp, columns, levels)] =
                            1.0e-12 * (1.0 + 0.01 * lev + 0.001 * sp);
        }
        require_status(
            catchem_state_bind_unified_chemistry_checked(state, concentration.data(), columns, levels, species),
            "bind concentration");
        // Process bridges consume the state clock's timestep.  Set it to the
        // same 10-second step used by the legacy parity driver.
        require_status(catchem_state_set_time_checked(state, 2026, 1, 1, 0, 0, 0, 1, args.dt), "set parity clock");
        for (int step = 0; step < args.steps; ++step)
            require_status(catchem_core_run_timestep(core, args.dt), "run timestep");
        write_json(args.snapshot, concentration, columns, levels, species, args.steps);
        catchem_core_destroy(core);
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "candidate parity runner: " << error.what() << '\n';
        return 1;
    }
}
