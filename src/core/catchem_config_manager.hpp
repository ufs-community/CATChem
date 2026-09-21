#pragma once
#include "catchem_physical_validation.hpp"
#include <map>
#include <string>
#include <string_view>
#include <vector>
#include <yaml-cpp/yaml.h>

namespace catchem {

    enum class ValidationSeverity { Warning, Error };

    struct ValidationIssue {
        ValidationSeverity severity = ValidationSeverity::Error;
        std::string category;
        std::string path;
        std::string message;
        std::string suggested_correction;
    };

    struct ValidationReport {
        std::vector<ValidationIssue> issues;
        bool has_errors() const;
        std::string format() const;
    };

    /// @brief Basic simulation file and verbosity settings from YAML.
    struct SimulationConfig {
        std::string name;
        std::string start_date;
        std::string end_date;
        std::string species_filename;
        std::string emission_filename;
        bool verbose_enabled = false;
        std::string log_level; ///< simulation/verbose/log_level, as written in the YAML (empty when unset)
    };

    struct RuntimeConfig {
        int nx = 1;
        int ny = 1;
        int nz = 1;
        double dt = 3600.0;
        int nsteps = 1;
    };

    /// @brief Grid dimensions and static grid settings from YAML.
    struct GridConfig {
        int number_of_levels = 1;
        int number_of_soil_layers = 0;
    };

    /// @brief Runtime timestep settings from YAML.
    struct TimestepConfig {
        int transport_timestep_in_s = 0;
        int chemistry_timestep_in_s = 0;
    };

    /// @brief Diagnostic output settings from YAML.
    struct DiagnosticOutputConfig {
        bool enabled = false;
        std::string directory;
        std::string prefix;
        int frequency = 0;
        std::string format;
        int compress_lev = 0;
        bool process_diagnostics = false;
        std::vector<std::string> diag_list;
        /// Run-level NetCDF global attributes (institution, references, ...).
        /// Written to every diagnostic file; user entries override core defaults
        /// on key collision (spec FR-011).
        std::map<std::string, std::string> attributes;
    };

    /// @brief Diagnostic collection settings from YAML.
    struct DiagnosticCollectionConfig {
        bool enabled = false;
        int buffer_size = 0;
    };

    /// @brief Top-level diagnostic settings from YAML.
    struct DiagnosticsConfig {
        DiagnosticOutputConfig output;
        DiagnosticCollectionConfig collection;
    };

    /// @brief Process activation and nested process settings from YAML.
    struct ProcessConfig {
    private:
        YAML::Node settings_node;
        friend class ConfigManager;

    public:
        bool activate = false;
        bool diagnostics = false;
        std::string scheme;
        std::vector<std::string> diag_species;

        void set_settings_node(const YAML::Node& node) { settings_node = node; }

        bool get_bool(std::string_view key, bool default_val = false) const;
        double get_double(std::string_view key, double default_val = 0.0) const;
        int get_int(std::string_view key, int default_val = 0) const;
        std::string get_string(std::string_view key, std::string_view default_val = "") const;
        std::vector<double> get_vector(std::string_view key) const;

        /// @brief List every nested option path as "<scheme>.<key>".
        ///
        /// Framework keys (activate, diagnostics, scheme, gas_scheme,
        /// aero_scheme, diag_species) are excluded because the process layer
        /// consumes them directly.  Process registration supplies a
        /// validator that checks these paths against the scheme's accepted
        /// options so a typo or removed parameter fails at initialization
        /// instead of silently keeping its compiled default.
        std::vector<std::string> option_paths() const;
    };

    /// @brief Species metadata loaded from a CATChem species YAML file.
    struct SpeciesConfig {
        std::string name;
        std::string long_name;
        std::string description;
        std::vector<std::string> aliases;
        std::vector<std::string> roles;

        bool is_gas = false;
        bool is_aerosol = false;
        bool is_tracer = false;
        bool is_advected = true;
        bool is_drydep = false;
        bool is_wetdep = false;
        bool is_photolysis = false;
        bool is_gocart_aero = false;
        bool is_dust = false;
        bool is_seasalt = false;
        bool is_hydrophilic = true;

        double molecular_weight_kg_mol = 0.0;
        double mw_g = 0.0;
        double density = 0.0;
        double radius = 0.0;
        double lower_radius = 0.0;
        double upper_radius = 0.0;
        double viscosity = 0.0;

        // Dry deposition parameters
        double dd_f0 = 0.0;
        double dd_hstar = 0.0;
        double dd_DvzAerSnow = 0.0;
        double dd_DvzMinVal_snow = 0.0;
        double dd_DvzMinVal_land = 0.0;

        // Wet deposition parameters
        double henry_k0 = 0.0;
        double henry_cr = 0.0;
        double henry_pKa = 0.0;
        double wd_retfactor = 0.0;
        bool wd_LiqAndGas = false;
        double wd_convfacI2G = 0.0;
        std::vector<double> wd_rainouteff = {0.0, 0.0, 0.0};
        double wd_reevap_frac = 0.5;

        // Chemical loss rate and background volume-mixing ratio
        double t_chem_loss = -1.0;
        double BackgroundVV = 1.0e-20;
        std::string mie_name;
    };

    /// @brief Mapping for one external emission source field.
    struct EmissionFieldMapping {
        std::string long_name;
        std::string units;
        std::vector<double> scale;
        std::vector<std::string> map;
    };

    /// @brief Emission mapping category from a CATChem emission YAML file.
    struct EmissionCategoryMapping {
        std::map<std::string, EmissionFieldMapping> fields;
    };

    /// @brief Aerosol optics (Mie) table inputs from the top-level "mie:" section.
    ///
    /// Mirrors the legacy upstream/develop configuration: "mie.directory" holds the
    /// table directory and "mie.files" maps an aerosol type code (SS, DU, BC, ...) to
    /// the NetCDF optics file supplying it.  Species bind to a type through their
    /// "__mie_name" attribute.  "files" keeps the YAML declaration order so every rank
    /// loads the same tables in the same sequence (determinism).
    struct MieConfig {
        std::string directory = "./";
        std::vector<std::pair<std::string, std::string>> files;
    };

    struct ConfigData {
        SimulationConfig simulation;
        RuntimeConfig runtime;
        GridConfig grid;
        TimestepConfig timesteps;
        DiagnosticsConfig diagnostics;
        std::string species_filename; ///< simulation:species_filename, as written in the YAML
        std::vector<std::string> active_processes;
        std::map<std::string, ProcessConfig> processes;
        std::vector<SpeciesConfig> species;
        std::string mechanism_identity;
        std::vector<std::string> mechanism_capabilities;
        PhysicalValidationPolicy physical_validation_policy = PhysicalValidationPolicy::Reject;
        std::map<std::string, EmissionCategoryMapping> emission_mappings;
        MieConfig mie;
    };

    class ConfigManager {
    private:
        YAML::Node root_node;

    public:
        // When false, load_from_file() suppresses the stdout echo of the parsed
        // YAML.  The NUOPC driver disables it on non-root PETs so the effective
        // configuration is printed once per run instead of once per PET.
        // Defaults to true: standalone and test runs keep the historical
        // behavior of printing exactly what the core parsed.
        static bool echo_config_to_stdout;

        ConfigData data;
        bool is_loaded = false;
        std::string config_file_path;
        ValidationReport validation_report;

        ConfigManager() = default;
        void load_from_file(const std::string& filename);
        void load_species_file(const std::string& filename);
        void load_emission_mapping_file(const std::string& filename);
        const ValidationReport& validate(bool strict = true);
        void validate_or_throw(bool strict = true);

        // Safe path-based queries
        bool get_bool(std::string_view path, bool default_val = false) const;
        double get_double(std::string_view path, double default_val = 0.0) const;
        int get_int(std::string_view path, int default_val = 0) const;
        std::string get_string(std::string_view path, std::string_view default_val = "") const;
        std::vector<std::string> get_string_list(std::string_view path) const;

        // Structured process / emission queries
        bool is_process_active(std::string_view process_name) const;
        bool is_category_active(std::string_view category_name) const;
        std::string find_process_file_setting(std::string_view process_name) const;
    };

} // namespace catchem
