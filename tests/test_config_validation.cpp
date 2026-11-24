#include "catchem_api.hpp"
#include "catchem_config_manager.hpp"
#include "catchem_test_config.hpp"
#include <cassert>
#include <stdexcept>
#include <string>

int main() {
    const std::string fixtures = std::string(catchem::test::TEST_DIR) + "/fixtures/";
    catchem::ConfigManager valid;
    valid.load_from_file(fixtures + "platform_integrity_valid.yml");
    valid.load_species_file(fixtures + "platform_integrity_species.yml");
    assert(!valid.validate().has_errors());
    assert(valid.data.species[0].name == "CO");
    assert(valid.data.species[1].name == "O3");

    // A watch species that is not in the mechanism is a warning, not an error.
    valid.data.simulation.watch_species = {"co", "nope"};
    const auto& watch_report = valid.validate();
    assert(!watch_report.has_errors());
    bool warned = false;
    for (const auto& issue : watch_report.issues)
        if (issue.severity == catchem::ValidationSeverity::Warning &&
            issue.message.find("watch species is absent from active mechanism: nope") != std::string::npos)
            warned = true;
    assert(warned);
    valid.data.simulation.watch_species.clear();

    // Scalar-parsing contract (strict_booleans.yml): booleans are true/false
    // only; null or absent means the default; a legacy spelling or a wrong
    // type is an error, never a silent fallback.
    const auto throws_invalid = [](const auto& call) {
        try {
            call();
        } catch (const std::invalid_argument& e) {
            return std::string(e.what()).find("Configuration value is not") != std::string::npos;
        }
        return false;
    };
    catchem::ConfigManager scalars;
    scalars.load_from_file(fixtures + "strict_booleans.yml");
    assert(scalars.get_bool("simulation/strict_true", false) == true);
    assert(scalars.get_bool("simulation/strict_false", true) == false);
    assert(throws_invalid([&] { scalars.get_bool("simulation/legacy_yes", false); }));
    assert(throws_invalid([&] { scalars.get_bool("simulation/legacy_one", false); }));
    assert(scalars.get_bool("simulation/null_value", true) == true);
    assert(scalars.get_bool("simulation/missing", true) == true);
    assert(scalars.get_int("simulation/nx", 7) == 2);
    assert(throws_invalid([&] { scalars.get_int("simulation/float_as_int", 7); }));
    assert(throws_invalid([&] { scalars.get_double("simulation/not_a_number", 2.5); }));
    assert(scalars.get_double("simulation/timestep", 2.5) == 60.0);
    assert(scalars.get_string("simulation/null_value", "dflt") == "dflt");
    assert(scalars.get_string("simulation/not_a_number", "dflt") == "abc");
    assert(throws_invalid([&] { scalars.get_string("simulation", "dflt"); }));
    const auto& dust = scalars.data.processes.at("dust");
    assert(dust.activate);
    assert(throws_invalid([&] { dust.get_bool("legacy_on", false); }));
    assert(throws_invalid([&] { dust.get_bool("legacy_zero", true); }));
    assert(dust.get_bool("strict_false", true) == false);
    assert(dust.get_bool("missing", true) == true);
    // The activate switch follows the same contract at parse time.
    catchem::ConfigManager legacy_activate;
    assert(throws_invalid([&] { legacy_activate.load_from_file(fixtures + "strict_booleans_bad_activate.yml"); }));

    catchem::ConfigManager invalid;
    invalid.load_from_file(fixtures + "platform_integrity_invalid.yml");
    catchem::SpeciesConfig first;
    first.name = "O3";
    catchem::SpeciesConfig duplicate;
    duplicate.name = "o3";
    invalid.data.species = {first, duplicate};
    catchem::EmissionFieldMapping mapping;
    mapping.map = {"MISSING", "O3"};
    mapping.scale = {1.0};
    invalid.data.emission_mappings["anthro"].fields["x"] = mapping;
    const auto& report = invalid.validate();
    assert(report.has_errors());
    const std::string text = report.format();
    for (const char* required : {"simulation/species_filename", "simulation/nx", "grid/number_of_levels",
                                 "diagnostics/output/frequency", "unknown process", "duplicate species",
                                 "different lengths", "absent from active mechanism", "unknown_section"})
        assert(text.find(required) != std::string::npos);

    catchem::ConfigManager met_mapping;
    met_mapping.load_from_file(fixtures + "platform_integrity_valid.yml");
    met_mapping.load_species_file(fixtures + "platform_integrity_species.yml");
    catchem::EmissionFieldMapping met_field;
    met_field.map = {"MET_CLAYFRAC"};
    met_field.scale = {1.0};
    met_mapping.data.emission_mappings["dust"].fields["clayfrac"] = met_field;
    assert(!met_mapping.validate().has_errors());

    // Every active emission sector must have a mapping.  Sector discovery is
    // deliberately generic: this catches fire today and any future configured
    // sector without adding names to the validator.
    catchem::ConfigManager missing_sector_mapping;
    missing_sector_mapping.load_from_file(fixtures + "platform_integrity_missing_emission_mapping.yml");
    missing_sector_mapping.load_species_file(fixtures + "platform_integrity_species.yml");
    const auto& missing_mapping_report = missing_sector_mapping.validate();
    const std::string missing_mapping_text = missing_mapping_report.format();
    assert(missing_mapping_text.find("processes/extemis/fire") != std::string::npos);
    assert(missing_mapping_text.find("active emission sector has no mapping entry: fire") != std::string::npos);

    // The production-style default configuration is the concrete all-sector
    // fixture.  Validation itself discovers sectors from YAML, so this test
    // expands automatically when a new enabled category is added.
    const std::string default_config_dir = std::string(catchem::test::TEST_DIR) + "/Configs/Default/";
    catchem::ConfigManager all_sector_config;
    all_sector_config.load_from_file(default_config_dir + "CATChem_new_config.yml");
    all_sector_config.load_species_file(default_config_dir + "CATChem_species.yml");
    all_sector_config.load_emission_mapping_file(default_config_dir + "CATChem_emission.yml");
    const auto& all_sector_report = all_sector_config.validate();
    assert(all_sector_report.format().find("active emission sector has no mapping entry") == std::string::npos);
    // The shipped default config must validate cleanly with the new emission
    // temporal options absent (they default to daily_hold=false and
    // monthly_anchor='climatological').
    assert(!all_sector_report.has_errors());

    // An unknown monthly_anchor value must fail loudly and list the allowed
    // values (feature 014, US4, FR-018).
    catchem::ConfigManager bogus_anchor;
    bogus_anchor.load_from_file(fixtures + "emission_bogus_monthly_anchor.yml");
    const auto& bogus_report = bogus_anchor.validate();
    assert(bogus_report.has_errors());
    const std::string bogus_text = bogus_report.format();
    assert(bogus_text.find("unknown monthly_anchor value: bogus") != std::string::npos);
    assert(bogus_text.find("climatological") != std::string::npos && bogus_text.find("file") != std::string::npos);

    // A valid 'file' anchor must not raise the enum error.
    catchem::ConfigManager file_anchor;
    file_anchor.load_from_file(fixtures + "emission_file_monthly_anchor.yml");
    assert(file_anchor.data.emission_categories.at("anthro1").monthly_anchor == "file");
    assert(file_anchor.data.emission_categories.at("anthro1").daily_hold == false);
    assert(file_anchor.validate().format().find("unknown monthly_anchor") == std::string::npos);

    void* core = reinterpret_cast<void*>(1);
    assert(catchem_core_create_from_config_checked((fixtures + "platform_integrity_invalid.yml").c_str(), &core) !=
           CATCHEM_SUCCESS);
    assert(core == nullptr);
    return 0;
}
