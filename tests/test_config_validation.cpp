#include "catchem_api.hpp"
#include "catchem_config_manager.hpp"
#include "catchem_test_config.hpp"
#include <cassert>
#include <string>

int main() {
    const std::string fixtures = std::string(catchem::test::TEST_DIR) + "/fixtures/";
    catchem::ConfigManager valid;
    valid.load_from_file(fixtures + "platform_integrity_valid.yml");
    valid.load_species_file(fixtures + "platform_integrity_species.yml");
    assert(!valid.validate().has_errors());
    assert(valid.data.species[0].name == "CO");
    assert(valid.data.species[1].name == "O3");

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

    void* core = reinterpret_cast<void*>(1);
    assert(catchem_core_create_from_config_checked((fixtures + "platform_integrity_invalid.yml").c_str(), &core) !=
           CATCHEM_SUCCESS);
    assert(core == nullptr);
    return 0;
}
