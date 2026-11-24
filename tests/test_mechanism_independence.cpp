#include "catchem_core_architecture_test_helpers.hpp"
#include <algorithm>
#include <cassert>
#include <stdexcept>
#include <vector>

int main() {
    auto mechanism = catchem::test::synthetic_mechanism({"xenon_a", "radical_b", "aerosol_c"}, "unfamiliar");
    mechanism.species[1].aliases = {"rb"};
    mechanism.species[1].roles = {"photolysis.ozone"};
    mechanism.capabilities.insert(catchem::canonical_species_name("photolysis"));
    mechanism.rebuild_index();
    assert(mechanism.index_of("RADICAL_B") == 1);
    assert(mechanism.index_of("rb") == 1);
    assert(mechanism.index_for_role("PHOTOLYSIS.OZONE") == 1);
    assert(mechanism.has_capability("PHOTOLYSIS"));

    std::reverse(mechanism.species.begin(), mechanism.species.end());
    mechanism.rebuild_index();
    assert(mechanism.species[mechanism.index_of("radical_b")].short_name == "radical_b");

    bool duplicate_failed = false;
    try {
        (void)catchem::test::synthetic_mechanism({"same", "SAME"});
    } catch (const std::invalid_argument&) {
        duplicate_failed = true;
    }
    assert(duplicate_failed);

    bool role_failed = false;
    try {
        (void)mechanism.index_for_role("missing.role");
    } catch (const std::out_of_range&) {
        role_failed = true;
    }
    assert(role_failed);

    std::vector<std::string> many;
    for (int i = 0; i < 2000; ++i)
        many.push_back("species_" + std::to_string(i));
    assert(catchem::test::synthetic_mechanism(std::move(many)).species.size() == 2000);
    return 0;
}
