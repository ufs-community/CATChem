#include "catchem_api.hpp"
#include "catchem_core.hpp"
#include "catchem_diagnostic_manager.hpp"
#include <cassert>
#include <stdexcept>
#include <vector>

int main() {
    void* core_handle = catchem_core_create(2, 3, 1);
    assert(core_handle);
    auto* core = static_cast<catchem::Core*>(core_handle);
    auto manager = core->get_diagnostic_manager();
    const std::vector<int> dims = {2, 1};
    const std::vector<catchem::SemanticAxis> axes = {catchem::SemanticAxis::Column, catchem::SemanticAxis::Singleton};
    manager->register_field_contract("instant", "instantaneous value", "1", catchem::DiagType::FIELD_2D, dims,
                                     catchem::DiagnosticPolicy::Instantaneous, 0.0, axes);
    manager->register_field_contract("accumulated", "timestep accumulation", "kg", catchem::DiagType::FIELD_2D, dims,
                                     catchem::DiagnosticPolicy::TimestepAccumulated, -1.0, axes);
    manager->register_field_contract("persistent", "persistent value", "m", catchem::DiagType::FIELD_2D, dims,
                                     catchem::DiagnosticPolicy::Persistent, 0.0, axes);

    static_cast<double*>(manager->get_host_write_pointer("instant"))[0] = 7.0;
    static_cast<double*>(manager->get_host_write_pointer("accumulated"))[0] = 8.0;
    static_cast<double*>(manager->get_host_write_pointer("persistent"))[0] = 9.0;
    manager->begin_timestep();
    assert(static_cast<const double*>(manager->get_host_read_pointer("instant"))[0] == 0.0);
    assert(static_cast<const double*>(manager->get_host_read_pointer("accumulated"))[0] == -1.0);
    assert(static_cast<const double*>(manager->get_host_read_pointer("persistent"))[0] == 9.0);
    assert(manager->get_field("instant")->latest_writer == catchem::LatestWriter::Synchronized);
    assert(manager->get_field("persistent")->latest_writer == catchem::LatestWriter::HostCurrent);

    for (int timestep = 2; timestep <= 3; ++timestep) {
        static_cast<double*>(manager->get_host_write_pointer("instant"))[0] = timestep;
        static_cast<double*>(manager->get_host_write_pointer("accumulated"))[0] += timestep;
        manager->begin_timestep();
        assert(manager->get_field("instant")->generation == static_cast<std::size_t>(timestep));
        assert(manager->get_field("persistent")->generation == static_cast<std::size_t>(timestep));
        assert(static_cast<const double*>(manager->get_host_read_pointer("persistent"))[0] == 9.0);
    }

    manager->register_field_contract("instant", "instantaneous value", "1", catchem::DiagType::FIELD_2D, dims,
                                     catchem::DiagnosticPolicy::Instantaneous, 0.0, axes);
    bool mismatch_rejected = false;
    try {
        manager->register_field_contract("instant", "different meaning", "1", catchem::DiagType::FIELD_2D, dims,
                                         catchem::DiagnosticPolicy::Instantaneous, 0.0, axes);
    } catch (const std::invalid_argument&) {
        mismatch_rejected = true;
    }
    assert(mismatch_rejected);

    // --- feature 013: unpack_labels contract (lenient subset) ---
    const std::vector<int> bin_dims = {2, 3};
    const std::vector<catchem::SemanticAxis> bin_axes = {catchem::SemanticAxis::Column,
                                                         catchem::SemanticAxis::Category};
    const std::vector<std::string> bin_labels = {"DUST1", "DUST2", "DUST3"};
    manager->register_field_contract("per_bin", "per bin", "kg", catchem::DiagType::FIELD_2D, bin_dims,
                                     catchem::DiagnosticPolicy::Instantaneous, 0.0, bin_axes, bin_labels);
    assert(manager->get_axes("per_bin") == bin_axes);
    assert(manager->get_unpack_labels("per_bin") == bin_labels);
    // A field with no packed dimension carries no labels.
    assert(manager->get_unpack_labels("instant").empty());
    // INV-9: re-registration with the same contract but differing labels is rejected.
    bool label_mismatch_rejected = false;
    try {
        manager->register_field_contract("per_bin", "per bin", "kg", catchem::DiagType::FIELD_2D, bin_dims,
                                         catchem::DiagnosticPolicy::Instantaneous, 0.0, bin_axes,
                                         {"DUST1", "DUST2", "WRONG"});
    } catch (const std::invalid_argument&) {
        label_mismatch_rejected = true;
    }
    assert(label_mismatch_rejected);
    // INV-8: the leading axis of every process diagnostic must be Column.
    bool non_column_rejected = false;
    try {
        manager->register_field_contract("bad_axis", "bad", "kg", catchem::DiagType::FIELD_2D, bin_dims,
                                         catchem::DiagnosticPolicy::Instantaneous, 0.0,
                                         {catchem::SemanticAxis::Level, catchem::SemanticAxis::Singleton});
    } catch (const std::invalid_argument&) {
        non_column_rejected = true;
    }
    assert(non_column_rejected);
    // --- feature 013 T021: strict packed-axis validation (INV-3..7) ---
    // Each rejected contract throws before mutating the manager, so the tests
    // are independent and share one manager.
    auto rejects = [&](const char* name, const std::vector<int>& d,
                      const std::vector<catchem::SemanticAxis>& a,
                      const std::vector<std::string>& l) {
        bool thrown = false;
        try {
            manager->register_field_contract(name, "d", "kg", catchem::DiagType::FIELD_2D, d,
                                             catchem::DiagnosticPolicy::Instantaneous, 0.0, a, l);
        } catch (const std::invalid_argument&) {
            thrown = true;
        }
        return thrown;
    };
    const std::vector<int> two = {2, 2};
    // INV-3: two packed axes cannot be unpacked by the writer.
    assert(rejects("inv3", {2, 2, 2},
                   {catchem::SemanticAxis::Column, catchem::SemanticAxis::Species, catchem::SemanticAxis::Category},
                   {"a", "b"}));
    // INV-4: a packed axis must carry exactly one label per slot.
    assert(rejects("inv4_short", two, {catchem::SemanticAxis::Column, catchem::SemanticAxis::Species}, {"a"}));
    assert(rejects("inv4_long", two, {catchem::SemanticAxis::Column, catchem::SemanticAxis::Species},
                   {"a", "b", "c"}));
    assert(rejects("inv4_none", two, {catchem::SemanticAxis::Column, catchem::SemanticAxis::Species}, {}));
    // INV-5: labels with nothing to label are a contract bug.
    assert(rejects("inv5", {2, 1}, {catchem::SemanticAxis::Column, catchem::SemanticAxis::Singleton}, {"a"}));
    // INV-6: labels must be NetCDF-name-safe.
    assert(rejects("inv6_digit", two, {catchem::SemanticAxis::Column, catchem::SemanticAxis::Species},
                   {"1bad", "ok"}));
    assert(rejects("inv6_dash", two, {catchem::SemanticAxis::Column, catchem::SemanticAxis::Species},
                   {"has-dash", "ok"}));
    assert(rejects("inv6_empty", two, {catchem::SemanticAxis::Column, catchem::SemanticAxis::Species},
                   {"", "ok"}));
    // INV-7: duplicate labels would unpack to the same variable name.
    assert(rejects("inv7", two, {catchem::SemanticAxis::Column, catchem::SemanticAxis::Species}, {"a", "a"}));
    // A fully valid packed contract still registers (positive control).
    assert(!rejects("inv_ok", two, {catchem::SemanticAxis::Column, catchem::SemanticAxis::Species}, {"a", "b"}));
    assert(manager->get_unpack_labels("inv_ok").size() == 2);
    // Determinism: get_registered_names() yields insertion order, not hash order.
    const auto names = manager->get_registered_names();
    assert(names.size() == 5);
    assert(names[0] == "instant" && names[1] == "accumulated" && names[2] == "persistent" && names[3] == "per_bin" &&
           names[4] == "inv_ok");

    assert(catchem_core_destroy_checked(core_handle) == CATCHEM_SUCCESS);
    return 0;
}
