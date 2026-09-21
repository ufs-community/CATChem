#include "catchem_physical_validation.hpp"
#include <algorithm>
#include <sstream>
namespace catchem {
    void PhysicalValidationReport::observe(std::string field, std::string rule, double value, std::size_t location,
                                           std::string correction) {
        auto found = std::find_if(issues_.begin(), issues_.end(),
                                  [&](const auto& issue) { return issue.field == field && issue.rule == rule; });
        if (found == issues_.end()) {
            issues_.push_back({std::move(field), std::move(rule), 0, value, value, {}, std::move(correction)});
            found = std::prev(issues_.end());
        }
        ++found->count;
        found->observed_min = std::min(found->observed_min, value);
        found->observed_max = std::max(found->observed_max, value);
        if (found->locations.size() < 16)
            found->locations.push_back(location);
    }
    std::string PhysicalValidationReport::format() const {
        std::ostringstream out;
        for (const auto& issue : issues_)
            out << issue.field << " [" << issue.rule << "]: " << issue.count << " invalid values, range "
                << issue.observed_min << ".." << issue.observed_max
                << (issue.correction.empty() ? "" : "; " + issue.correction) << '\n';
        return out.str();
    }
} // namespace catchem
