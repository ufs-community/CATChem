#pragma once
#include <cstddef>
#include <limits>
#include <string>
#include <vector>
namespace catchem {
    enum class PhysicalValidationPolicy { Reject, WarnAndClamp, CountAndContinue };
    struct PhysicalIssue {
        std::string field, rule;
        std::size_t count = 0;
        double observed_min = std::numeric_limits<double>::infinity();
        double observed_max = -std::numeric_limits<double>::infinity();
        std::vector<std::size_t> locations;
        std::string correction;
    };
    class PhysicalValidationReport {
    public:
        void clear() { issues_.clear(); }
        void observe(std::string field, std::string rule, double value, std::size_t location,
                     std::string correction = {});
        bool empty() const noexcept { return issues_.empty(); }
        std::size_t issue_count() const noexcept { return issues_.size(); }
        const std::vector<PhysicalIssue>& issues() const noexcept { return issues_; }
        std::string format() const;

    private:
        std::vector<PhysicalIssue> issues_;
    };
} // namespace catchem
