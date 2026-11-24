#pragma once

#include <algorithm>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace catchem {

    enum class IssueSeverity { Warning, Error };

    struct ValidationIssueDetail {
        IssueSeverity severity = IssueSeverity::Error;
        std::string category;
        std::string path;
        std::string message;
        std::string correction;
    };

    class ValidationIssueReport {
    public:
        void add(ValidationIssueDetail issue) { issues_.push_back(std::move(issue)); }
        bool has_errors() const noexcept {
            return std::any_of(issues_.begin(), issues_.end(),
                               [](const auto& issue) { return issue.severity == IssueSeverity::Error; });
        }
        bool empty() const noexcept { return issues_.empty(); }
        const std::vector<ValidationIssueDetail>& issues() const noexcept { return issues_; }
        std::string format() const {
            std::ostringstream output;
            for (const auto& issue : issues_) {
                output << (issue.severity == IssueSeverity::Error ? "error" : "warning") << " [" << issue.category
                       << "] " << issue.path << ": " << issue.message;
                if (!issue.correction.empty())
                    output << " (" << issue.correction << ")";
                output << '\n';
            }
            return output.str();
        }

    private:
        std::vector<ValidationIssueDetail> issues_;
    };

} // namespace catchem
