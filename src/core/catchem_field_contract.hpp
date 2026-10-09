#pragma once

#include <algorithm>
#include <cctype>
#include <cstddef>
#include <string>
#include <vector>

namespace catchem {

    enum class DataflowStatus : int {
        Success = 0,
        NullArgument = 1,
        MissingField = 2,
        RankMismatch = 3,
        ExtentMismatch = 4,
        InvalidIndex = 5,
        StaleGeneration = 6,
        DuplicateMapping = 7,
        InvalidState = 8,
        InternalError = 9
    };

    enum class SemanticAxis { Column, Level, Interface, SoilLayer, Species, Category, Singleton };
    enum class PersistencePolicy { Timestep, Persistent };
    enum class AvailabilityState { Unavailable, Current, Stale };
    enum class LatestWriter { Uninitialized, HostCurrent, DeviceCurrent, Synchronized };
    enum class AccessIntent { Read, Write, ReadWrite };
    enum class ExecutionSpaceIntent { Host, Device, Either };
    enum class FieldRequirement { Required, Optional };

    inline std::string canonicalize_field_identity(std::string name) {
        std::transform(name.begin(), name.end(), name.begin(),
                       [](unsigned char value) { return static_cast<char>(std::toupper(value)); });
        return name;
    }

    struct FieldContract {
        std::string canonical_name;
        std::string units;
        std::vector<std::size_t> extents;
        std::vector<SemanticAxis> axes;
        std::vector<std::string> aliases;
        PersistencePolicy persistence = PersistencePolicy::Timestep;

        bool structurally_valid() const noexcept {
            if (canonical_name.empty() || units.empty() || extents.empty() || extents.size() != axes.size())
                return false;
            for (const auto extent : extents) {
                if (extent == 0)
                    return false;
            }
            return true;
        }

        bool compatible_with(const FieldContract& other) const noexcept {
            return canonicalize_field_identity(canonical_name) == canonicalize_field_identity(other.canonical_name) &&
                   units == other.units && extents == other.extents && axes == other.axes &&
                   persistence == other.persistence;
        }
    };

    struct FieldAccessContract {
        std::string canonical_name;
        std::string units;
        std::vector<SemanticAxis> axes;
        PersistencePolicy persistence = PersistencePolicy::Timestep;
        FieldRequirement requirement = FieldRequirement::Required;
        AccessIntent access = AccessIntent::Read;
        ExecutionSpaceIntent execution_space = ExecutionSpaceIntent::Either;
        bool produced = false;

        bool reads() const noexcept { return access == AccessIntent::Read || access == AccessIntent::ReadWrite; }
        bool writes() const noexcept { return access == AccessIntent::Write || access == AccessIntent::ReadWrite; }
    };

} // namespace catchem
