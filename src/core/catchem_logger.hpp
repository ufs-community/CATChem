#pragma once
#include "catchem_state_manager.hpp"
#include <initializer_list>
#include <optional>
#include <string_view>
#include <utility>

namespace catchem {

    class Logger {
    public:
        using ContextList = std::initializer_list<std::pair<std::string_view, std::string_view>>;

        // Minimum severity actually emitted.  An explicit override set with
        // set_level() (the YAML configuration path, simulation/verbose/log_level)
        // always wins; otherwise the CATCHEM_LOG_LEVEL environment variable
        // (DEBUG | INFO | WARN | ERROR, case-insensitive) is used, and an unset
        // or unrecognized value selects INFO so production runs stay quiet.
        enum class Level { Debug = 0, Info = 1, Warn = 2, Error = 3 };

        static void debug(const StateManager* state, std::string_view message, ContextList context = {});
        static void info(const StateManager* state, std::string_view message, ContextList context = {});
        static void warn(const StateManager* state, std::string_view message, ContextList context = {});
        static void error(const StateManager* state, std::string_view message, ContextList context = {});

        // Returns true when a message of the given level would be emitted.
        // Call sites with expensive context preparation should check this
        // first and skip the whole block when it returns false.
        static bool enabled(Level level);

        // Parse a level name case-insensitively; "warning" is accepted as an
        // alias for "warn".  Returns std::nullopt for unrecognized names.
        static std::optional<Level> level_from_string(std::string_view name);

        // Force the emission threshold, overriding CATCHEM_LOG_LEVEL.
        static void set_level(Level level);

        // Drop the explicit override; CATCHEM_LOG_LEVEL applies again.
        static void clear_level();

    private:
        static void log(const StateManager* state, std::string_view level, Level level_id, std::string_view message,
                        ContextList context);
        static bool should_color(int fd);
        static Level threshold();
        static int& override_level();
    };

} // namespace catchem
