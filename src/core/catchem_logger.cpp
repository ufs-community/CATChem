#include "catchem_logger.hpp"
#include <algorithm>
#include <cctype>
#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <unistd.h>

namespace catchem {

    bool Logger::should_color(int fd) {
        const char* no_color = std::getenv("NO_COLOR");
        if (no_color && no_color[0] != '\0') {
            return false;
        }
        return isatty(fd);
    }

    void Logger::debug(const StateManager* state, std::string_view message, ContextList context) {
        log(state, "DEBUG", Level::Debug, message, context);
    }

    void Logger::info(const StateManager* state, std::string_view message, ContextList context) {
        log(state, "INFO ", Level::Info, message, context);
    }

    void Logger::warn(const StateManager* state, std::string_view message, ContextList context) {
        log(state, "WARN ", Level::Warn, message, context);
    }

    void Logger::error(const StateManager* state, std::string_view message, ContextList context) {
        log(state, "ERROR", Level::Error, message, context);
    }

    std::optional<Logger::Level> Logger::level_from_string(std::string_view name) {
        std::string value(name);
        std::transform(value.begin(), value.end(), value.begin(),
                       [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
        if (value == "debug")
            return Level::Debug;
        if (value == "info")
            return Level::Info;
        if (value == "warn" || value == "warning")
            return Level::Warn;
        if (value == "error")
            return Level::Error;
        return std::nullopt;
    }

    int& Logger::override_level() {
        // Function-local static avoids static-initialisation-order problems and
        // stays valid across shared-library teardown.  -1 means "no override".
        static int override_value = -1;
        return override_value;
    }

    void Logger::set_level(Level level) {
        override_level() = static_cast<int>(level);
    }

    void Logger::clear_level() {
        override_level() = -1;
    }

    Logger::Level Logger::threshold() {
        // An explicit override (YAML configuration) always wins over the
        // environment variable, which is resolved once and cached.
        if (const int forced = override_level(); forced >= 0)
            return static_cast<Level>(forced);
        static const Level env_level = [] {
            const char* raw = std::getenv("CATCHEM_LOG_LEVEL");
            if (raw != nullptr && raw[0] != '\0') {
                if (const auto parsed = level_from_string(raw))
                    return *parsed;
            }
            return Level::Info;
        }();
        return env_level;
    }

    bool Logger::enabled(Level level) {
        return level >= threshold();
    }

    void Logger::log(const StateManager* state, std::string_view level, Level level_id, std::string_view message,
                     ContextList context) {
        if (!enabled(level_id))
            return;

        // 1. Get exact current UTC Timestamp
        auto now = std::chrono::system_clock::now();
        std::time_t now_time = std::chrono::system_clock::to_time_t(now);
        std::tm* utc_tm = std::gmtime(&now_time);

        std::ostringstream ss;
        ss << std::put_time(utc_tm, "%Y-%m-%d %H:%M:%S");
        std::string timestamp = ss.str();

        // 2. Format Level with ANSI Coloring
        int fd = (level == "ERROR") ? fileno(stderr) : fileno(stdout);
        bool color = should_color(fd);

        std::string colored_level(level);
        if (color) {
            if (level == "DEBUG")
                colored_level = "\033[36mDEBUG\033[0m"; // Cyan
            else if (level == "INFO ")
                colored_level = "\033[32mINFO \033[0m"; // Green
            else if (level == "WARN ")
                colored_level = "\033[33mWARN \033[0m"; // Yellow
            else if (level == "ERROR")
                colored_level = "\033[31mERROR\033[0m"; // Red
        }

        // 3. Assemble Service Name (exactly 15 chars, left-justified)
        std::string service = "catchem";
        service.append(15 - service.length(), ' ');

        // 4. Assemble Trace ID (exactly 8 chars)
        std::string trace = (state && !state->runtime_trace_id().empty()) ? state->runtime_trace_id() : "global  ";
        if (trace.length() < 8)
            trace.append(8 - trace.length(), ' ');

        // 5. Build full golden prefix
        std::ostringstream out;
        out << "[" << timestamp << "] [" << colored_level << "] [" << service << "] [" << trace << "] " << message;

        // 6. Append Key-Value Context dictionary
        if (context.size() > 0) {
            out << " |";
            for (const auto& [key, value] : context) {
                out << " " << key << "=" << value;
            }
        }
        out << "\n";

        // 7. Stream out cleanly
        if (level == "ERROR") {
            std::cerr << out.str() << std::flush;
        } else {
            std::clog << out.str() << std::flush;
        }
    }

} // namespace catchem
