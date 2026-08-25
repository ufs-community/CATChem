#include "catchem_logger.hpp"
#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
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
        log(state, "DEBUG", message, context);
    }

    void Logger::info(const StateManager* state, std::string_view message, ContextList context) {
        log(state, "INFO ", message, context);
    }

    void Logger::warn(const StateManager* state, std::string_view message, ContextList context) {
        log(state, "WARN ", message, context);
    }

    void Logger::error(const StateManager* state, std::string_view message, ContextList context) {
        log(state, "ERROR", message, context);
    }

    void Logger::log(const StateManager* state, std::string_view level, std::string_view message, ContextList context) {
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
        std::string trace = (state && !state->trace_id.empty()) ? state->trace_id : "global  ";
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
