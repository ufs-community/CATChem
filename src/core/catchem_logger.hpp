#pragma once
#include "catchem_state_manager.hpp"
#include <initializer_list>
#include <string_view>
#include <utility>

namespace catchem {

    class Logger {
    public:
        using ContextList = std::initializer_list<std::pair<std::string_view, std::string_view>>;

        static void debug(const StateManager* state, std::string_view message, ContextList context = {});
        static void info(const StateManager* state, std::string_view message, ContextList context = {});
        static void warn(const StateManager* state, std::string_view message, ContextList context = {});
        static void error(const StateManager* state, std::string_view message, ContextList context = {});

    private:
        static void log(const StateManager* state, std::string_view level, std::string_view message,
                        ContextList context);
        static bool should_color(int fd);
    };

} // namespace catchem
