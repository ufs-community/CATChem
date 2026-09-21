# Design Specification: Standardized Columnar Logging Architecture

This document specifies the design, implementation, and integration of a standardized, human-readable logging architecture across CATChem.

---

## 1. Architectural Philosophy

We explicitly reject JSON logging in favor of strict, columnar plain-text logging to maximize human readability via standard terminal streams and log files. Clear visual alignment and searchability (via tools like grep) are primary design goals.

### The Golden Log Format
Every single log emitted by this application MUST strictly adhere to this exact padded template:
```
[TIMESTAMP] [LEVEL] [SERVICE] [TRACE_ID] Message | key=value
```

*Example:*
```
[2026-07-28 10:14:02] [INFO ] [core-api      ] [req-12b4] User login successful | user_id=42
```

---

## 2. Core Requirements

1. **Timestamp**: `YYYY-MM-DD HH:MM:SS` format in UTC time.
2. **Level**: Exactly 5 characters, space-padded (`DEBUG`, `INFO `, `WARN `, `ERROR`).
3. **Service**: Left-justified, padded to a fixed width of exactly 15 characters (e.g. `"core           "`).
4. **Trace ID**: A 6-8 character alphanumeric correlation ID managed in the `StateManager` execution lifecycle (e.g., `req-12b4`).
5. **Message**: A clear, static string describing the event.
6. **Context (Key=Value)**: If a developer passes key-value pairs, they are appended to the end of the log line separated by a ` | `, formatted exactly as `key1=value1 key2=value2`.

---

## 3. Logger Central Wrapper Interface (`catchem_logger.hpp`)

The logger will reside in the `catchem` namespace. It will accept the state manager pointer, a format message, and an initializer list of string-view pairs for context tagging:

```cpp
#pragma once
#include "catchem_state_manager.hpp"
#include <string_view>
#include <initializer_list>
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
        static void log(const StateManager* state, std::string_view level, std::string_view message, ContextList context);
        static bool should_color();
    };

} // namespace catchem
```

---

## 4. StateManager Integration & Trace ID Propagation

A public `trace_id` property will be added to `catchem::StateManager`. The trace ID is generated automatically as a short 8-character random alphanumeric string when the state is instantiated, representing the simulation timeline or runtime execution context:

```cpp
// Sourced dynamically from StateManager
namespace catchem {
    class StateManager {
    public:
        std::string trace_id; // e.g. "req-12b4"
        // ... rest of fields ...
    };
}
```

If `state` is null, the logger will default the Trace ID column to `[global  ]`.

---

## 5. ANSI Color Rules and Stripping
* Visual colors are applied to the `[LEVEL]` tag on local interactive terminal streams.
* Colors are completely stripped if the `NO_COLOR` environment variable is set or if standard output streams are redirected to log files (`!isatty(fileno(stdout))` and `!isatty(fileno(stderr))`).

---

## 6. Verification and Testing Strategy
* **Unit Tests**: Add tests verifying formatting layout, zero-allocation context processing, ANSI color stripping, and Trace ID retrieval.
* **CTests Integration**: Confirm that standard CTest runs pass with the newly integrated logger calls.
