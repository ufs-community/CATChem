# Process Development

This section covers developing new processes and schemes for CATChem. See the **[CATChem User Guide](../../user-guide/index.md#process-documentation)** for a description of all processes available in CATChem.

## Quick Links

- **[Architecture Overview](architecture.md)** - Understanding process design patterns
- **[Process Generator](process-generator.md)** - Complete guide to using the automated process generator
- **[Creating Custom Processes](creating.md)** - Manual process development
- **[Templates and Patterns](templates.md)** - Code templates and best practices
- **[Testing Processes](testing.md)** - Testing strategies and frameworks

## Overview

CATChem processes are modular components that implement specific atmospheric transport, chemical, emission, or loss schemes. Each process follows a standardized interface and lifecycle.

## Process Architecture

All processes inherit from C++ `catchem::ProcessInterface` and must implement `init`, `run`, and `finalize`. Process C++ wrappers interface with pure Fortran science schemes using C-interoperable Fortran ScienceBridges (`BIND(C)`).

```cpp
#include "catchem_process_interface.hpp"

namespace catchem {

    class MyProcessProcess : public ProcessInterface {
    public:
        std::string active_scheme;

        MyProcessProcess();

        std::string get_name() const override { return "myprocess"; }
        void init(std::shared_ptr<StateManager> state) override;
        void run(std::shared_ptr<StateManager> state) override;
        void finalize() override;
    };

} // namespace catchem
```

## Creating New Processes

To create a new process, follow the [Creating Custom Processes](creating.md) guide.

## Scheme Development

Schemes implement specific algorithms within a process. See the `seasalt` process for an example of how to structure schemes.

## Testing Processes

All new processes must include unit and integration tests. See the [Testing Processes](testing.md) guide for more information.

## Documentation

All new processes must be documented using Doxygen-style comments.

## See Also

- [Process Architecture](architecture.md)
- [Testing Processes](testing.md)
- [Process Templates](templates.md)
