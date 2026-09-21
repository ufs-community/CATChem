#include "catchem_error.hpp"
#include <cassert>
#include <stdexcept>
#include <string>

int main() {
    int sentinel = 42;

    // Should not throw for a valid pointer.
    catchem::require_field_pointer("UnitTest", "VALID_FIELD", &sentinel);

    bool threw = false;
    try {
        catchem::require_field_pointer("UnitTest", "MISSING_FIELD", nullptr);
    } catch (const std::runtime_error& error) {
        threw = true;
        const std::string message = error.what();
        assert(message == "FATAL ERROR: UnitTest process missing required field MISSING_FIELD");
    }

    assert(threw && "require_field_pointer must throw for null pointers");
    return 0;
}
