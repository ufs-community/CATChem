#include "catchem_api.hpp"
#include "catchem_core.hpp"
#include "catchem_process_registry.hpp"
#include <cassert>
#include <iostream>

extern "C" {
void catchem_register_settling_cpp();
void catchem_register_drydep_cpp();
void catchem_register_seasalt_cpp();
void catchem_register_dust_cpp();
void catchem_register_wetdep_cpp();
void catchem_register_so4chem_cpp();
void catchem_register_carbchem_cpp();
}

int main(int argc, char* argv[]) {
    Kokkos::initialize(argc, argv);
    {
        std::cout << "\n==========================================" << std::endl;
        std::cout << "RUNNING TEST: ProcessRegistry and C-API Process Registration" << std::endl;
        std::cout << "==========================================\n" << std::endl;

        // 1. Register processes via C-API handlers
        catchem_register_settling_cpp();
        catchem_register_drydep_cpp();
        catchem_register_seasalt_cpp();
        catchem_register_dust_cpp();
        catchem_register_wetdep_cpp();
        catchem_register_so4chem_cpp();
        catchem_register_carbchem_cpp();

        auto& registry = catchem::ProcessRegistry::get_instance();

        // 2. Verify all 7 core physical processes are registered in the ProcessRegistry
        assert(registry.has_process("settling") && "ProcessRegistry should contain 'settling'");
        assert(registry.has_process("drydep") && "ProcessRegistry should contain 'drydep'");
        assert(registry.has_process("seasalt") && "ProcessRegistry should contain 'seasalt'");
        assert(registry.has_process("dust") && "ProcessRegistry should contain 'dust'");
        assert(registry.has_process("wetdep") && "ProcessRegistry should contain 'wetdep'");
        assert(registry.has_process("so4chem") && "ProcessRegistry should contain 'so4chem'");
        assert(registry.has_process("carbchem") && "ProcessRegistry should contain 'carbchem'");

        std::cout << "SUCCESS: All 7 core physical processes verified in ProcessRegistry!\n";

        // 3. Test factory creation of process instances
        auto settling_proc = registry.create("settling");
        assert(settling_proc != nullptr && "Factory creation of 'settling' process failed");
        assert(settling_proc->get_name() == "settling");

        auto drydep_proc = registry.create("drydep");
        assert(drydep_proc != nullptr && "Factory creation of 'drydep' process failed");
        assert(drydep_proc->get_name() == "drydep");

        std::cout << "SUCCESS: Process instance factory creation verified!\n";

        // 4. Test C-API process addition to core
        void* core_ptr = catchem_core_create(4, 5, 2);
        assert(core_ptr != nullptr);

        catchem_core_add_process_by_name(core_ptr, "settling");
        catchem_core_add_process_by_name(core_ptr, "drydep");

        catchem_core_destroy(core_ptr);
        std::cout << "SUCCESS: C-API catchem_core_add_process_by_name verified!\n";
    }
    Kokkos::finalize();
    return 0;
}
