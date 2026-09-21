// src/core/catchem_fortran_process.hpp
#pragma once
#include "catchem_process_interface.hpp"
#include <memory>
#include <string>

namespace catchem {

    // C-linkage declarations matching Fortran bridge callbacks
    extern "C" {
    typedef void (*FortranBridgeCallback)(void* state_mgr);
    }

    class FortranProcess : public ProcessInterface {
    private:
        std::string name;
        FortranBridgeCallback bridge_callback;

    public:
        FortranProcess(const std::string& process_name, FortranBridgeCallback callback)
            : name(process_name), bridge_callback(callback) {}

        std::string get_name() const override { return name; }

        void init([[maybe_unused]] std::shared_ptr<StateManager> state) override {
            // Initial setup if required
        }

        void run(std::shared_ptr<StateManager> state) override {
            // 1. Sync device Views to host unified memory
            state->sync_to_host();

            // 2. Invoke the Fortran bridging callback
            if (bridge_callback) {
                bridge_callback(static_cast<void*>(state.get()));
                if (state->chemistry().conc)
                    state->chemistry().conc->mark_host_modified();
            }

            // 3. Sync modified host buffers back to device Views
            state->sync_to_device();
        }

        void finalize() override {
            // Cleanup if required
        }
    };

} // namespace catchem
