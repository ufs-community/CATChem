#pragma once
#include "catchem_state_manager.hpp"
#include <memory>
#include <string>
#include <string_view>

namespace catchem {

    struct ProcessNames {
        static constexpr std::string_view GasChem = "gaschem";
        static constexpr std::string_view Photolysis = "photolysis";
        static constexpr std::string_view Settling = "settling";
        static constexpr std::string_view Dust = "dust";
        static constexpr std::string_view SeaSalt = "seasalt";
        static constexpr std::string_view CarbChem = "carbchem";
        static constexpr std::string_view SO4Chem = "so4chem";
    };

    class ProcessInterface {
    public:
        virtual ~ProcessInterface() = default;
        virtual std::string get_name() const = 0;
        virtual void init(std::shared_ptr<StateManager> state) = 0;
        virtual void run(std::shared_ptr<StateManager> state) = 0;
        virtual void finalize() = 0;
    };

} // namespace catchem
