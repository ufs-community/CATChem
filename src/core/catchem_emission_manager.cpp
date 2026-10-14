#include "catchem_emission_manager.hpp"

namespace catchem {

    void EmissionManager::add_category(const EmissionCategory& category) {
        categories[category.category_name] = category;
        is_loaded = true;
    }

    bool EmissionManager::has_category(const std::string& name) const {
        return categories.find(name) != categories.end();
    }

    const EmissionCategory* EmissionManager::get_category(const std::string& name) const {
        auto it = categories.find(name);
        return it != categories.end() ? &it->second : nullptr;
    }

} // namespace catchem
