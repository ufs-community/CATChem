#pragma once
#include "catchem_error.hpp"
#include <condition_variable>
#include <cstdint>
#include <memory>
#include <mutex>
#include <unordered_map>
#include <utility>
#include <vector>

namespace catchem {
    enum class HandleType : std::uint32_t { Core = 1, State = 2, Diagnostic = 3, Configuration = 4, Time = 5 };
    enum class HandleOwnership : std::uint8_t { Owned, Borrowed };
    struct HandleRecord {
        HandleType type = HandleType::Core;
        std::uint32_t contract_version = 1;
        std::uint64_t generation = 0;
        const void* owner = nullptr;
        HandleOwnership ownership = HandleOwnership::Borrowed;
    };

    class AdmissionLease {
        struct Control {
            explicit Control(HandleRecord value) : record(std::move(value)) {}
            HandleRecord record;
            std::mutex mutex;
            std::condition_variable drained;
            std::size_t active = 0;
            bool closing = false;
        };
        friend class HandleRegistry;
        std::vector<std::shared_ptr<Control>> controls_;
        explicit AdmissionLease(std::vector<std::shared_ptr<Control>> value) : controls_(std::move(value)) {}

    public:
        AdmissionLease() = default;
        AdmissionLease(const AdmissionLease&) = delete;
        AdmissionLease& operator=(const AdmissionLease&) = delete;
        AdmissionLease(AdmissionLease&&) noexcept = default;
        AdmissionLease& operator=(AdmissionLease&& other) noexcept {
            if (this != &other) {
                release();
                controls_ = std::move(other.controls_);
            }
            return *this;
        }
        ~AdmissionLease() { release(); }
        explicit operator bool() const noexcept { return !controls_.empty(); }
        void release() noexcept {
            for (auto iterator = controls_.rbegin(); iterator != controls_.rend(); ++iterator) {
                const auto& control = *iterator;
                std::lock_guard<std::mutex> lock(control->mutex);
                if (control->active > 0)
                    --control->active;
                if (control->closing && control->active == 0)
                    control->drained.notify_all();
            }
            controls_.clear();
        }
    };

    class HandleRegistry {
        using Control = AdmissionLease::Control;

    public:
        static HandleRegistry& instance() {
            static HandleRegistry registry;
            return registry;
        }
        bool add(const void* handle, HandleRecord record) {
            if (!handle)
                return false;
            std::lock_guard<std::mutex> lock(mutex_);
            return records_.emplace(handle, std::make_shared<Control>(std::move(record))).second;
        }
        std::pair<BoundaryStatus, AdmissionLease> acquire(const void* handle, HandleType expected,
                                                          const void* owner = nullptr,
                                                          std::uint64_t generation = 0) noexcept {
            if (!handle)
                return {BoundaryStatus::NullArgument, AdmissionLease{}};
            std::vector<std::shared_ptr<Control>> controls;
            {
                std::lock_guard<std::mutex> registry_lock(mutex_);
                const auto found = records_.find(handle);
                if (found == records_.end())
                    return {BoundaryStatus::InvalidHandle, AdmissionLease{}};
                const auto& record = found->second->record;
                if (record.type != expected)
                    return {BoundaryStatus::WrongHandleType, AdmissionLease{}};
                if (owner && record.owner != owner)
                    return {BoundaryStatus::InvalidHandle, AdmissionLease{}};
                if (generation && record.generation != generation)
                    return {BoundaryStatus::StaleGeneration, AdmissionLease{}};
                if (record.owner) {
                    const auto parent = records_.find(record.owner);
                    if (parent == records_.end())
                        return {BoundaryStatus::InvalidHandle, AdmissionLease{}};
                    controls.push_back(parent->second);
                }
                controls.push_back(found->second);
                for (const auto& control : controls)
                    control->mutex.lock();
                bool closing = false;
                for (const auto& control : controls)
                    closing = closing || control->closing;
                if (!closing)
                    for (const auto& control : controls)
                        ++control->active;
                for (auto iterator = controls.rbegin(); iterator != controls.rend(); ++iterator)
                    (*iterator)->mutex.unlock();
                if (closing)
                    return {BoundaryStatus::InvalidHandle, AdmissionLease{}};
            }
            return {BoundaryStatus::Success, AdmissionLease(std::move(controls))};
        }
        BoundaryStatus validate(const void* handle, HandleType expected, const void* owner = nullptr,
                                std::uint64_t generation = 0) const noexcept {
            auto acquired = const_cast<HandleRegistry*>(this)->acquire(handle, expected, owner, generation);
            return acquired.first;
        }
        BoundaryStatus close_and_wait(const void* handle, HandleType expected) noexcept {
            std::shared_ptr<Control> control;
            {
                std::lock_guard<std::mutex> registry_lock(mutex_);
                const auto found = records_.find(handle);
                if (found == records_.end())
                    return handle ? BoundaryStatus::InvalidHandle : BoundaryStatus::NullArgument;
                if (found->second->record.type != expected)
                    return BoundaryStatus::WrongHandleType;
                control = found->second;
                {
                    std::lock_guard<std::mutex> lock(control->mutex);
                    control->closing = true;
                }
                records_.erase(found);
            }
            std::unique_lock<std::mutex> lock(control->mutex);
            control->drained.wait(lock, [&] { return control->active == 0; });
            return BoundaryStatus::Success;
        }
        void remove(const void* handle) noexcept {
            std::lock_guard<std::mutex> lock(mutex_);
            records_.erase(handle);
        }
        void invalidate_children(const void* owner) noexcept {
            std::lock_guard<std::mutex> lock(mutex_);
            for (auto iterator = records_.begin(); iterator != records_.end();) {
                if (iterator->second->record.owner == owner)
                    iterator = records_.erase(iterator);
                else
                    ++iterator;
            }
        }
        std::size_t size() const noexcept {
            std::lock_guard<std::mutex> lock(mutex_);
            return records_.size();
        }

    private:
        mutable std::mutex mutex_;
        std::unordered_map<const void*, std::shared_ptr<Control>> records_;
    };
} // namespace catchem
