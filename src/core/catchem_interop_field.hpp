#pragma once
#include "catchem_field_contract.hpp"
#include "catchem_kokkos_compat.hpp"
#include <array>
#include <memory>
#include <stdexcept>
#include <type_traits>
#include <vector>
#ifdef CATCHEM_ENABLE_KOKKOS
#include <mdspan/mdspan.hpp>
#endif

namespace catchem {

    // --- SFINAE template helpers for standard Kokkos::mdspan mapping defined at namespace scope ---
    template <typename DataType, int Rank> struct MdspanTypeHelper;

    template <typename DataType> struct MdspanTypeHelper<DataType, 1> {
        using type = Kokkos::mdspan<DataType, Kokkos::extents<int, Kokkos::dynamic_extent>, Kokkos::layout_left>;
    };

    template <typename DataType> struct MdspanTypeHelper<DataType, 2> {
        using type = Kokkos::mdspan<DataType, Kokkos::extents<int, Kokkos::dynamic_extent, Kokkos::dynamic_extent>,
                                    Kokkos::layout_left>;
    };

    template <typename DataType> struct MdspanTypeHelper<DataType, 3> {
        using type =
            Kokkos::mdspan<DataType,
                           Kokkos::extents<int, Kokkos::dynamic_extent, Kokkos::dynamic_extent, Kokkos::dynamic_extent>,
                           Kokkos::layout_left>;
    };

    /**
     * @class InteropField
     * @brief High-performance wrapper mapping host-allocated contiguous pointers to device-capable views.
     *
     * This class translates standard raw pointer buffers across language boundaries (e.g. Fortran to C++).
     * With Kokkos enabled it provides unmanaged LayoutLeft Views plus a managed device mirror; in host-only
     * builds the same interface is served directly by Kokkos::mdspan (standalone library) over the host
     * buffer, and the device-sync operations are no-ops.
     *
     * @tparam DataType Numeric type of the grid elements (typically double).
     * @tparam Rank Dimensional dimensionality of the field (1, 2, or 3).
     */
    template <typename DataType, int Rank> class InteropField {
    public:
        using MdspanType = typename MdspanTypeHelper<DataType, Rank>::type;
        std::array<std::size_t, Rank> immutable_extents{};
        FieldContract contract;
        std::size_t generation = 0;
        AvailabilityState availability = AvailabilityState::Current;
        LatestWriter latest_writer = LatestWriter::HostCurrent;
        std::size_t host_to_device_sync_count = 0;
        std::size_t device_to_host_sync_count = 0;

        void mark_host_modified() { latest_writer = LatestWriter::HostCurrent; }
        void mark_device_modified() { latest_writer = LatestWriter::DeviceCurrent; }
        void set_generation(std::size_t value) {
            generation = value;
            availability = AvailabilityState::Current;
        }
        void invalidate() { availability = AvailabilityState::Unavailable; }
        bool is_current(std::size_t value) const {
            return availability == AvailabilityState::Current && generation == value;
        }
        std::size_t extent(int axis) const { return immutable_extents.at(static_cast<std::size_t>(axis)); }

        const DataType* host_read() {
            sync_to_host();
            return host_data();
        }

        DataType* host_write() {
            sync_to_host();
            mark_host_modified();
            return host_data();
        }

        auto device_read() {
            sync_to_device();
            return view();
        }

        auto device_write() {
            sync_to_device();
            mark_device_modified();
            return view();
        }

    private:
        void validate_dimensions(const std::vector<int>& dims) {
            if (dims.size() != static_cast<std::size_t>(Rank))
                throw std::invalid_argument("InteropField dimension rank mismatch");
            for (int r = 0; r < Rank; ++r) {
                if (dims[r] <= 0)
                    throw std::invalid_argument("InteropField dimensions must be positive");
                immutable_extents[r] = static_cast<std::size_t>(dims[r]);
            }
            contract.extents.assign(immutable_extents.begin(), immutable_extents.end());
        }

    public:
#ifdef CATCHEM_ENABLE_KOKKOS
        using HostSpace = Kokkos::HostSpace;
        using DeviceSpace = Kokkos::DefaultExecutionSpace::memory_space;

        template <typename T, int R, typename Space, bool Unmanaged> struct ViewType;

        template <typename T, typename Space, bool Unmanaged> struct ViewType<T, 1, Space, Unmanaged> {
            using type = typename std::conditional_t<
                Unmanaged, Kokkos::View<T*, Kokkos::LayoutLeft, Space, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,
                Kokkos::View<T*, Kokkos::LayoutLeft, Space>>;
        };

        template <typename T, typename Space, bool Unmanaged> struct ViewType<T, 2, Space, Unmanaged> {
            using type = typename std::conditional_t<
                Unmanaged, Kokkos::View<T**, Kokkos::LayoutLeft, Space, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,
                Kokkos::View<T**, Kokkos::LayoutLeft, Space>>;
        };

        template <typename T, typename Space, bool Unmanaged> struct ViewType<T, 3, Space, Unmanaged> {
            using type = std::conditional_t<
                Unmanaged, Kokkos::View<T***, Kokkos::LayoutLeft, Space, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,
                Kokkos::View<T***, Kokkos::LayoutLeft, Space>>;
        };

        using HostViewType = typename ViewType<DataType, Rank, HostSpace, true>::type;
        using DeviceViewType = typename ViewType<DataType, Rank, DeviceSpace, false>::type;

        HostViewType host_view;     ///< Contiguous host layout left view mapping raw pointer.
        DeviceViewType device_view; ///< Dedicated device view managed by Kokkos execution space.
        bool is_gpu_target;         ///< Flag indicating if device-host data sync is required.

        /**
         * @brief Constructs InteropField and maps raw host memory.
         * @param ptr Raw host pointer (must be non-null).
         * @param dims Vector indicating grid bounds.
         * @throws std::invalid_argument if ptr is null or dimensions are invalid.
         */
        InteropField(DataType* ptr, const std::vector<int>& dims) {
            if (ptr == nullptr) {
                throw std::invalid_argument("InteropField constructor failed: input pointer is null.");
            }
            validate_dimensions(dims);

            is_gpu_target = !std::is_same_v<HostSpace, DeviceSpace>;

            if constexpr (Rank == 1) {
                host_view = HostViewType(ptr, dims[0]);
                if (is_gpu_target)
                    device_view = DeviceViewType("dev_field_1d", dims[0]);
            } else if constexpr (Rank == 2) {
                host_view = HostViewType(ptr, dims[0], dims[1]);
                if (is_gpu_target)
                    device_view = DeviceViewType("dev_field_2d", dims[0], dims[1]);
            } else if constexpr (Rank == 3) {
                host_view = HostViewType(ptr, dims[0], dims[1], dims[2]);
                if (is_gpu_target)
                    device_view = DeviceViewType("dev_field_3d", dims[0], dims[1], dims[2]);
            }
        }

        /** @brief Copy host data buffer to the selected Kokkos device view. */
        void sync_to_device() {
            if (latest_writer != LatestWriter::HostCurrent)
                return;
            if constexpr (!std::is_same_v<HostSpace, DeviceSpace>) {
                Kokkos::deep_copy(device_view, host_view);
            }
            ++host_to_device_sync_count;
            latest_writer = LatestWriter::Synchronized;
        }

        /** @brief Updates the host raw pointer and updates the unmanaged host view. */
        void update_host_pointer(DataType* ptr) {
            if (ptr == nullptr) {
                throw std::invalid_argument("InteropField::update_host_pointer failed: input pointer is null.");
            }
            mark_host_modified();
            if constexpr (Rank == 1) {
                host_view = HostViewType(ptr, host_view.extent(0));
            } else if constexpr (Rank == 2) {
                host_view = HostViewType(ptr, host_view.extent(0), host_view.extent(1));
            } else if constexpr (Rank == 3) {
                host_view = HostViewType(ptr, host_view.extent(0), host_view.extent(1), host_view.extent(2));
            }
        }

        /** @brief Sync calculated outputs from Kokkos device back to host buffer. */
        void sync_to_host() {
            if (latest_writer != LatestWriter::DeviceCurrent)
                return;
            if constexpr (!std::is_same_v<HostSpace, DeviceSpace>) {
                Kokkos::deep_copy(host_view, device_view);
            }
            ++device_to_host_sync_count;
            latest_writer = LatestWriter::Synchronized;
        }

        /** @brief Retrieves the active Kokkos view mapped to current execution space. */
        auto view() const {
            if constexpr (!std::is_same_v<HostSpace, DeviceSpace>) {
                return device_view;
            } else {
                return host_view;
            }
        }

        /** @brief Raw pointer to the host-side buffer. */
        DataType* host_data() const { return host_view.data(); }

        /** @brief Constructs a modern, standard-conforming mdspan referencing active data memory. */
        MdspanType mdspan() const {
            auto v = view();
            if constexpr (Rank == 1) {
                return MdspanType(v.data(), v.extent(0));
            } else if constexpr (Rank == 2) {
                return MdspanType(v.data(), v.extent(0), v.extent(1));
            } else if constexpr (Rank == 3) {
                return MdspanType(v.data(), v.extent(0), v.extent(1), v.extent(2));
            }
        }
#else
        DataType* data_ptr;         ///< Raw host buffer (not owned).
        std::array<int, Rank> dims; ///< Grid bounds per rank.
        bool is_gpu_target = false; ///< Host-only builds never target a device.

        InteropField(DataType* ptr, const std::vector<int>& dim_vec) : data_ptr(ptr) {
            if (ptr == nullptr) {
                throw std::invalid_argument("InteropField constructor failed: input pointer is null.");
            }
            validate_dimensions(dim_vec);
            for (int r = 0; r < Rank; ++r)
                dims[r] = dim_vec[r];
        }

        void sync_to_device() {
            if (latest_writer == LatestWriter::HostCurrent) {
                ++host_to_device_sync_count;
                latest_writer = LatestWriter::Synchronized;
            }
        }
        void sync_to_host() {
            if (latest_writer == LatestWriter::DeviceCurrent) {
                ++device_to_host_sync_count;
                latest_writer = LatestWriter::Synchronized;
            }
        }

        /** @brief Updates the host raw pointer. */
        void update_host_pointer(DataType* ptr) {
            if (ptr == nullptr) {
                throw std::invalid_argument("InteropField::update_host_pointer failed: input pointer is null.");
            }
            data_ptr = ptr;
            mark_host_modified();
        }

        /** @brief Raw pointer to the host-side buffer. */
        DataType* host_data() const { return data_ptr; }

        /** @brief Host mdspan over the buffer (the only "view" in host-only builds). */
        MdspanType mdspan() const {
            if constexpr (Rank == 1) {
                return MdspanType(data_ptr, dims[0]);
            } else if constexpr (Rank == 2) {
                return MdspanType(data_ptr, dims[0], dims[1]);
            } else if constexpr (Rank == 3) {
                return MdspanType(data_ptr, dims[0], dims[1], dims[2]);
            }
        }

        /** @brief Same as mdspan(); kernels index views and mdspans identically via operator(). */
        MdspanType view() const { return mdspan(); }
#endif
    };

} // namespace catchem
