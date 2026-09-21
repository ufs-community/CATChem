# C++20 Kokkos::mdspan Shared State Integration Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Integrate the standard-conforming C++20 Kokkos backport of `mdspan` directly into `catchem::InteropField`, enabling parenthesis-based multidimensional span indexing over active shared memory states.

**Architecture:** We will include the Kokkos-bundled header `<mdspan/mdspan.hpp>` in `catchem_interop_field.hpp`. We will add template helpers inside `InteropField` to map Rank (1, 2, 3) to the corresponding dynamic `Kokkos::mdspan` types using `Kokkos::layout_left` layout. We will define an `.mdspan()` method that constructs and returns this mdspan over the active Kokkos view buffer.

**Tech Stack:** C++20, Kokkos v5.1.1, CMake

## Global Constraints

- Target C++20 utilizing the Kokkos backport of mdspan (`Kokkos::mdspan`), avoiding direct dependency on C++23 `<mdspan>`.
- Retain Fortran column-major storage layout (`Kokkos::layout_left` / `Kokkos::LayoutLeft`) across the pointer boundary to achieve zero-copy execution on CPU targets.

---

### Task 1: Integrate `Kokkos::mdspan` into `InteropField`

**Files:**
- Modify: `src/core/catchem_interop_field.hpp`

**Interfaces:**
- Produces: `catchem::InteropField::MdspanType` type alias
- Produces: `catchem::InteropField::mdspan()` const method

- [ ] **Step 1: Write the failing compilation test**

Create a temporary test to verify that `mdspan()` is declared:
```cpp
// tests/test_mdspan_compilation.cpp
#include "catchem_interop_field.hpp"
int main() {
    double data[6] = {0.0};
    catchem::InteropField<double, 2> field(data, {2, 3});
    auto m = field.mdspan();
    return 0;
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `g++ -std=c++20 -Isrc/core -I/usr/local/include -c tests/test_mdspan_compilation.cpp`
Expected: FAIL (no member `mdspan`)

- [ ] **Step 3: Modify `catchem_interop_field.hpp`**

Insert the dynamic namespace macros, standard headers, SFINAE template helpers, and the `.mdspan()` getter into `src/core/catchem_interop_field.hpp`:

```cpp
#pragma once
#include <Kokkos_Core.hpp>
#include <mdspan/mdspan.hpp>
#include <vector>
#include <memory>
#include <type_traits>

namespace catchem {

template <typename DataType, int Rank>
class InteropField {
public:
    using HostSpace = Kokkos::HostSpace;
    using DeviceSpace = Kokkos::DefaultExecutionSpace::memory_space;

    template <typename T, int R, typename Space, bool Unmanaged>
    struct ViewType;

    template <typename T, typename Space, bool Unmanaged>
    struct ViewType<T, 1, Space, Unmanaged> {
        using type = typename std::conditional_t<Unmanaged,
            Kokkos::View<T*, Kokkos::LayoutLeft, Space, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,
            Kokkos::View<T*, Kokkos::LayoutLeft, Space>>;
    };

    template <typename T, typename Space, bool Unmanaged>
    struct ViewType<T, 2, Space, Unmanaged> {
        using type = typename std::conditional_t<Unmanaged,
            Kokkos::View<T**, Kokkos::LayoutLeft, Space, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,
            Kokkos::View<T**, Kokkos::LayoutLeft, Space>>;
    };

    template <typename T, typename Space, bool Unmanaged>
    struct ViewType<T, 3, Space, Unmanaged> {
        using type = std::conditional_t<Unmanaged,
            Kokkos::View<T***, Kokkos::LayoutLeft, Space, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,
            Kokkos::View<T***, Kokkos::LayoutLeft, Space>>;
    };

    using HostViewType = typename ViewType<DataType, Rank, HostSpace, true>::type;
    using DeviceViewType = typename ViewType<DataType, Rank, DeviceSpace, false>::type;

    // --- SFINAE type-bound helpers for standard Kokkos::mdspan mapping ---
    template <int R>
    struct MdspanTypeHelper;

    template <>
    struct MdspanTypeHelper<1> {
        using type = Kokkos::mdspan<DataType, Kokkos::extents<int, Kokkos::dynamic_extent>, Kokkos::layout_left>;
    };

    template <>
    struct MdspanTypeHelper<2> {
        using type = Kokkos::mdspan<DataType, Kokkos::extents<int, Kokkos::dynamic_extent, Kokkos::dynamic_extent>, Kokkos::layout_left>;
    };

    template <>
    struct MdspanTypeHelper<3> {
        using type = Kokkos::mdspan<DataType, Kokkos::extents<int, Kokkos::dynamic_extent, Kokkos::dynamic_extent, Kokkos::dynamic_extent>, Kokkos::layout_left>;
    };

    using MdspanType = typename MdspanTypeHelper<Rank>::type;

    HostViewType host_view;
    DeviceViewType device_view;
    bool is_gpu_target;

    InteropField(DataType* ptr, const std::vector<int>& dims) {
        is_gpu_target = !std::is_same_v<HostSpace, DeviceSpace>;

        if constexpr (Rank == 1) {
            host_view = HostViewType(ptr, dims[0]);
            if (is_gpu_target) device_view = DeviceViewType("dev_field_1d", dims[0]);
        } else if constexpr (Rank == 2) {
            host_view = HostViewType(ptr, dims[0], dims[1]);
            if (is_gpu_target) device_view = DeviceViewType("dev_field_2d", dims[0], dims[1]);
        } else if constexpr (Rank == 3) {
            host_view = HostViewType(ptr, dims[0], dims[1], dims[2]);
            if (is_gpu_target) device_view = DeviceViewType("dev_field_3d", dims[0], dims[1], dims[2]);
        }
    }

    void sync_to_device() {
        if constexpr (!std::is_same_v<HostSpace, DeviceSpace>) {
            Kokkos::deep_copy(device_view, host_view);
        }
    }

    void sync_to_host() {
        if constexpr (!std::is_same_v<HostSpace, DeviceSpace>) {
            Kokkos::deep_copy(host_view, device_view);
        }
    }

    auto view() const {
        if constexpr (!std::is_same_v<HostSpace, DeviceSpace>) {
            return device_view;
        } else {
            return host_view;
        }
    }

    // --- standard C++20 non-owning mdspan mapping ---
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
};

} // namespace catchem
```

- [ ] **Step 4: Verify compilation inside Docker**

Run: `docker run --rm -v $(pwd):/workspace -w /workspace/build-test cece-dev:latest bash -c "apt-get update -y && apt-get install -y python3 && cmake .. -DENABLE_KOKKOS=ON -DENABLE_TESTING=ON && make test_catchem_interop"`
Expected: PASS

- [ ] **Step 5: Clean up and commit**

```bash
rm -f tests/test_mdspan_compilation.cpp
git add src/core/catchem_interop_field.hpp
git commit -m "feat(core): introduce type-safe C++20 Kokkos::mdspan accessors inside InteropField"
```

---

### Task 2: Validate `mdspan` Indexing in Integration Suite

**Files:**
- Modify: `tests/test_catchem_interop.cpp`

**Interfaces:**
- Consumes: `catchem::InteropField::mdspan()`

- [ ] **Step 1: Write integration assertions**

Add **TEST 8: Standard C++20 mdspan Representation and Indexing** inside `tests/test_catchem_interop.cpp` right before `return 0;`:
```cpp
        // ==========================================
        // TEST 8: Standard C++20 mdspan Representation and Indexing
        // ==========================================
        {
            int n_cols = 4;
            int n_levels = 5;
            int n_species = 2;

            void* core = catchem_core_create(n_cols, n_levels, n_species);
            auto* state_obj = static_cast<catchem::StateManager*>(catchem_core_get_state_manager(core));

            // 1. Bind mock temperature 3D array (using layout left column-major)
            std::vector<double> temp_array(n_cols * n_levels, 298.15);
            temp_array[0 + 0 * n_cols] = 273.15; // Bottom-left level 0
            temp_array[1 + 2 * n_cols] = 300.00; // Col 1, Level 2

            catchem_state_bind_met_3d(state_obj, "T", temp_array.data());
            catchem_state_sync_to_device(state_obj);

            // 2. Extract standard mdspan accessor
            auto temp_mds = state_obj->met.T->mdspan();

            // 3. Assert dimensions and layout access
            assert(temp_mds.extent(0) == n_cols);
            assert(temp_mds.extent(1) == n_levels);
            assert(temp_mds(0, 0, 0) == 273.15);
            assert(temp_mds(1, 2, 0) == 300.00);

            std::cout << "INFO: mdspan dimension 0 (cols) = " << temp_mds.extent(0) << "\n";
            std::cout << "INFO: mdspan dimension 1 (levels) = " << temp_mds.extent(1) << "\n";
            std::cout << "INFO: Verified mdspan(0,0,0) = " << temp_mds(0, 0, 0) << " K\n";
            std::cout << "INFO: Verified mdspan(1,2,0) = " << temp_mds(1, 2, 0) << " K\n";

            catchem_core_destroy(core);
            std::cout << "SUCCESS: C++20 Kokkos::mdspan Multidimensional Access Validation Passed!\n";
        }
```

- [ ] **Step 2: Build and execute full suite inside Docker**

Run: `docker run --rm -v $(pwd):/workspace -w /workspace/build-test cece-dev:latest bash -c "apt-get update -y && apt-get install -y python3 && cmake .. -DENABLE_KOKKOS=ON -DENABLE_TESTING=ON && make test_catchem_interop && ./tests/test_catchem_interop"`

Expected: PASS, printing:
```text
SUCCESS: C++20 Kokkos::mdspan Multidimensional Access Validation Passed!
```

- [ ] **Step 3: Commit**

```bash
git add tests/test_catchem_interop.cpp
git commit -m "test(interop): add standard mdspan indexing and layout assertions"
```
