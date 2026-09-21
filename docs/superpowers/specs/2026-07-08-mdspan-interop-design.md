# Spec: C++20 Kokkos::mdspan Shared State Integration

## Overview
This specification details the design for incorporating standard-conforming C++20 multidimensional spans (`Kokkos::mdspan`) directly into CATChem's memory interop layer (`InteropField`). This enables standard-conforming, multi-dimensional array slicing and parenthesis-based indexing inside C++ physics solvers and host model interfaces, while fully preserving underlying column-major alignments (`Kokkos::layout_left`) and GPU/CPU synchronizations.

## 1. Objectives & Benefits
* **Standard-Conforming Semantics:** Provides an idiomatic C++20 standard multidimensional indexing interface (`std::mdspan` model) directly for all meteorological and chemical state fields.
* **Unified Memory Access:** Allows developers to write generic, non-Kokkos C++ physics kernels using standard `mdspan` while retaining native zero-copy alignments across language boundaries.
* **Seamless CPU/GPU Targets:** Integrates transparently with `InteropField` automatic dual-space host/device mirrors, binding standard `mdspan` over whichever memory buffer is active on the physical hardware target.

## 2. Foundational Memory Layer (`InteropField`)
We will include the Kokkos-bundled reference `mdspan` header `<mdspan/mdspan.hpp>` inside `src/core/catchem_interop_field.hpp`.

To support compile-time type-safety, we will declare type-safe aliasing helpers inside `InteropField` using template specializations for Rank 1, 2, and 3:

```cpp
template <typename DataType, int Rank>
class InteropField {
public:
    // ... Existing View space allocations ...

    // SFINAE / template helpers to map Kokkos::mdspan dynamically per Rank
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
```

We will implement the `.mdspan()` accessor to return the `mdspan` mapped directly over the raw data pointer of the active Kokkos View:

```cpp
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
```

## 3. Interoperability and Legacy Support
* **Dual-Space Coherence:** Any updates performed via the `.mdspan()` view will modify the underlying View's active execution memory directly. Subsequent synchronizations (`sync_to_host()`, `sync_to_device()`) continue to function without any changes, maintaining absolute bitwise consistency.
* **Direct Integration with MetState/ChemState:** Since `MetState` and `ChemState` hold structured `InteropField` pointers, every field (e.g., `met.T`, `chem.conc`) immediately inherits `.mdspan()` support.
