/**
 * @file catchem_kokkos_compat.hpp
 * @brief Single include point for the Kokkos / host-only build split.
 *
 * With CATCHEM_ENABLE_KOKKOS defined this pulls in Kokkos proper. Without
 * it, the standalone kokkos/mdspan library provides the Kokkos::mdspan data
 * types (same header, same namespace), and the Kokkos function-annotation
 * macros are masked with plain `inline` so shared kernel code compiles
 * host-only unchanged.
 */
#pragma once

// mdspan's multi-argument operator[] needs C++23; enable the reference
// implementation's operator() so kernels can index views and mdspans
// identically (Kokkos::View uses operator() natively). This must be set
// before ANY <mdspan/mdspan.hpp> include, in both the Kokkos and host-only
// flavors: under -std=gnu++23 the standalone header sees std::mdspan and
// drops its operator(), breaking `view(i, j)` indexing across the core.
#ifndef MDSPAN_USE_PAREN_OPERATOR
#define MDSPAN_USE_PAREN_OPERATOR 1
#endif

#ifdef CATCHEM_ENABLE_KOKKOS

#include <Kokkos_Core.hpp>

#else

// The vendored single-header mdspan (src/external/mdspan) defaults its
// namespace to std; pin it to Kokkos so both build flavors use the same
// Kokkos::mdspan spellings.
#ifndef MDSPAN_IMPL_STANDARD_NAMESPACE
#define MDSPAN_IMPL_STANDARD_NAMESPACE Kokkos
#endif
#include <mdspan/mdspan.hpp>

#ifndef KOKKOS_INLINE_FUNCTION
#define KOKKOS_INLINE_FUNCTION inline
#endif
#ifndef KOKKOS_FUNCTION
#define KOKKOS_FUNCTION inline
#endif
#ifndef KOKKOS_LAMBDA
#define KOKKOS_LAMBDA [=]
#endif

// Lifecycle no-ops so hosts and tests can call the Kokkos runtime
// entry points unconditionally.
namespace Kokkos {
    inline void initialize(int, char*[]) {}
    inline void initialize() {}
    inline void finalize() {}
} // namespace Kokkos

#endif
