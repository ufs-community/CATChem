# Vendored mdspan single header

`include/mdspan/mdspan.hpp` is a generated single-header build of the
[kokkos/mdspan](https://github.com/kokkos/mdspan) reference
implementation of C++23 `std::mdspan`.

It provides the `Kokkos::mdspan` data layer for host-only CATChem
builds (`CATCHEM_ENABLE_KOKKOS=OFF`). It is vendored — not fetched at
configure time — because UFS regression tests compile on
network-restricted compute nodes. With `CATCHEM_ENABLE_KOKKOS=ON` this
copy is not used; Kokkos ships its own bundled mdspan.

## Provenance

- Upstream: https://github.com/kokkos/mdspan
- Branch/commit: `stable` @ `884f17a24301955d47cbb22318f06b8d8bee7ca3`
  (same commit previously pinned by the FetchContent path)
- Generated: 2026-08-06
- License: Apache-2.0 WITH LLVM-exception (header retained in the file)

## How to regenerate / update

From a checkout of kokkos/mdspan at the desired commit:

```bash
git clone https://github.com/kokkos/mdspan && cd mdspan
git checkout <commit>
python3 make_single_header.py include/experimental/mdspan \
  | sed "s|$(pwd)/||g" \
  > <catchem>/src/external/mdspan/include/mdspan/mdspan.hpp
```

The `sed` strips the absolute checkout path that the generator embeds
in its `//BEGIN_FILE_INCLUDE:` comments, so the vendored file diffs
reproducibly. After regenerating, update the commit and date above and
rebuild the `CATCHEM_ENABLE_KOKKOS=OFF` Docker image
(`docker buildx build --platform linux/amd64 -f docker/Dockerfile
--build-arg CATCHEM_ENABLE_KOKKOS=OFF .`) — its ctest run is the gate.

## Local modification policy

None. Never hand-edit `mdspan.hpp`; regenerate from upstream instead.

## Notes for consumers

- The single header defaults `MDSPAN_IMPL_STANDARD_NAMESPACE` to
  `std`; CATChem's `src/core/catchem_kokkos_compat.hpp` defines it to
  `Kokkos` before inclusion so the same `Kokkos::mdspan` spellings
  work in both build flavors.
- C++17 or later is required (CATChem host-only builds use C++17).
