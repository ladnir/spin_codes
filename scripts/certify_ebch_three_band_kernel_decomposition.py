#!/usr/bin/env python3
"""Exact kernel-decomposition certificate for the three-band EBCH layout.

Let ``P_b : F_2^64 -> F_2^{n_b}`` be the committed EBCH encoder restricted
to fixed band ``b`` of sizes 42, 43, and 43.  This verifier checks

    ker(P_0) direct-sum ker(P_1) direct-sum ker(P_2) = F_2^64.

The identity is the local algebraic input for the exponent-two finite-field
Brascamp--Lieb inequality used by the packet-profile outer enumerator.
Everything here is exact binary linear algebra.
"""

from __future__ import annotations

import hashlib

from certify_bch_band_projections import BANDS, row_reduce
from probe_bch_forward_xor_circuit import DIMENSION, target_outputs


def kernel_basis(forms: list[int]) -> list[int]:
    reduced, pivots = row_reduce(forms, DIMENSION)
    pivot_set = set(pivots)
    basis = []
    for free in range(DIMENSION):
        if free in pivot_set:
            continue
        vector = 1 << free
        for row, pivot in zip(reduced, pivots):
            if row >> free & 1:
                vector |= 1 << pivot
        if any((form & vector).bit_count() & 1 for form in forms):
            raise SystemExit("three-band kernels: invalid nullspace basis")
        basis.append(vector)
    return basis


def main() -> None:
    outputs = target_outputs()
    kernels = []
    projection_ranks = []
    for begin, end in BANDS:
        forms = outputs[begin:end]
        rank = len(row_reduce(forms, DIMENSION)[1])
        kernel = kernel_basis(forms)
        if rank + len(kernel) != DIMENSION:
            raise SystemExit("three-band kernels: rank-nullity mismatch")
        projection_ranks.append(rank)
        kernels.append(kernel)

    joined = [vector for kernel in kernels for vector in kernel]
    joined_rank = len(row_reduce(joined, DIMENSION)[1])
    if sum(map(len, kernels)) != DIMENSION or joined_rank != DIMENSION:
        raise SystemExit("three-band kernels: kernels do not form a direct sum")

    digest = hashlib.sha256()
    for band, kernel in enumerate(kernels):
        for vector in kernel:
            digest.update(bytes((band,)))
            digest.update(vector.to_bytes(8, "little"))

    print("exact EBCH three-band kernel decomposition")
    print(f"band_sizes={','.join(str(end-begin) for begin, end in BANDS)}")
    print(f"projection_ranks={','.join(map(str, projection_ranks))}")
    print(f"kernel_dimensions={','.join(str(len(kernel)) for kernel in kernels)}")
    print(f"joined_kernel_rank={joined_rank}")
    print(f"kernel_basis_sha256={digest.hexdigest()}")
    print("rank_nullity=PASS")
    print("kernel_direct_sum=PASS")
    print("status=EXACT_INTEGER_THREE_BAND_KERNEL_DECOMPOSITION")


if __name__ == "__main__":
    main()
