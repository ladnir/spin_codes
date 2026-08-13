#!/usr/bin/env python3
"""Exact rational case ledger for composing fixed-band collision components.

The local envelopes are proved by:

* ``certify_three_band_pair_core.py``;
* ``certify_three_band_composition_local.py``; and
* ``certify_three_band_degree3_transfer_outward.py``.

This file performs only the final positive generating-function arithmetic.
The uniform r=4..64 conditional continuation premise is certified by
``certify_three_band_degree4_64_transfer_outward.py``; every number below is
an exact rational consequence of the named local certificates.
"""

from __future__ import annotations

import math
from fractions import Fraction

from outward_log2 import log2_fraction
from certify_three_band_exact_length import (
    graph_factor_upper,
    worst_puncture_block_adjustment,
)


PHYSICAL_BLOCKS = 1 << 14
MAXIMUM_COMPONENT_BLOCKS = 323


def report(name: str, value: Fraction) -> None:
    print(f"{name}_log2_upper={log2_fraction(value).hi}")


def main() -> None:
    pole = Fraction(1, 20)
    cutoff_conversion = Fraction(181, 50) ** 106
    graph_factor = graph_factor_upper(pole)
    residual_block = Fraction(1, 1 << 54)
    residual_all = (1 + residual_block) ** PHYSICAL_BLOCKS
    residual_nonempty = residual_all - 1
    residual_at_least_two = (
        residual_nonempty - PHYSICAL_BLOCKS * residual_block
    )
    residual_at_least_three = (
        residual_at_least_two
        - math.comb(PHYSICAL_BLOCKS, 2) * residual_block**2
    )

    decorated_root = (
        Fraction(1, 1 << 125)
        * 768
        * cutoff_conversion
        * graph_factor
    )
    decorated_root_degree_at_least_five = (
        Fraction(1, 1 << 165)
        * 768
        * cutoff_conversion
        * graph_factor
    )
    decorated_child = Fraction(1, 1 << 112)
    decorated_children_nonempty = decorated_child / (1 - decorated_child)

    no_high_high_two_batches = (
        decorated_root * decorated_children_nonempty * residual_all
    )
    no_high_high_one_batch_three_residual = (
        decorated_root * residual_at_least_three
    )
    no_high_high_degree_five_two_residual = (
        decorated_root_degree_at_least_five * residual_at_least_two
    )
    no_high_high_degree_five_one_residual = Fraction(1, 1 << 50)

    adjacent_base = Fraction(3, 4 * (1 << 39))
    degree_three_transfer = Fraction(11, 16 * (1 << 30))
    conditional_high_child = (
        degree_three_transfer
        * MAXIMUM_COMPONENT_BLOCKS
        * 378
        * math.comb(63, 2)
        * worst_puncture_block_adjustment(pole) ** 3
    )
    # The companion outward sweep proves every r=4..64 row is smaller.
    extension_activity = conditional_high_child + residual_nonempty
    if extension_activity >= 1:
        raise SystemExit("composition ledger: extension activity does not contract")
    adjacent_continuation = (
        adjacent_base * extension_activity / (1 - extension_activity)
    )

    case_bounds = {
        "no_high_high_two_or_more_batches": no_high_high_two_batches,
        "no_high_high_one_batch_three_or_more_residual": (
            no_high_high_one_batch_three_residual
        ),
        "no_high_high_degree_ge_five_two_or_more_residual": (
            no_high_high_degree_five_two_residual
        ),
        "no_high_high_degree_ge_five_one_residual": (
            no_high_high_degree_five_one_residual
        ),
        "adjacent_high_batches_nonempty_continuation": adjacent_continuation,
    }
    for name, value in case_bounds.items():
        report(name, value)
    high_connected = sum(case_bounds.values(), Fraction(0))
    report("high_connected_case_sum", high_connected)
    if high_connected > Fraction(1, 1 << 41):
        raise SystemExit("composition ledger: high connected sum exceeds 2^-41")
    print("high_connected_case_sum_le_2^-41=PASS")
    print("r_ge_4_conditional_transfer_envelope=OUTWARD_VERIFIED")
    print("status=EXACT_RATIONAL_COMPOSITION_LEDGER")


if __name__ == "__main__":
    main()
