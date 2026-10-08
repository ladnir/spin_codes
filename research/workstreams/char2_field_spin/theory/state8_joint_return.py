"""Joint selected-output-zero/return rank cover for symmetric state8 maps.

The full numerical equation matrix is exposed once. Its return rows and
selected output rows are not charged as independent events. This is a
local upper cover, not a whole-code certificate.
"""
from fractions import Fraction
from functools import lru_cache
from math import comb

from state_rank_cap_transfer import return_flat_factor
from state_rank_cap_transfer import return_minimum_weight as minimum16
from state8_width32_reference import return_minimum_weight as minimum32


def minimum_return_weight(weight, width):
    if width == 16:
        return minimum16(weight)
    if width == 32:
        return minimum32(weight)
    raise ValueError("only the checked widths16 and32 are supported")


@lru_cache(maxsize=48)
def joint_return_terms(weight, width=16):
    if not 1 <= weight <= width:
        raise ValueError("joint-return terms require a nonempty input")
    minimum = minimum_return_weight(weight, width)
    # v is the number of UNSELECTED output coordinates, not exact weight.
    return tuple((v, comb(width, v), max(8, 1+weight-v))
                 for v in range(minimum, width+1))


@lru_cache(maxsize=48)
def split_joint_return_terms(weight, width=16):
    """(active zeros i, inactive zeros z, multiplicity, joint rank charge).

    Inactive selected zeros satisfy z<h<=8. Their equations and active
    identity pivots have rank i+z. The invertible full map also gives rank
    >=h+w-v>=w-width+i+2z+1. All bounds concern the SAME joint matrix.
    """
    if not 1 <= weight <= width:
        raise ValueError("joint-return terms require a nonempty input")
    minimum = minimum_return_weight(weight, width)
    return tuple((i, z, comb(weight, i)*comb(width-weight, z),
                  max(8, i+z, weight-width+i+2*z+1))
                 for i in range(weight+1)
                 for z in range(min(7, width-weight)+1)
                 if width-i-z >= minimum)


def return_bounds(weight, y, tau, width=16, cap=8):
    """Return (pure, split-joint) bounds, including the weighted state flats.

    cap bounds the incoming state-form rank but requires actual nonzero
    witness state. For 1<=tau<=255, arbitrary common-history source-form
    dependencies are permitted. All current diagonals must remain fresh.
    """
    if width not in (16, 32) or not 0 <= weight <= width:
        raise ValueError("invalid checked width or input weight")
    if not (0 < y <= 1 and 1 <= tau <= 255 and 1 <= cap <= 8):
        raise ValueError("invalid output mark, rank tilt, or state-rank cap")
    y = Fraction(y) if isinstance(y, int) else y
    tau = Fraction(tau) if isinstance(tau, int) else tau
    if not weight:
        return y*0, y*0
    flats = return_flat_factor(weight, tau, cap)
    pure = flats*y**minimum_return_weight(weight, width)/tau**8
    joint = sum(sum(comb(8, k)/tau**k
                    for k in range(min(weight-1, cap-1, 7-z)+1))
                *count*y**(width-i-z)*(1-y)**(i+z)/tau**charge
                for i, z, count, charge in split_joint_return_terms(weight, width))
    return pure, joint


def joint_return_bound(weight, y, tau, width=16, cap=8):
    """Minimum of two COMPLETE valid covers, never their product."""
    return min(return_bounds(weight, y, tau, width, cap))
