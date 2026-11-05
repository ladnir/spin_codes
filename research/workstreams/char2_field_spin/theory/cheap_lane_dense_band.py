r"""Directed 35% bound for actual source supports of size 7688 through 8192.

This is a class certificate for the cheap lane-refresh construction, not a
whole-code certificate. It covers arbitrary nonzero GF(2^128) messages in
the class, including messages with full GF256 lane span. It does not count
those messages as though they were independent GF256 messages.

Proof, with q=256, m=65280, P=32768 packets, and N=64P:

Fix an actual active source-group set T, |T|=r. Fix the route, D1, the
state-mask sequence H_0,...,H_(P-1), H_0 empty, and requested output-zero
sets Z_t. For each nonfinal packet set

    M = [A[Z,H]; (H16 D1)[H_next^c,H]],
    N = [I[Z,:]; (H16 D1)[H_next^c,:] C].

Omit the lower rows in the final packet. Preimpose the E-valued equations
L N X=0, where L spans the left kernel of M, on F^(128r). Their kernel U
is defined over E=GF256. If a group projection of U is zero, the proposed
actual support T is impossible. Otherwise, for uniform nonzero u in U,
each group is missing with probability less than 1/q. Markov's inequality
with at most 32 missing groups therefore gives

    P_q(dim U) <= sum_{S subset T, |T\S|<=32} n_S / ((q-1)c),
    c = 1-r/(33q) >= 1/33,

where n_S counts the E-vectors in U with exact group support S. This is an
initial-potential comparison only. It is not a claim about later F-valued
quotients. Fresh lambda pivots contract P_q(d)=(q^d-1)/(q-1) by at most
1/q: independent pivots reduce dimension; dependent pivots retain the
required nonzero form for at most one of m actual coefficients, and m>=q.

Use the positive identity z^wt(Y)=sum_{Z subset zeros(Y)}
z^(N-|Z|)(1-z)^|Z| before averaging the outer. The preserved scalar outer
density beta=q^160/(q-1)^32 applies once to each group in S. It yields
the local scalar auxiliary kernel

    K_S = sum_Z z^(64-|Z|)(1-z)^|Z| q^-rank[M,N_(I(S))].

These auxiliary scalar inputs and states need not retain the actual state
support; this proof imposes no support feasibility on them. Removing one
matrix column lowers rank by at most one. Hence, for j=|T\S|, replacing
S by T costs at most q^(256j) over the complete stream. The total source
factor is bounded by beta^r C32, where

    C32 = sum_{j=0}^32 binom(8192,j) (q^256/beta)^j.

Padding the actual allowed input positions I(T) to all 64 positions costs
q^(N-256r). This holds for every route; no route distribution is changed.

For full input, substitute Y=X+A_H u. Since CA=0, the next pre-refresh
state is u+C Y. Elimination of Y_Z gives the exact identity

    rank[M,N_all] = |Z| + rank R_(H_next^c)(C_(Z^c), I_H),
    R = H16 D1.

Superregularity of H16 and nonzero D1 imply that the last rank is at least
min(16-j,h), where h=|H|, j=|H_next|. Thus K_full is at most
g(z)^64 q^-min(16-j,h), with g(z)=z+(1-z)/q. The final packet, which has
no next-state constraint, contributes g(z)^64 without an extra transition.

After counting destination masks, define

    B[h,j] = binom(16,j) q^-min(16-j,h),
    v[h] = q^-h + q^-8,
    rho = (1+1/q)^16 + 65536 q^-8.

The inequality q^-min(a,h)<=q^-a+q^-h gives Bv<=rho v. Also
v[0]/min_h v[h]=q^8. Consequently the sum over all state-mask paths is at
most g(z)^N q^8 rho^(P-1). There are exactly P-1 charged transitions;
the unused last refresh is not charged or assigned a dummy state mask.

Take L=floor(7N/20)=734003, z=7/3315, and r0=7688. Since beta/q^256<1,
union over all actual supports of size at least r0, using at most 2^8192
sets, gives the explicit upper bound

    U = 2^8192 * 33/(q-1) * q^N * (beta/q^256)^r0 * C32
        * q^8 * rho^(P-1) * g(z)^N / z^L.

All averaging and scalar-density steps occur inside positive sums. No
joint two-shadow density, uniform-F coefficient distribution, or heuristic
state distribution is used. This helper evaluates U with Arb balls and
checks U<2^-98 directly. It has no optimizer and writes no result files.
"""

import argparse
from hashlib import sha256
import json
from math import comb
from pathlib import Path
import shutil


def verify(precision):
    """Replay the fixed class bound with directed Arb arithmetic."""
    if precision not in (192, 256):
        raise ValueError("precision must be 192 or 256 bits")
    source = Path(__file__).resolve()
    free_bytes = shutil.disk_usage(source.parent).free
    if free_bytes < 20 * 1024**3:
        raise RuntimeError("The proof run requires at least 20 GiB free")

    from flint import arb, ctx

    ctx.prec = precision
    q = arb(256)
    groups = 8192
    first_group_count = 7688
    packets = 32768
    records = 64 * packets
    output_limit = 7 * records // 20
    beta = q**160 / arb(255)**32
    missing_factor = sum(
        arb(comb(groups, j)) * (q**256 / beta)**j for j in range(33)
    )
    rho = (1 + 1 / q)**16 + 65536 * q**(-8)
    z = arb(7) / 3315
    g = z + (1 - z) / q
    bound = (
        arb(2)**groups * 33 / (q - 1)
        * q**records * (beta / q**256)**first_group_count
        * missing_factor * q**8 * rho**(packets - 1)
        * g**records / z**output_limit
    )
    passed = bool(bound < arb(2)**(-98))
    if not passed:
        raise ArithmeticError("The directed bound did not prove 98 bits")
    return {
        "schema": "cheap-lane-dense-band-v1",
        "scope": "actual active source-group counts 7688 through 8192",
        "whole_code_certificate": False,
        "full_field_messages": True,
        "directed": True,
        "precision_bits": precision,
        "distance": "7/20",
        "output_weight_limit": output_limit,
        "z": "7/3315",
        "source_groups": groups,
        "minimum_active_source_groups": first_group_count,
        "maximum_missing_shadow_groups": 32,
        "charged_refreshes": packets - 1,
        "target_margin_bits": 98,
        "proved_below_target": passed,
        "log2_upper_interval": str(bound.log() / arb(2).log()),
        "source_sha256": sha256(source.read_bytes()).hexdigest(),
        "free_bytes_at_start": free_bytes,
        "generated_file_bytes": 0,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--precision", type=int, choices=(192, 256), default=192)
    args = parser.parse_args()
    print(json.dumps(verify(args.precision), indent=2))
