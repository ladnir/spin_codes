"""Local positive output/feedback rank tags; not a global certificate.

Fix allowed input positions I, incoming/outgoing state masks H,H', and
requested output zeros Z. Stack the output and feedback equations as
[M,N_I]. Let k=rank(M), d=rank([M,N_I])-k, h=|H|, and let kf be the
feedback M rank. A supplied tag df must lower-bound the feedback input
rank. Put p=|Z intersect I|, e=rank(A[Z minus I,H]), and
ea=rank(A[Z,H]). Conservative LOWER bounds for e,ea are also valid.

Appending rows cannot decrease the input rank difference: if the M rank
increases, the joint rank increases too. Consequently d>=df. Eliminating
the occupied output identities gives output joint rank p+e. Thus
k>=k0=max(kf,ea), k+d>=R0=max(kf+df,p+e), and k<=h.

For tau>1, the weighted point mass tau^-k delta_d has CDF at most
  F(t)=0                         if t<df or max(k0,R0-t)>h;
       tau^-max(k0,R0-t)          otherwise.
F is nondecreasing. Its positive differences therefore dominate the
actual point mass against every nonnegative NONINCREASING continuation
in the rank debit. This follows by writing the continuation as a positive
sum of prefix indicators. It covers the best-32-column rank potential,
its exact full-rank cutoff, and fixed-threshold positive markers x<=1.
One must not interpret these tags as probabilities or normalize them.

F decreases pointwise as df increases. A stochastically lower feedback
rank measure, or separate upper bounds for its rank-class masses, can
therefore be averaged through F. Output ranks do not depend on fresh D1.
Activation feasibility must remain in any activation tag measure; this
helper's numerical screen considers only LIVE incoming states.

The output-zero expansion has positive weights
z^(64-|Z|)(1-z)^|Z|. The screen sums those weights using conservative
cardinality ranks of the actual four RS[16,4] blocks. It does not assume
that the coefficient-basis emission matrix is superregular. At h=16
these cardinality ranks are exact. Other state-block profiles are local
diagnostics, not a maximum over all masks of a given total state size.
The old floor-only expression is reported ONLY for comparison. It does
not bound the same raw output-subset CDF: the latter can retain extinct
state vectors. Neither multiplying the old/new expressions nor taking
their local minimum is justified here. Their separately completed
global proofs could be compared, but no global proof is evaluated here.

--exclude-impossible drops a zero subset when an active RS block has
at least four requested zero positions outside I. These positions have
no input contribution. Any four RS rows have full column rank, so they
would force every state coordinate in that block to zero, contradicting
the specified exact active state. Every genuine zero subset survives
this test. The remaining positive weights are NOT renormalized.

--gamma splits a fixed global mark z=a*b, a=z^gamma, b=z^(1-gamma).
For every genuine output, z^wt <= a^floor*b^wt. Expand only b^wt and
multiply its rank-tag measure by a^floor. A global proof would keep the
prefactor z^-L, not b^-L. Gamma0 gives the new bound; gamma1 recovers the
old bound exactly. The screen reports gamma0,.5,1 separately when asked;
it never takes their pointwise rank-CDF minimum. A later operator may
fix gamma per local I,H,H' before the D1 and zero-subset averages.
The global z must be fixed before the output-zero expansion. New rank
tags depend on Z, so independently retuning z by the later rank-threshold
class is NOT justified by the older feedback-only threshold argument.

--self-test checks small synthetic full matrices exactly over GF2, with
rational tau, plus arbitrary monotone continuations and a best-32 cutoff.
--screen is FLOAT diagnostic only, not directed arithmetic. It checks
local caps, not route/profile averages, source counts, or boundary terms.
No files or caches are written. Each invocation requires20 GiB free and
has a50-second budget. There is no global search or encoder benchmark.
"""

import argparse
from collections import defaultdict
from fractions import Fraction
import itertools
import json
import math
from pathlib import Path
import random
import shutil
import sys
import time

sys.dont_write_bytecode = True
DEADLINE = None
TAU = 65280


def guard():
    if DEADLINE is not None and time.monotonic() > DEADLINE:
        raise RuntimeError("bounded50-second invocation exhausted")


def rank_cdf(h, kf, df, p, e, ea, ell, tau):
    """Positive weighted rank-CDF envelope; tau may be an exact Fraction."""
    if not (0 <= kf <= h and 0 <= df <= ell and 0 <= p <= ell
            and 0 <= e <= ea <= h and tau > 1):
        raise ValueError("invalid rank envelope arguments")
    k0, r0 = max(kf, ea), max(kf+df, p+e)
    result = []
    for t in range(ell+1):
        required = max(k0, r0-t)
        result.append(0 if t < df or required > h else tau**(-required))
    return result


def differences(cdf):
    """Exact inputs give exact positive differences; FLOAT is diagnostic."""
    previous, tags = 0, []
    for current in cdf:
        assert current >= previous-1e-14*max(1, abs(previous))
        tags.append(max(0, current-previous))
        previous = current
    return tags


def rank2(rows):
    """Rank of bit-packed rows over GF2, for exact synthetic tests only."""
    pivots = {}
    for row in rows:
        while row:
            pivot = row.bit_length()-1
            if pivot in pivots:
                row ^= pivots[pivot]
            else:
                pivots[pivot] = row
                break
    return len(pivots)


def exact_tests():
    rng = random.Random(31057)
    cases = inequalities = continuations = 0
    for sample in range(120):
        guard()
        h, nout, nfb = rng.randrange(4), rng.randrange(1, 6), rng.randrange(4)
        ell = rng.randrange(min(3, nout)+1)
        inputs = sorted(rng.sample(range(nout), ell))
        index = {point: i for i, point in enumerate(inputs)}
        emission = [rng.randrange(1 << h) for _ in range(nout)]
        mf = [rng.randrange(1 << h) for _ in range(nfb)]
        nf = [rng.randrange(1 << ell) for _ in range(nfb)]
        jf = [a | b << h for a, b in zip(mf, nf)]
        kf, df_actual = rank2(mf), rank2(jf)-rank2(mf)
        for zero_mask in range(1 << nout):
            zeros = [i for i in range(nout) if zero_mask >> i & 1]
            mout = [emission[i] for i in zeros]
            nout_rows = [(1 << index[i]) if i in index else 0 for i in zeros]
            joint_out = [a | b << h for a, b in zip(mout, nout_rows)]
            p = sum(i in index for i in zeros)
            e = rank2([emission[i] for i in zeros if i not in index])
            ea = rank2(mout)
            k = rank2(mf+mout)
            d = rank2(jf+joint_out)-k
            assert d >= df_actual
            assert rank2(joint_out) == p+e
            assert k >= max(kf, ea) and k+d >= max(kf+df_actual, p+e)
            cases += 1
            for tau in (Fraction(2), Fraction(3), Fraction(17)):
                exact_weight = tau**(-k)
                previous_df = None
                for df in range(df_actual+1):
                    cdf = rank_cdf(h, kf, df, p, e, ea, ell, tau)
                    assert all(a <= b for a, b in zip(cdf, cdf[1:]))
                    assert all(exact_weight*(d <= t) <= cdf[t]
                               for t in range(ell+1))
                    if previous_df is not None:
                        assert all(a <= b for a, b in zip(cdf, previous_df))
                    previous_df = cdf
                    tags = differences(cdf)
                    tests = [[Fraction(t <= cutoff) for t in range(ell+1)]
                             for cutoff in range(ell+1)]
                    tests += [[tau**(-t) for t in range(ell+1)],
                              [Fraction(ell+1-t) for t in range(ell+1)]]
                    for function in tests:
                        assert sum(a*b for a, b in zip(tags, function)) >= exact_weight*function[d]
                        continuations += 1
                    inequalities += ell+1
    # A genuine best32 continuation remains nonincreasing even at its cutoff.
    fixed = [i % 5 for i in range(63)]
    dim, tau = 124, Fraction(3)
    values = []
    for debit in range(33):
        rank = sum(sorted(fixed+[debit], reverse=True)[:32])
        values.append((tau**(dim-rank)-1)/(tau-1) if rank < dim else 0)
    assert all(a >= b for a, b in zip(values, values[1:]))
    assert values[0] > 0 and values[-1] == 0
    # Mixtures with an actual source-rank tag are valid without normalization.
    for ell in range(1, 5):
        for df in range(ell+1):
            cdf = rank_cdf(2, 1, df, min(ell, 2), 1, 2, ell, Fraction(3))
            assert all(a <= b for a, b in zip(cdf, cdf[1:]))
    return {"status": "exact tests passed", "synthetic_zero_set_cases": cases,
            "rational_cdf_inequalities": inequalities,
            "monotone_continuation_checks": continuations,
            "best32_cutoff_monotonicity": "passed",
            "scope": "local algebra and order properties, not a distance certificate"}


def feedback_tags(occ, h, j):
    """Existing conservative D1 rank-class upper masses; never normalize."""
    c, m = sum(min(4, x) for x in occ), 16-j
    k = max(0, c-h)
    if j == 0:
        return [(k, 1.0)]
    if h >= m or not k:
        return [(0, 1.0)]
    checks, slots = m-h, 16-h
    result = []
    for d in range(max(0, checks+k-slots), min(checks, k)+1):
        cap = (1.0 if d == k else
               min(1.0, ((256**k-1)//255)*255.0**(-checks)
                   /((256**(k-d)-1)//255)))
        result.append((d, cap))
    return result


def output_histogram(occ, state_blocks, z, exclude_impossible=False):
    """Positive zero expansion, grouped by (p,e_lower,ea_lower).

    For v requested zeros, rank(A_block[zeros,H_b]) is at least
    max(0,h_b-max(0,4-v)). It is exact when h_b is0 or4. This cap is
    uniform over point sets and within-block state masks, including the
    paired route's shared points; it does not require independent masks.
    """
    histogram = {(0, 0, 0): 1.0}
    for ell, h in zip(occ, state_blocks):
        guard()
        block = defaultdict(float)
        for p in range(ell+1):
            for v in range(17-ell):
                if exclude_impossible and h and v >= 4:
                    continue
                e = max(0, h-max(0, 4-v))
                ea = max(0, h-max(0, 4-p-v))
                weight = (math.comb(ell, p)*math.comb(16-ell, v)
                          * z**(16-p-v)*(1-z)**(p+v))
                if weight:
                    block[p, e, ea] += weight
        following = defaultdict(float)
        for (p, e, ea), left in histogram.items():
            for (q, f, fa), right in block.items():
                following[p+q, e+f, ea+fa] += left*right
        histogram = following
        if len(histogram) > 65536:
            raise RuntimeError("local output histogram exceeds bounded capacity")
    total = sum(histogram.values())
    assert 0 < total <= 1+2e-12
    if not exclude_impossible:
        assert math.isclose(total, 1.0, rel_tol=2e-12)
    return histogram


def local_tags(occ, state_blocks, j, z, histogram, floor_tilt=1.0):
    """Sum positive tags directly, not differences of averaged FLOAT CDFs."""
    h, ell = sum(state_blocks), sum(occ)
    kf = min(16-j, h)
    tags = feedback_tags(occ, h, j)
    active_blocks = (h+3)//4
    floor = max(0, 13*active_blocks-sum(sorted(occ, reverse=True)[:active_blocks]))
    old = [0.0]*(ell+1)
    for df, cap in tags:
        old[df] += cap*TAU**(-kf)*z**floor
    new = [0.0]*(ell+1)
    for (p, e, ea), weight in histogram.items():
        for df, cap in tags:
            cdf = rank_cdf(h, kf, df, p, e, ea, ell, float(TAU))
            for t, entry in enumerate(differences(cdf)):
                new[t] += weight*cap*entry
    new = [entry*floor_tilt**floor for entry in new]
    return old, new


def moment(tags, rho):
    return sum(weight*rho**d for d, weight in enumerate(tags))


def screen(zs, incoming, exclude_impossible=False, gammas=(0.0,)):
    patterns = {16: (4, 4, 4, 4), 8: (4, 4, 0, 0), 4: (4, 0, 0, 0)}
    occupations = [(1, 1, 0, 0), (2, 2, 0, 0), (2, 2, 1, 1), (8, 8, 8, 8)]
    reports = []
    for occ, h, z, gamma in itertools.product(occupations, incoming, zs, gammas):
        guard()
        state_blocks = patterns[h]
        a, b = z**gamma, z**(1-gamma)
        histogram = output_histogram(occ, state_blocks, b, exclude_impossible)
        entries = []
        for j in (8, 12, 15, 16):
            guard()
            if h+j < 17-4*sum(x != 0 for x in occ):
                continue
            old, new = local_tags(occ, state_blocks, j, z, histogram, a)
            if gamma == 1:
                assert all(math.isclose(x, y, rel_tol=2e-12, abs_tol=0)
                           for x, y in zip(old, new))
            moments = []
            for name, rho in (("unweighted", 1.0),
                              ("half_rank", 1/math.sqrt(TAU)),
                              ("full_rank", 1/TAU)):
                old_moment, new_moment = [moment(tags, rho) for tags in (old, new)]
                moments.append({"continuation": name,
                                "old_log2": math.log2(old_moment),
                                "new_log2": math.log2(new_moment),
                                "raw_gain_bits": math.log2(old_moment/new_moment)})
            entries.append({"outgoing": j, "moments": moments,
                            "full_rank_tag_mass": new[-1],
                            "deficient_tag_mass": sum(new[:-1])})
        reports.append({"occupancies": occ, "incoming_state_blocks": state_blocks,
                        "z": z, "gamma": gamma, "a": a, "b": b,
                        "positive_output_features": len(histogram),
                        "retained_zero_subset_weight": sum(histogram.values()),
                        "entries": entries})
    return {"status": "FLOAT local diagnostic; not a band proof",
            "reports": reports, "whole_code_certificate": False,
            "old_new_comparison": "different proof weights; no local minimum or product",
            "exclude_impossible_zero_sets": exclude_impossible,
            "state_profile_scope": "h16 is unique; h8/h4 are specified clustered profiles only",
            "rank_scope": "positive source-rank tags preserved; no global composition evaluated"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument("--screen", action="store_true")
    parser.add_argument("--exclude-impossible", action="store_true")
    parser.add_argument("--gamma", type=float, nargs="+", default=[0.0])
    parser.add_argument("--z", type=float, nargs="+", default=[.99, .995, .998])
    parser.add_argument("--incoming", type=int, nargs="+", default=[16])
    args = parser.parse_args()
    if not args.self_test and not args.screen:
        parser.error("choose --self-test and/or --screen")
    if not (1 <= len(args.z) <= 3 and all(.1 <= z < 1 for z in args.z)
            and 1 <= len(args.incoming) <= 3
            and all(h in (4, 8, 16) for h in args.incoming)
            and 1 <= len(args.gamma) <= 3
            and all(g in (0.0, .5, 1.0) for g in args.gamma)):
        parser.error("use at most3 marks in [.1,1), incoming sizes4,8,16, gamma0,.5,1")
    free = shutil.disk_usage(Path(__file__).resolve().anchor).free
    if free < 20*(1 << 30):
        parser.error("at least20 GiB free required")
    global DEADLINE
    start = time.monotonic()
    DEADLINE = start+50
    result = {"whole_code_certificate": False, "free_bytes_at_start": free}
    if args.self_test:
        result["self_test"] = exact_tests()
    if args.screen:
        result["screen"] = screen(args.z, args.incoming, args.exclude_impossible,
                                  args.gamma)
    result["elapsed_seconds"] = time.monotonic()-start
    result["generated_data_bytes"] = 0
    print(json.dumps(result, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
