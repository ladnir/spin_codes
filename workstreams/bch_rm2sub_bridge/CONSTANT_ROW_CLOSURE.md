# Closing occupancies by separating the all-one row

Update: the complete K=2^16 instance is now certified at 53.9443672720 bits.
See [T128_S19_M16_CLOSURE.md](T128_S19_M16_CLOSURE.md) for all-occupancy
coverage, the exact ledger, and reproduction. The anchor results below record
the earlier partial checkpoint and explain the split used in the full proof.

For the selected t128_s19 map at K=2^16, the all-one-row split closes the
previously weak full-occupancy endpoint. It also closes the first sparse
gap, Q=150, and the intermediate anchors Q=200,256,384. These results use
the same fixed BCH-256 outer and inner map as the implementation study.
They do not change the encoder.

Each new certificate includes every possible number of all-one rows at
its stated occupancy. The remaining work is coverage between occupancies,
not missing all-one-row cases within these anchors.

## The bound

Let L=512 be the number of outer rows, and fix the number Q of nonzero
rows. Write h for the number of all-one rows and d=Q-h for the remaining
nonzero rows. An all-one row supplies exactly one bit to every region;
it needs neither a Bernoulli majorant nor a choice among weight bands.

The ordinary nonzero rows use the 12 nonconstant bands already present in
the sparse certificate engine. For each band g, fix a Bernoulli probability
p_g and its exact row-counting factor Gamma_g. The BCH spectrum upper caps
determine Gamma_g, so no unknown spectrum is replaced by a heuristic.

For a fixed positive tilt lambda, let R_j be the outward four-state
transfer for a region containing j input ones. The four classes retain
zero, arbitrary nonzero, uniform nonzero, and bounded-density nonzero
states. The exact character bound from `syndrome_density_v1.py` controls
the state created at activation. Its inequalities are implemented in
`activation_density_arb.py` with Arb arithmetic.

Define a matrix sequence V_0(j)=R_j. Each ordinary row applies the
entrywise envelope

\[
V_{r+1}(j)=\max_g\Gamma_g^{1/256}
 \bigl((1-p_g)V_r(j)+p_gV_r(j+1)\bigr).
\]

For the case with d ordinary rows and h all-one rows, use V_d(h).
This matrix dominates every assignment of the d ordinary rows to bands.
The row-permutation counting inequality and the 256 independent region
permutations give the following bound on that case's expected bad-message
count:

\[
\binom LQ\binom Qh\,12^d\,e^{\lambda H}
 e_Z^\mathsf T V_d(h)^{256}\mathbf1,
\qquad H=13107.
\]

Here e_Z selects the zero initial state. The state continues between
regions; the matrix product does not reset it. The binomial factors count
the occupied row positions and which occupied rows are all-one.

Sum these bounds over all h=0,...,Q. Different fixed witnesses may be
used for different cases. Taking the smallest certified bound for each
case is valid before the final sum.

The previous adaptive bound included the all-one word among its 13 band
choices. Fixing h removes that deterministic word from the adaptive
maximum and from the ordinary-band label count. This change, rather than
the activation-density refinement alone, produced the endpoint closure.

## Verified results

The producer uses 256-bit Arb arithmetic and the existing directed scaled
recurrence. An independent 512-bit run replays each fixed witness. All
case sums below are exact rational sums of retained dyadic bounds.

| Occupancy Q | All-one-row cases included | Retained margin for this occupancy |
|---:|---:|---:|
| 150 | 151 | 72.7616 bits |
| 200 | 201 | 72.3489 bits |
| 256 | 257 | 71.9944 bits |
| 384 | 385 | 71.4113 bits |
| 512 | 513 | 70.9972 bits |

The evaluator retains at most 80 bits per case. These margins are therefore
conservative stored bounds, not estimates of the actual margin or its curve.

Together with the sparse ledger, the certificates cover Q=1,...,150,
Q=162,...,170, and Q=200,256,384,512. The total covered contribution
corresponds to 48.7642 bits. The decrease from the earlier 54-bit partial
ledger comes from adding conservatively bounded sparse occupancies, not
from a changed code or a measured loss. This is not yet a full-code margin:
all other occupancies remain unclosed.

`summarize_constant_coverage.py` authenticates the sparse and constant-split
replay receipts and recomputes the combined exact union. Its output is
`generated/t128_s19_m16_combined_constant_coverage_v1.json`.

The successful endpoint receipt is
`generated/t128_s19_m16_q512_constant_v2.json`; the first-gap receipt is
`generated/t128_s19_m16_q150_constant_v2.json`. The three intermediate
anchors use `v1`. Each has a matching `_replay.json` receipt.

The bridge regression run passed 18 tests, including exact small checks of
every band sequence and constant-row offset in the new directed fold.
The seven inherited activation-density tests also passed. They exercise
exact syndrome counts and exhaustive small-field state transitions.

## Other approaches tested

The resumed sparse search first extended coverage to Q=1,...,149 and
Q=162,...,170. The activation-density refinement alone improved the old
endpoint bound only slightly, leaving its margin near -8824 bits.

A complete three-category type cover stalled near -7262 bits. Its bulk
row majorant at p=1/2 pays about 149.22 bits per row because of the upper
caps at weights 38 and 218. Combining all nontrivial weights therefore
loses too much information. Separating those tails into five categories
improved the complete binary64 cover from roughly -94804 to -26229 bits,
but did not close it within the bounded search. Those covers are diagnostics,
not additional certificates.

The new `certify_density_types.py` provides an outward checker for such
covers. It recomputes exact counting costs and normalizes witness probabilities
as rationals. No full passing type cover has been produced with that checker.

## Next step

The all-one-row split has now been automated across the complete remaining
range, replayed, and merged with the sparse ledger. The next target is
K=2^18, followed by K=2^20. Existing frozen producers and numerical receipts
have not been modified.
