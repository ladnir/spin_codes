# Case 7: bounded searches and an exact Boolean reduction

Follow-up: `KERNEL_COSET_PAYOFF.md` tests the suggested coset constraints.
Their gain is below 0.005 bits, with a checked ceiling below 0.0057 bits.
This auxiliary route is now parked in favor of the main BCH/estimator work.

Case 7 remains open. This turn found neither a rank-31 witness nor a
weight-30 counterexample. The full P-dual distance lower bound is still 30,
and the unconditional K=2^28 BCH/RM2Sub margin is still 42.510028 bits.
The additional 1.075004 bits remain conditional on excluding weight 30.

Two independently checked weight-36 words give an explicit upper bound
d(P dual)<=36. These are words of the auxiliary dual code, not counterexamples
to the BCH outer or to the SPIN distance certificate.

## New structured description of the remaining case

Use P dual, field coordinates, and F7 from `PDUAL_RANK_PROGRESS.md`. Regard a
binary word as a Boolean function on the eight coordinates of F256 in the
polynomial basis for modulus 0x14D. Its algebraic normal form is the unique
multilinear polynomial representing that function over F2.

An exact Boolean transform on every generator proves the following:

- Every word in P dual has Boolean degree at most five.
- Its degree-five part is an injective linear function of the eight bits
  of F7. These degree-five parts form an eight-dimensional space.
- Consequently, the subcode E=ker(F7) from the previous turn is exactly
  the intersection of P dual with RM(4,8). It has dimension 117.

Here RM(4,8) means all Boolean functions of degree at most four on eight
variables. E is a specific subspace, not the whole Reed-Muller code.
The checker eliminates the paired values (F7, degree-five part) on a basis
of P dual and checks that the resulting map has rank eight in both images.

Coordinate scaling by a nonzero field element a multiplies F7 by a^7 and
preserves P dual and weight. Since gcd(7,255)=1, a -> a^7 permutes F256's
nonzero elements. Any potential weight-30 word in case 7 can therefore be
normalized to F7=1. The earlier kernel theorem has already excluded F7=0.

Fix any word g in P dual with F7(g)=1; surjectivity guarantees one exists.
The remaining proof question is precisely whether the affine coset g+E
contains a word of weight 30. The normalized degree-five part is fixed,
not an additional unknown. Its 25 monomials are recorded by their eight-bit
variable masks in `generated/pdual_quintic_structure_v1.json`.
This reduction does not establish a distance bound for g+E.

## A stronger final rank step was tested

Under F7=1, all values F(7*2^j) equal one. Thus a sum of two corresponding
support vectors has coordinate sum zero. The old rank witnesses used only
individual monomials with zero Fourier value and individual nonzero pivots.

The new test takes a checked rank-29 prefix. After a valid transformation,
it replaces the next append v7 by v7+vz for an exponent z in the zero set.
That appended vector still has coordinate sum one and is independent of
the preceding zero-sum span. A final transformation could send its two terms
to equal known unit values, while sending the other 29 rows into the zero
set. Appending v7 again would then prove rank at least 31.

The implementation checks each of these premises explicitly. It tested
10,112 rank-30 candidates from a bounded beam and found no such extension.
This rules out neither other prefixes nor more general polynomial witnesses.
No rank improvement is claimed from this run.

## Counterexample searches

| Search | Bounded work | Outcome |
|---|---|---|
| SAT with exact weight 30 and F7=1 | 45 seconds | UNKNOWN: timeout |
| One-/two-row information-set search | 533,159 randomized bases; 45 seconds | Best weight 36 |
| Four-row meet-in-the-middle search | 225,151 randomized bases; 45 seconds | Best weight 36 |

The SAT model imposes orthogonality to a generator basis of P0, weight 30,
and eight binary equations for F7=1 on the 255 punctured coordinates. The
normalization is safe: the existing affine argument lets a zero coordinate
be moved to the extension position, and scaling preserves that position.
No UNSAT certificate was produced or assumed.

The information-set searches operate in the full extended P dual. Each
randomized systematic basis is obtained solely by binary row operations.
The second search also combines pairs from disjoint halves of the basis
when their projections onto ten nonpivot coordinates agree. The hot kernels
use fixed four-limb XOR and population-count operations and fixed-size arrays.

The two searches evaluated approximately 4.20 billion and 2.59 billion
candidates, respectively; 812 million of the latter were four-row candidates.
These counts include repetitions. Candidates are not uniform independent
samples of the low shell. Neither count supplies a confidence interval or
a quantitative upper bound on A30.

Both reported words passed independent checks of all 131 orthogonality
equations, their exact weights, and their F7 values. The respective F7 values
are 159 and 225, so both lie in the residual case rather than E. The compact
receipts contain their full hexadecimal words. Thus the checked interval
from these results is 30<=d(P dual)<=36, not a proof of distance 32 or 36.

## Where the literature check stops

The classical Kasami--Tokura--Azumi classification cited here covers weights
below 2.5 times the minimum distance. For RM(5,8), that threshold is 20,
so this particular result does not classify our weight-30 target.
See [the original 1976 paper](https://www.sciencedirect.com/science/article/pii/S0019995876903557).
We did not substitute a theorem for the larger range or assume that the
fixed degree-five part alone forces weight at least 32.

## Next bounded proof effort

Do not spend the next turn merely increasing these search budgets. The
useful new object is the fixed affine coset g+E. Derive a small set of spectrum
or coset constraints from E, then test whether they exclude the existing bad
LP spectrum before investing in a larger solve. Preserve the restriction to
the actual E; arbitrary degree-four perturbations give a different problem.

If a rank approach is retained, it needs a genuinely more general witness
than the single two-term final append tested here. Any claimed improvement
must still close all residual case-7 branches. The completed case-15 cover
and the distance-32 kernel theorem remain intact.

## Artifacts and verification

The new sources are `pdual_normalized_extension.py`, `pdual_weight30_sat.py`,
`pdual_low_weight.cpp`/`.py`, `pdual_stern_search.cpp`/`.py`, and
`pdual_quintic_structure.py`. Compact receipts use the corresponding `_v1.json`
names under `generated/`. Executables and bulk search state are not proof
inputs; the reported words are checked directly.

```text
python -B pdual_low_weight.py --output generated/pdual_low_weight_v1.json --verify
python -B pdual_quintic_structure.py --output generated/pdual_quintic_structure_v1.json --verify
python -B -m unittest test_pdual_case7_search test_pdual_refinement test_low_shell_roi test_curve_spectrum
```

No existing certificate source or other workstream was modified. The CWC
presentation keeps exact reductions and explicit witnesses separate from
unsuccessful searches and unproved distance claims.
