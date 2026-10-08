# Grouped-L4 cut follow-up

This bounded follow-up finds no new obstruction to 40% distance. Its closest
stationary diagnostic is −0.00245038074, improving the previous
−0.00261722887 without changing its sign. The search now includes genuinely
heterogeneous partial groups. It does not prove distance, bound finite
setup-failure probability, or justify a 40-bit claim. The reported 8.33 ms
native control therefore remains an uncertified cost control.

The files added here are `theory/hd_grouped_cut_search.py` and
`theory/hd_grouped_cut_tests.py`. All earlier helpers and native files remain
unchanged. No native build, benchmark, or remote command was run.

## Exact row-set data and the diagnostic objective

Fix R=32768 source rows, b=64 outer coordinates, and N=bR=2097152 output
records. Each source row uses GRS[64,32]. The placement permutes groups of
four rows by the frozen odd strides, preserving the lane inside each group.
All offsets are zero. There is one weighted pre-prefix, an independent
uniform permutation in each region, and the two-additive-state post-inner.
Both scans carry their states across region boundaries.

For a fixed selected source-row set U, write v[g,l] for its indicator in
group g and lane l, with 0<=l<4. Let a=|U| and M=R/4. Define

    W = sum_g (v[g,0]v[g,1] + v[g,1]v[g,2] + v[g,2]v[g,3]),
    E_c = sum_g v[g,3]v[g+a_c,0].

Group indices are modulo M. The exact number of safe pre-prefix inverse
columns in region c is A_c=W+E_c. Only the last region has a correction:

    A_63 = W+E_63 + v[-a_63,3](1-v[0,0]).

This is the final singleton inverse column, not a reset of the prefix.
The selected row-check space has dimension h=32a.

For a union of output intervals, let m be its number of safe positions,
u its number of unsafe positions, and t its number of intervals. The
existing rank-tail cut gives

    G=m-2t,    C=u+2t,    rank(H_S)<=h+C.

If G>h and C<=e-h-1, where e=ceil(2N/5)-1, it produces a nonzero codeword
of weight at most e. This implication is finite and deterministic for
every allowed choice of arithmetic coefficients. It is unchanged here.

The search objective is weaker: replace each fixed-cardinality regional
safe mask by a stationary Bernoulli mask of density p_c=A_c/R. For a common
positive mark eta, let F_eta(p) be the mean optimum of G-eta*C per position
over all interval unions. The floating diagnostic is

    D_eta(U) = (1/64) sum_c F_eta(A_c/R)
               + (2/5)eta - (1+eta)a/(2R).

The code evaluates the same clipped interval-score chain as the earlier
helper. With eta=n/d, its integer difference state obeys

    D_next = w + clip(D,-2(n+d),0),

where w=d at a safe position and w=-n otherwise. The outside optimum
increases by max(0,D). An independent exhaustive test compares the finite
recurrence with every output subset on six positions.

A negative value at even one mark separates the desired gain/cost point
in this stationary model. Mixing interval policies cannot evade that
particular separating mark: the DP already optimizes their marked score,
and randomized mixtures preserve the same upper bound on its mean.
This observation does **not** transfer the stationary calculation into a
finite-probability upper bound for uniform fixed-cardinality masks.

## Why a whole-group bit search misses partial patterns

Only first and last lanes affect E_c. For fixed endpoint bits and a fixed
number of selected interior lanes, the best possible within-group counts
are as follows; these are exact maxima over the four interior choices.

| First, last | Zero interior lanes | One | Two |
|---|---:|---:|---:|
| 0,0 | 0 | 0 | 1 |
| 0,1 or 1,0 | 0 | 1 | 2 |
| 1,1 | 0 | 1 | 3 |

For a nonnegative weighted sum sum_c w_c A_c, put B=sum_c w_c. Starting
from full or empty groups, removing any selected row loses at least B
from an internal group edge. Adding a row to an empty group gains at most
B from cross-group edges; interior lanes gain none. If the exchanged rows
share an edge, the actual gain is no larger. The final singleton preserves
these inequalities. Hence no single-row exchange improves this weighted
objective from a full/empty-group row set.

This is a search barrier, not a theorem that full groups are globally
optimal. The follow-up therefore uses genuinely partial initial patterns,
then permits arbitrary single-row exchanges as well as full-group moves.
Each proposed move uses exact safe-count differences. The pressure tangent
only shortlists candidates; the actual floating pressure must increase
before a move is accepted. There are at most 128 moves per phase, and only
the eight best removal/addition candidates under that tangent are tested.

## Bounded heterogeneous families and results

First, all 65536 sixteen-bit row patterns were repeated over the R rows.
This enumerates every heterogeneous pattern periodic across four L4 groups.
They collapse to 387 distinct exact profiles (a,A_1,A_3,terminal), since
every stride is congruent to 1 or 3 modulo four. Every profile either fails
an exact prerequisite or has a negative stationary mark. The closest
eligible profile has a=8192, A_c=6144 except for one terminal increment,
and diagnostic −0.02255563066 at eta=1/4.

Second, a directed endpoint construction exploits the fact that every
stride is odd: select first lanes in even groups and last lanes in odd
groups. Fill interior lanes inward from those endpoints. For
R/4<=a<=3R/4, all regions then have the exact identity

    A_c = a-R/8.

Thus all cross-group safe edges are present, even though the selected
endpoint sets differ. At source densities 0.3, 0.4, 0.5, 0.6, and 0.7,
the best separating marks in the retained menu give gaps approximately
−0.060833, −0.068588, −0.078921, −0.123003, and −0.204578.

Third, the four-lane Fourier symbol supplies non-periodic partial seeds.
For mode j on M groups, its unnormalized Hermitian matrix has unit entries
on the three internal path edges and cross entry

    z_j = (1/64) sum_c exp(2 pi i j a_c/M).

Thresholding the phases of its leading eigenvector selects different row
sets in the four lanes. The retained modes 87, 619, and 1599 are targeted
near the previous closest scalar-group candidates. All have a=3276.
Each seed receives up to 128 single-row moves, 128 full-group moves, and
128 further single-row moves at eta=1/12.

| Seed mode | Initial partial groups | Final partial groups | Final diagnostic |
|---|---:|---:|---:|
| 87 | 18 | 6 | −0.00245038074 |
| 619 | 122 | 32 | −0.00296251750 |
| 1599 | 208 | 56 | −0.00301828268 |

The best row set has safe counts between 2454 and 3159, total 167495,
and check capacity 104832. Its separating mark remains 1/12 among marks
1/8 through 1/18. For comparison, 128 pressure-directed full-group moves
from the original mode-1599 seed give −0.00254247646 and no partial groups.
These are local search outcomes, not optima over arbitrary source-row sets.

## Validation and the remaining question

Six local tests pass. They exhaust all eight-bit row sets for individual
flip deltas, cardinality-preserving swap deltas, and weighted marginals;
verify the interior table; compare the interval DP with exhaustive subset
selection; verify the full-group single-swap barrier including the terminal;
and check the opposite-parity endpoint identity at seven row counts.
Reproduce the bounded search with

    python -B theory/hd_grouped_cut_search.py --followup --moves 128
    python -B theory/hd_grouped_cut_search.py --moves 128
    python -B theory/hd_grouped_cut_tests.py

The commands are run from this workstream directory. They are mathematical
diagnostics, not encoder benchmarks.

The next useful cut work is finite, rather than a larger Fourier sweep.
Small row sets, fixed-cardinality effects, and row sets chosen after seeing
the regional routes remain unbounded. In particular, a negative stationary
mean does not rule out a cut event whose probability exceeds 2^-40. A finite
failure lower bound for a specified row set and interval policy would
reject the candidate; conversely, these searches supply no upper bound on
that probability. Until one of those tasks is resolved, the appropriate
status is **no new obstruction found; no certificate**.
