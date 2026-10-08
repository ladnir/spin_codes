# Bounded cuts for fixed affine placement

No counterexample was found in this bounded cut audit. The ungrouped
placement passes the retained numerical screens. Permuting common groups
of two or four rows also passes the retained screens, including partial
groups. These are negative diagnostics, not a distance certificate or a
proof that no other cut works. A four-row cost control is therefore a
reasonable next measurement, with its distance status still unresolved.

## Construction and exact safe counts

Fix the 64 odd strides in `theory/hd_two_state_regional_strides.py`.
Let R=32768. The outer code applies GRS[64,32] independently to R message
rows. Coordinate c is placed in the source-row order

    pi_c(j) = stride_c*j mod R.

Every affine offset is zero. The encoder then applies one input-weighted
prefix scan, independent uniform permutations within the 64 coordinate
regions, and the original two-additive-state postscan. Both scans carry
state globally. The pre-scan weights are nonzero. The postscan is the map
in `HD_ADDITIVE_TWO_STATE.md`, including its stated coefficient law.

Fix a set U of a source rows. Its check-coordinate space has dimension
h=32a. One inverse-prefix column has nonzero entries at two successive
positions in the placed stream; the last global column has only one entry.
Call a column safe when its image lies in this check-coordinate space.
For this placement, the two source rows at every nonterminal inverse
column are distinct, so safety is exactly membership of both rows in U.

For each region c, define its cyclic internal-edge count

    E_c(U) = sum_x 1_U(x)*1_U(x+stride_c).

The omitted cyclic edge joins pi_c(R-1)=-stride_c to row zero. At an
inter-region boundary, the actual next source row is also zero. Thus
the cyclic edge and the global-carry edge have the same safety condition.
Only the final global column differs:

    A_c = E_c(U)                       for c<63,
    A_63 = E_63(U) + 1_U(-stride_63)*(1-1_U(0)).

These are exact integer counts, not expected densities. In particular,
all-zero offsets require at most one terminal allowance, not 64 boundary
allowances. After regional routing, the safe mask in region c is a uniform
A_c-subset. The 64 masks are independent for this fixed U.

## Uniform exclusion for small source-row sets

The existing directed spectral receipt proves that every nonconstant
eigenvalue of the symmetric stride average is below 323/1000. Therefore,
for alpha=a/R,

    sum_c A_c <= floor(64R*(alpha^2+(323/1000)*alpha*(1-alpha)))+1.

Exact integer comparisons show that this upper bound is below 32a for
every 1<=a<=8567. This statement holds for every row set of each such size,
not merely for the tested examples. At a=8192, the upper bound is 258081,
against 262144 check coordinates.

The mixed-interval cut requires gain greater than h, and gain cannot
exceed the number of safe columns. Thus this prerequisite cannot hold
for those source-row-set sizes, under any realized regional routing.
The spectral test stops excluding the prerequisite at a=8568. At the other
end, a>=26215 already gives h>ceil(0.4N)-1, so the same rank-plus-cost
witness cannot produce a forbidden support through its stated bound.
Neither observation is a minimum-distance lower bound.

## The retained mixed-interval diagnostic

The inverse postscan has at most two tail dimensions. For a union S of
t output intervals, let m and u be its safe and unsafe position counts.
The rank argument in `HD_TWO_STATE_CUTS.md` gives

    rank(H_S) <= h+u+2t.

Define G=m-2t and C=u+2t. The sufficient failure event is
G>h and C<=ceil(0.4N)-h-2. Unsafe holes are allowed within intervals.
The intervals do not reset encoder states.

For a positive rational mark eta, an interval dynamic program maximizes
G-eta*C. Its clipped score difference has finitely many states. The helper
computes its stationary mean F_eta(p) under Bernoulli safe bits of density p.
It then reports the numerical diagnostic

    gap_eta = average_c F_eta(A_c/R)
              + eta*(0.4-a/(2R)) - a/(2R).

This calculation retains all 64 regional densities. It does not replace
them by their mean. A successful finite cut would exceed the corresponding
score threshold, so negative gaps are useful evidence against the tested
row set. However, a negative stationary expectation is not a finite
probability bound. The calculation does not certify fixed-cardinality
tails, finite boundary losses, or adaptive choices of U.

The ungrouped menu contains 92 row sets:

- 21 cyclic intervals, with densities from 0.275 through 0.775.
- 26 specified unions of low-bit fibers, including parity.
- 40 cosine-threshold sets from eight frequencies, including the five largest positive spectral modes.
- Five deterministic, cardinality-preserving hillclimbs, each limited to 512 swaps and initialized at mode1082.

The hillclimb maximizes the number of internal stride edges, not the
interval score. These families are bounded probes, not an exhaustive
search over low-bit unions or row sets.

Of the 92 sets, 39 fail the safe-count prerequisite exactly and one fails
the check-capacity prerequisite. All remaining stationary diagnostics are
negative. The closest is mode15931 with a=13107: it has 460580 safe columns,
419424 checks, and gap approximately -0.06132919447 at eta=1/3.
No row set in this menu justified an expensive finite failure calculation.

## Grouping two or four rows for store locality

For group size L in {2,4}, the proposed placement is

    pi_c(j) = L*(stride_c*floor(j/L) mod (R/L)) + (j mod L).

It preserves the order inside each common group while affinely permuting
the groups. This changes the construction. The ungrouped safe-count
exclusion above must not be transferred unchanged to it.

The safe counts also have a simple exact form for arbitrary partial groups.
Let v_(b,l) indicate whether row Lb+l belongs to U. Define

    W = sum_b sum_(l=0)^(L-2) v_(b,l)*v_(b,l+1),
    first_b = v_(b,0),       last_b = v_(b,L-1).

Then the cyclic count in region c is

    W + sum_b last_b*first_(b+stride_c).

The same single terminal correction applies in the final region.
For U consisting of m whole groups, this becomes

    (L-1)*m + internal_edges_c(selected groups),

again with at most one final correction. The earlier deterministic ceiling
for a scalar bidiagonal postscan does not apply to this dense-inverse,
two-state postscan merely because the placement preserves groups.

The retained grouped menus contain 69 row sets for L=2 and 129 for L=4.
Each menu includes cyclic whole-group sets, top-five Fourier thresholds
on the quotient row space, and every nonempty proper within-group pattern
at five selected-group densities. The last family explicitly tests partial
groups; it does not cover every assignment of different patterns to groups.

For L=2, ten sets fail the safe-count prerequisite and all other diagnostics
are negative. For L=4, 51 sets fail that prerequisite. Three coarse mark
diagnostics were initially positive; testing the retained mark eta=1/12
makes all three negative:

| L=4 quotient mode | Source rows a | Refined gap at eta=1/12 |
|---|---:|---:|
| 1599 | 3276 | -0.002617228871 |
| 87 | 3276 | -0.002685068991 |
| 619 | 3276 | -0.002683936511 |

Every retained L=4 set therefore has a negative diagnostic at one of the
tested marks. This permits a cost-only prototype without a presently known
cut counterexample. It does not establish a 40-bit failure guarantee.

## Artifacts, checks, and next step

`theory/hd_affine_cut.py` contains the exact safe counts, bounded menus,
and numerical pressure screen. Its separate `exact_profile_moments` path
can evaluate the fixed 3G-C block policy using exact without-replacement
moments for heterogeneous regional counts. That path was not used to
claim a finite failure in this audit.

All nine tests pass. They check literal inverse-prefix supports on every
small bit string, both ordinary and grouped placements, the terminal
correction, the exact spectral cardinality range, and the grouped identity
for whole groups. They also cross-check the stationary chain and retain
the three refined grouped diagnostics as explicitly numerical regressions.
No native build, benchmark, remote command, or existing-file edit was made.

The next useful measurement is the L=4 cost control. The next proof task
is a finite tail bound for grouped cuts and a source-rank argument for
the fixed placement. General check subspaces, route-adaptive row sets,
and heterogeneous partial-group patterns remain open; negative screens
must not be presented as a distance certificate.

```text
python -B research/workstreams/char2_field_spin/theory/hd_affine_cut.py
python -B research/workstreams/char2_field_spin/theory/hd_affine_cut.py --group 2
python -B research/workstreams/char2_field_spin/theory/hd_affine_cut.py --group 4
python -B -m unittest discover -s research/workstreams/char2_field_spin/theory -p hd_affine_cut_tests.py -v
```
