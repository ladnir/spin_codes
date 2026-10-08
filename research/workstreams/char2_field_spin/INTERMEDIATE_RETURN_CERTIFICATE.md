# A complete intermediate-weight class from return constraints

The [32,22] candidate now has a finite bound for every bad word whose
inner input has weight at least 0.19N. This combines the new interval
0.19N <= W < 0.24N with the earlier dense-input theorem.
It is **not a full-code distance certificate**: nonzero inputs below
0.19N remain the main open class.

The new proof charges equations that force the inner state back to zero.
The earlier two-state support cover counted these returns but omitted
their rank contribution. No stronger assumption on random state mixing
is introduced.

## Construction and finite statement

Use the fixed GRS[32,22] outer over GF256, 32 independent uniform regional
permutations, and independent nonzero GF256 coordinate diagonals.
The inner has sixteen-symbol packets, four state symbols, and period-64
injection labels. Initial state is zero; final state is discarded.
The setup components have the independence stated in
[DUAL_MDS_PHASE.md](DUAL_MDS_PHASE.md).

The theorem applies both to that dual-MDS inner and to its separate
[fixed emission-scale variant](DUAL_MDS_OBSTRUCTIONS.md).
The proof uses the following structural properties: C A=0, minimum
nonzero emission weight thirteen, independent sets of up to four
injection columns, and invertibility of the one-output deterministic
update. Both constructions satisfy them.

For requested K=2^ell, put L=64 ceil(K/(22*64)), N=32L, and
e=ceil(N/10)-1. The probability is over ideal setup, and concerns
existence of any nonzero message in the indicated weight class.

| Requested size | Actual K | N | Covered inner-input weights | Class margin lower bound |
|---|---:|---:|---:|---:|
| K16 | 66176 | 96256 | 18289 through 23101 | 646 bits |
| K18 | 263296 | 382976 | 72766 through 91914 | 3616 bits |
| K20 | 1048960 | 1525760 | 289895 through 366182 | 15647 bits |

Each bound excludes output weight at most e in its entire listed class,
not just selected source rows, windows, or sampled input weights.
The independent dense theorem covers W >= ceil(0.24N), with conservative
margins 578/3381/14745 bits. A union bound therefore covers all W >=
ceil(0.19N) with margins at least 577/3380/14744 bits.
These numbers must not be reported as the code's distance margin.

## Fresh pivots at a zero-state checkpoint

Fix a route, regional orientations, and exact input/output supports I,J.
On the shortened outer space, retain one common span of all equations.
Discard a branch when a required nonzero input form lies in that span.
This is the alive-rank potential of
[REGIONAL_SPARSE_COVER.md](REGIONAL_SPARSE_COVER.md).

The zero-output phase reveals diagonals only on I minus J. In its
resulting quotient, the next state is

    S' = F_J S + sum_(j in I intersect J) c_j d_j ell_j.

All diagonals in this sum remain independently uniform and nonzero.
If the support cover requires the next state to be zero, append S'=0
to the common equation span. Do not merely require the formal state
map to have been zero already: the support cover describes witness
values, not identities of linear forms.

Let a=|I intersect J| and m=min(a,4). Select m injection columns as pivots;
they are independent by the MDS property. Reveal the other a-m diagonals.
A left inverse of the pivot columns produces m equations

    d_j ell_j + offset_j = 0.

Each offset is independent of every pivot diagonal. Apply the conditional
alive-rank lemma sequentially. For 1 <= tau <= 255, each equation
contracts expected alive potential by at most tau^-1. Ignoring any
remaining state equations is conservative. A return therefore contributes
an additional factor tau^-min(a,4), using no diagonal twice.

For a nonzero zero-state-to-zero-state packet, I=J and |I|>=5, so this
factor is tau^-4. For one output and no input/output overlap, F_J is
an invertible transvection; a nonzero incoming state cannot return.
The scaled variant forbids all no-overlap returns, but the certificate
does not use that stronger restriction.

## Positive transfer covering all support patterns

Let x mark input weight and y mark output weight. Take x>0 and 0<y<=1.
Write g=tau^-1 and index states by zero/nonzero witness state.

    T_00 = 1 + tau^-4 sum_(w=5)^16 binom(16,w) (xy)^w,
    T_01 = sum_(w=1)^16 binom(16,w) (xy)^w,
    T_11 = sum_(u=13)^16 binom(16,u) (gx+y+xy)^u.

For T_10, sum over support pairs with nonempty input and output, union
size at least thirteen, and not (|J|=1 and |I intersect J|=0).
A pair contributes

    x^|I| y^|J| tau^(-|I minus J|-min(|I intersect J|,4)).

The implementation expands this as a positive multinomial sum.
No cancellation-prone subtraction is used.
Start in zero state and accept either final state. With G=N/16 packets,
the generating-function upper bound is e_0 T^G (1,1)^T.
All actual support covers are included, including multiple excursions,
different row supports, and excursions crossing regional boundaries.

## Route averaging and the complete finite expression

For an exact input support with row occupancies r_i, the shortened outer
dimension is D=sum_i max(0,r_i-10). Occupancies one through ten are
infeasible for an exact nonzero row support.
Fix the total weight W and put alpha=W/N. For every regional support
profile and every positive lambda, the row-average upper bound is

    E_route[1_feasible tau^D]
      <= (L+1)^32 lambda^-W A(alpha)^L,

    A(alpha) = (1-alpha)^32
               + tau^-10 (1+(tau*lambda-1)*alpha)^32.

To prove it, temporarily replace each regional fixed-size subset by
independent Bernoulli indicators with the same mean. Conditioning each
column on its exact size costs at most L+1. Insert lambda to the total
weight; after conditioning its value is lambda^W. Independence of rows
gives a product of row polynomials. Bound the empty row and the relaxed
nonempty row separately, then apply arithmetic--geometric mean across
the 32 regional densities. This gives the displayed profile-independent
expression. Extra empty/short-row terms only increase the bound.

The conditional diagonal bound is uniform over routes, orientations,
and preceding equation spans. It can therefore be averaged separately
for each support pattern. A nonzero witness leaves equation rank at most
D-1. Combining all factors gives, for every fixed W,

    Pr[bad output with inner-input weight W]
      <= (L+1)^32 A(W/N)^L e_0 T^(N/16) 1
         / (tau * (x*lambda)^W * y^e).

No factor counts field-valued messages individually. The bound already
concerns existence of any message realizing the support pattern.
The constant-state or exact-state-distribution assumptions of a
fixed-message Markov analysis are not used.

## Directed replay and audit

The replay stores six rational witnesses and uses outward-rounded
80-digit interval arithmetic for positive matrix powers and all factors.
It partitions [0.19,0.24) into 500 rational density cells. Integer
endpoints cover each input weight exactly once. For a cell, the two
positive terms in A use their respective worst endpoints; the factor
(x*lambda)^-W also uses its correct endpoint. Multiply by the number
of weights in that cell, then sum every cell.

Floating optimization proposed the witnesses but is not part of the
default replay or the validity check. A separate exact rational test
checks the transfer using weight/intersection counts, and another
checks the positive matrix-power routine. Exhaustive-diagonal small-field
tests cover the return-pivot inequality for arbitrary dependent forms.
An independent reviewer audited the probability conditioning, return
pivots, complete weight coverage, and endpoint directions.

```sh
python -B -m unittest discover -s research/workstreams/char2_field_spin/theory -p test_return_cover.py
python -B research/workstreams/char2_field_spin/theory/return_cover_certificate.py
```

The next problem is the remaining input range 0<W<0.19N. The return-only
relaxation still gives positive exponents there. The joint rank/freshness
operator and the fixed emission scale address that gap; they are not
assumed by this certificate.
