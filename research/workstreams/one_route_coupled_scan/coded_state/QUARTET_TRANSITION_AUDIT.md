# Why the partial-state moment remains large

The fixed quartet construction and its 45%/40-bit goal are unchanged.
This audit distinguishes interpolation loss from overcounting in the
[averaged-threshold routing bound](QUARTET_BURST_CUT.md).
It does not establish a whole-code distance certificate or a bad codeword.

## Exact cells: interpolation is not the only loss

For a support history S, write A for total selected state coordinates,
C for covered inactive inputs, and D for deaths. The routing score is
F(S)=A/4-C-(37/50)D. At r=2100 inactive rows, the existing field argument
would close if every support satisfied F(S)<=18171.
The burst certificate reduces this sufficient condition to bounding

    G(S)=sum_t[L_t max(u_t,v_t)/512-C_t]
         -(37/50)[D-sum_t max(u_t-v_t,0)] <=14469.

Here u_t,v_t are consecutive state weights and L_t is the actual inactive
input count at update t. All randomness belongs to the once-sampled setup.

`quartet_transition_audit.py` computes directed positive-coefficient
intervals for individual weight transitions. Rational tilts only improve
conditioning; integer recurrences and rational logarithm bounds verify the
result. No floating-point optimizer value is trusted.

At p=525/2048 and log2(b)=8, the log2 kernel intervals and old chord values
are approximately:

| Incoming/outgoing weight | Directed cell interval | Chord value |
|---|---:|---:|
| 128/128 | [125.6142,125.6152] | 257.0789 |
| 256/256 | [202.7374,202.7388] | 202.7381 |
| 384/384 | [141.2564,141.2574] | 212.1988 |
| 128/384 | [611.3904,611.3920] | 655.2114 |
| 384/128 | [451.7275,451.7295] | 507.9875 |

The table decimals are summaries; the helper emits directed rational
endpoints. The chord values are FLOAT screens. The actual quartet-erasure
filter is retained in every cell. At weight256 the chord is already an
anchor, so interpolation cannot explain that large entry.

A joint tilt also imposes the incoming weight, outgoing weight, and regional
inactive quota. Sixteen equal-weight transitions are combined before
conditioning on exactly2100 inactive inputs among8192 positions. The
directed log2 upper bounds for this region are approximately1924.83,
2958.46,2266.63 at weights128,256,384 and log2(b)=8. At log2(b)=4 they are
1597.82,1544.62,1210.28. These remain loose positive bounds. They are not
exact conditioned moments or probabilities of bad codewords.

## An exact harmless subcount

We next lower-bound part of the moment to test whether a substantial count
can arise without a routing violation. Fix ONE inactive-row set of size2100.
Use weight256 at every update of a complete region. For a quartet, let
(a,v,t,c,d) denote incoming weight, outgoing weight, inactive labels,
covered inactive labels, and deaths. The helper enumerates the integer
multiplicity m(a,v,t,c,d), including the existing local erasure filter.

Two integer histograms over128 quartets have totals

    type0: (256,256,131,33,44),
    type1: (256,256,132,34,43).

Use type0 in twelve specified updates and type1 in four specified updates.
The region then has A=4096, exactly2100 inactive labels, C=532, D=700.
Consequently

    G_region=2100/2-532-(37/50)700=0,
    F_region=4096/4-532-(37/50)700=-26.

These histories contribute weight one to b^G for EVERY b>0. Their routing
score is negative, not a violation.

For histogram counts n_j, its integer coefficient is

    H(n)=128! product_j m_j^n_j / product_j n_j!.

To obtain a substantial subcount, the frozen witness also permits disjoint
histogram moves. Each move replaces two local types by two others with the
same five statistics. Its exact factorial ratios are summed over every
feasible integer displacement. Disjoint sets of type indices make these
sums multiply without overlap. The helper verifies the statistics,
nonnegativity, disjointness, and integrality using exact arithmetic.

Let H0,H1 be the resulting coefficient sums. For any fixed incoming support
of weight256, the expected number of these regional continuations is

    Z=H0^12 H1^4 / [binom(512,256)^16 binom(8192,2100)].

The first normalization averages the fresh uniform quartet partitions;
permutation symmetry makes it independent of the particular incoming
support. The second normalization imposes the regional routing quota.
The directed result is

    274061728965881/1099511627776 <= log2 Z
      <=274061728965883/1099511627776,

or about249.258 bits. This is an expectation of a COUNT, not a probability.
Different histories may occur under the same setup.

The witness extends to valid full-length support histories. In the first
region, keep the state empty for15 updates and select256 coordinates in
the last update. Require all2100 inactive labels to occur in the preceding
7680 positions. This prefix has C=D=G=0 and F=64; its expected count is

    binom(512,256) binom(7680,2100)/binom(8192,2100).

Append255 regions from the subcount. Independence across regions and the
incoming-support symmetry justify multiplication. The terminal weight is
256, with no terminal constraint. Every resulting full history has
G=0 and F=64-255*26=-6566. Its expected count has directed log2 lower bound

    35096826360308999/549755813888 >63840.

Thus, at log2(b)=4, the unpruned first-moment expression
b^-14469 E[sum_S b^G(S)] already exceeds2^5964, even for ONE row set.
This proves a limitation of that particular expression. It does not rule
out larger markers, selective counts, or other proof methods. The witness
also has C+D=314160, above the remaining bad-codeword turnover threshold
255463 at r=2100. The geometry-only count does not use that exclusion.

## A selective formulation: erase nonpositive components

Fix the complete routing and quartet partitions. Regard each selected
state coordinate (t,j) as a vertex. Connect two selected vertices when
they are consecutive instances of the same coordinate, or when they both
cover one inactive input. The latter input has at most five state neighbors:
its outgoing coordinate and the four incoming coordinates of its quartet.
Let S_1,...,S_m be the connected components of this selected graph.

The score is exactly additive:

    F(S)=sum_i F(S_i).

To see this, let P count selected consecutive-coordinate pairs and let w
be terminal weight. The empty initial state gives D=A-w-P. Therefore
F=(1/4-37/50)A+(37/50)(w+P)-C. Vertex terms and consecutive-pair terms
belong to single components. Every covered inactive input also belongs to
exactly one component, by the definition of the edges. This proves the
identity, including the free terminal boundary.

Choose a global maximizer of F with the smallest number of selected
vertices. Every nonempty component has strictly positive score: deleting
a negative component improves F, and deleting a zero component preserves
F with fewer vertices. Moreover, F is supermodular, since coverage is
submodular and selected-pair terms are supermodular. The intersection of
all maximizers is consequently the unique inclusion-minimal maximizer.

This is a deterministic reduction, not a tail bound. It suggests counting
only positive connected witnesses or analyzing the canonical minimum cut.
One cannot replace a collection of positive components by one component
with score above18171: many smaller positive scores can accumulate.
Nor can one multiply component probabilities without proving the required
dependence statement under the regional quotas.

## Reproduction and next discriminating test

`quartet_neutral_type.py` defaults to exact replay of the frozen histograms
and moves. Optional `--propose` is a bounded NumPy search; proposal values
never establish a bound. The default replay needs no SciPy or optimizer.
It also checks4096 local masks,65536 component identities,2048
supermodularity inequalities, and component positivity in exhaustive toy
maximizers. The transition helper checks49 directed cells against exact
small polynomials and two integer root inequalities.

All runs were serial, remote, file-free, and resource-bounded. A missing
optional SciPy installation was handled by removing that dependency, not
installing packages. The frozen prior proof helper remains unchanged.

Next derive a quota-aware bound for collections of positive components, or
for the deficit of the constructive flow. Test whether that formulation
avoids the explicit harmless subcount before attempting full-size numerical
closure. Preserve the field exclusions and all-row union throughout. The
entire r=2100 class and the other unresolved row counts remain open.
