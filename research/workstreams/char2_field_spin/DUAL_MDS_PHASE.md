# Dual MDS phase maps with the existing four-state cost

A fixed reversal of the four state coordinates repairs the period-64 emission geometry.
Both the injection columns and emission rows then have rank four on every four distinct labels.
For an unchanged nonzero state, the 64-symbol emission has minimum weight 61, rather than the archived map's minimum weight 58.
The change adds no field products and retains one regional route and one orientation per region.
It does not establish a whole-code distance certificate for the proposed GRS[32,22] profile.

## The invariant pairing

Use \(F=\mathrm{GF}(256)\) in the AES representation and its embedded subfield \(H=\mathrm{GF}(16)\).
For \(u\in F\), define

\[
v(u)=(1,u,u^2,u^3)^T,
\quad
J=\begin{pmatrix}
0&0&0&1\\
0&0&1&0\\
0&1&0&0\\
1&0&0&0
\end{pmatrix}.
\]

The matrix \(J\) is symmetric, alternating, and invertible, with \(J^2=I\).
Its pairing on the cubic curve is

\[
v(a)^TJv(b)=a^3+a^2b+ab^2+b^3=(a+b)^3.
\]

Let \(T_\beta\) be the existing translation matrix,

\[
T_\beta=
\begin{pmatrix}
1&0&0&0\\
\beta&1&0&0\\
\beta^2&0&1&0\\
\beta^3&\beta^2&\beta&1
\end{pmatrix}.
\]

It satisfies \(T_\beta v(u)=v(u+\beta)\), \(T_\beta^{-1}=T_\beta\), and

\[
T_\beta^TJT_\beta=J,
\qquad JT_\beta=T_\beta^TJ.
\]

The pairing identity proves the first equality: translating both labels leaves \((a+b)^3\) unchanged.
Four distinct curve vectors span \(F^4\), so equality on curve vectors determines the matrix identity.

Set

\[
C_0=(v(u))_{u\in H},\qquad A_0=C_0^T.
\]

As before, \(C_0A_0=0\), since all subfield power sums of degrees zero through six vanish.
The revised local emission map is \(E_0=A_0J\); the local injection map remains \(C_0\).
Thus \(C_0E_0=0\).

## Global and local interfaces

For region \(r\), fix its invertible orientation \(B_r\in\mathrm{GL}(4,F)\).
Use the existing four phases \((0,2,6,4)\), whose cosets of \(H\) are disjoint.
At phase \(\beta\), define

\[
C_{r,\beta}=B_rT_\beta C_0,
\qquad
E_{r,\beta}=A_0JT_\beta B_r^{-1}.
\]

For a sixteen-symbol input packet \(X\) and incoming global state \(Q\), compute

\[
Y=X+E_{r,\beta}Q,
\qquad Q'=Q+C_{r,\beta}X.
\]

An optional input diagonal is applied before this step, as in the archived construction.
The product \(C_{r,\beta}E_{r,\beta}\) is zero.
No transpose relation in the ordinary coordinate pairing is assumed between these two maps.

The existing local coordinates \(S=T_\beta B_r^{-1}Q\) give

\[
Y=X+A_0JS,\qquad S'=S+C_0X.
\]

Changing phase from \(\beta\) to \(\gamma\) still applies \(T_{\beta+\gamma}\).
At a region boundary, the prepared map is still
\(B_{r+1}^{-1}B_rT_\beta\), where \(\beta\) is the preceding region's last phase.
The initial conversion remains \(B_0^{-1}\); the final conversion remains \(B_{\mathrm{last}}T_{\beta_{\mathrm{last}}}\).
Thus the previous phase, regional-boundary, and external-boundary schedules remain valid.

There is also an equivalent formulation with literal dual-translation emission \(A_0T_\beta^T\).
Take injection \(JT_\beta C_0=T_\beta^TJC_0\) and local coordinates \(T_\beta^TQ\).
Its local injection is \(JC_0\), and its local emission is \(A_0\).
The first formulation is preferable when preserving the existing forward phase-change schedule.

## Rank and support guarantees within one label period

Write \(a=u+\beta\) for a label in the 64-element union of cosets.
The corresponding injection column and emission row are

\[
B_rv(a),\qquad v(a)^TJB_r^{-1}.
\]

Multiplication by the fixed invertible matrices does not change ranks.
Every set of \(m\) distinct labels therefore gives injection rank and emission rank
\(\min(m,4)\).
For four labels, the Vandermonde determinant is the product of their six nonzero pairwise differences.
Both period maps are consequently MDS.

For a fixed nonzero state, the emission is evaluation of a nonzero polynomial of degree at most three at 64 distinct points.
It has at most three zeros, and its weight is at least 61.
Equality is attained by a cubic vanishing at any three of the labels.

More generally, fix a state subspace \(W\subseteq F^4\) of dimension \(d\), where \(1\le d\le4\).
At most \(4-d\) labels can vanish for every vector in \(W\).
Hence the union of the emission supports of \(W\) has size at least

\[
64-4+d=60+d.
\]

This bound is attained by the subspace of polynomials divisible by the product of any \(4-d\) selected linear factors.
It is a deterministic subspace-support statement, not an assumption that the state is random.
For an arbitrary fixed state subspace, the restriction to \(m\) selected labels has rank at least
\(\max(0,d-4+\min(m,4))\); its exact rank can still depend on that subspace.

These guarantees apply across consecutive empty input packets, when the global state remains unchanged within a region.
A run containing complete label periods contributes at least 61 nonzero outputs per complete period.
For an additional run of \(r<64\) distinct consecutive labels, the corresponding lower bound is \(\max(0,r-3)\).
Active packets change the state and cannot be merged into one unchanged-state polynomial evaluation.
The rank guarantee also concerns distinct labels: the same 64 labels repeat in later periods.

If a state subspace of dimension \(d\) is fixed before a fresh uniform regional orientation, an additional exact probability statement is available.
For \(m\le4\) selected labels, the probability that all their emission functionals vanish on the oriented subspace is

\[
\frac{{4-m\brack d}_{256}}{{4\brack d}_{256}},
\]

where the numerator is zero for \(d>4-m\), and the brackets denote Gaussian binomial coefficients.
The formula counts uniformly oriented \(d\)-subspaces contained in a fixed \((4-m)\)-dimensional kernel.
It must not be applied after conditioning on history that makes the subspace depend on that same regional orientation.

## Adjoint and cost

The new local step is not self-adjoint.
For output dual \(Z\) and next-state dual \(H\), its transpose is

\[
X^*=Z+A_0H,\qquad S^*=H+JC_0Z.
\]

The forward step reverses the state before emission; the transpose reverses the moment output before adding it to the state dual.
Reusing a forward local-step routine as the adjoint would be incorrect.
Internal phase adjoints remain \(T_\delta^T\), and the existing prepared boundary matrices still use their ordinary transposes.

The reversal is a fixed permutation of four state records.
With separate state registers, it can be expressed by static argument or register naming.
A packed state layout may require a lane permutation, but no additional field product is needed.
The field-product counts therefore remain:

- twelve GFNI512 products for both width16/four-state local maps;
- four more GFNI512 products if an input diagonal is used;
- four GFNI128 products per internal phase change;
- four GFNI512 products per prepared dense regional boundary.

For the previously costed GRS[32,22] outer, the leading input-diagonal-inclusive total remains
\(177/88\) GFNI512 processed-byte equivalents per source symbol.
Its separate counts are \(169/88\) GFNI512 and \(4/11\) GFNI128 per source symbol, before the finite boundary correction.
These counts inherit the earlier exclusions for lane shuffles, reductions, routing, and stores.
They are not measured throughput.

## Which previous proof bounds carry over

Keep the GRS[32,22] outer, route, nonzero input diagonals, independent uniform regional GL4 orientations, and period64 rounding unchanged.
The following previously proved classes survive the change to \(J\); none is a whole-code guarantee.

The [one-row return-count argument](RATE_DESIGN_REVIEW.md#exact-one-row-bound) uses only three local properties.
Each nonzero state emits at least thirteen symbols in an empty packet.
A single nonzero input has a nonzero injection.
Before that region's only active packet, a fresh regional orientation makes a nonzero incoming state uniform on \(F^4\setminus\{0\}\).
All three properties are retained: reversal preserves the packet emission code, injection is unchanged, and phase changes are invertible.
The conditional return probability remains \(1/(256^4-1)\).
Dropping the active packet still leaves at least \(13(L/16-1)\) outputs in every intermediate active region.
Consequently the existing conservative one-row bound is unchanged; it does not yet exploit the stronger 61-per-period emission bound.

Replaying `return_count_upper` with length32, dimension22, state4, and period64 gives:

| Nominal size | Rounded K | Encoded N | Approximate one-row margin |
|---|---:|---:|---:|
| K16 | 66176 | 96256 | 74.627302 bits |
| K18 | 263296 | 382976 | 72.634996 bits |
| K20 | 1048960 | 1525760 | 70.640794 bits |

The failure upper bounds are exact rational numbers; the displayed negative logarithms are decimal approximations.

The [dense-class theorem](REGIONAL_GLOBAL_ROUTE.md) allows any invertible inner independent of the route and input diagonals.
The revised inner remains causally invertible: from a known incoming state, recover the current input from its identity-feedthrough output, then update the state.
Its independent setup law is unchanged.
Thus the existing GRS[32,22] dense cutoff \(\alpha\ge0.24\) and conservative margins 578/3381/14745 bits remain valid for K16/K18/K20.

The [unequal common-support window event](REGIONAL_CAPACITY_COMPARISON.md) uses only the route, outer shortening dimensions, and four state constraints per required return.
Those quantities are unchanged.
Its K20 event-family upper-bound margin therefore remains approximately 90.85 bits.
This is a bound only for rows sharing one outer support and fitting one window per selected region.
It does not cover all sparse bad words, all support profiles, or multiple disjoint windows per region.

Together, the one-row, dense-class, and specified-window results still leave an unproved sparse multirow class.
The dual-MDS rank improvement is an additional input to that missing argument, not a justification for silently enlarging the scope of any carried-over bound.

## Verification and remaining gate

A read-only scalar check used the actual embedded GF16 geometry from `region64_reference.py`.
It verified the symplectic identity at all 256 offsets, all 4096 pairings of the 64 labels, and their distinctness.
It also checked 43 packet instances of local/global equivalence and the extended adjoint, including nonzero initial and terminal states and varied region boundaries.
All checks passed.
The replayable reference and eight-test suite now supplement those original checks:

```sh
python -B research/workstreams/char2_field_spin/theory/dual_mds_phase.py --check
python -B -m unittest discover -s research/workstreams/char2_field_spin/theory -p test_dual_mds_phase.py -v
```

The suite verifies all translation offsets, the local zero product, direct translated-label maps, 36 varied region/diagonal equivalence and extended-adjoint cases, and the tight generalized-support examples.
It exhausts all subsets of sizes one through four from a sixteen-label set spanning all phases and tests representative oriented minors from the full period.
The nonzero pairwise-difference product proves all distinct-label four-minors; the tested subset is not being used as a substitute for that proof.
A negative regression deliberately substitutes the forward local step into the adjoint and detects the resulting failure.

This repair strengthens the deterministic emission geometry while retaining the injection geometry and multiplication counts.
It does not change the state dimension or the outer shortening dimensions.
Existing whole-subspace obstructions based only on those dimensions therefore remain relevant.
The next gate is to use the exact distinct-label rank and subspace-support statements in the whole-code analysis of GRS[32,22].
No existing reference, production kernel, route, or benchmark was changed.
