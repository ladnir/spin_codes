# Larger-outer state-eight screen

The original screen below gave positive floating bounds at five active-row counts for GRS\([256,176]\). The follow-up directed replay now covers all counts 1 through 600 with 119 bits; the separate [aggregate certificate](FIELD256_K20_CERTIFICATE.md) closes the whole K20 code. The historical floating results are retained below. No native implementation or performance measurement was performed.

## Construction and finite geometry

Use the symmetric maps \(C=A^T\) from [STATE_CAPACITY_REPAIR.md](STATE_CAPACITY_REPAIR.md), with packet width \(16\), state dimension \(8\), and \(\eta=2\). Each packet has fresh independent nonzero input and state diagonals. The route consists of one independently permuted coordinate region for each outer coordinate.

At requested source size \(K=2^{20}\), the existing padding rule gives

\[
L=6016,\qquad N=256L=1540096,\qquad 176L=1058816.
\]

The full padded outer code has rate \(176/256=11/16\). Restricting to the requested source dimension still gives rate above \(2/3\). The bad-output threshold is

\[
e=\lceil N/10\rceil-1=154009.
\]

The screen bounds the full padded code, so it does not rely on a favorable choice of padding constraints.

## Executed two-state screen

For each exact active-row count \(a\), [superregular_return_screen.py](theory/superregular_return_screen.py) evaluates the previously audited rank cover. Write \(T_w(y,\tau)\) for its packet matrix at input weight \(w\), and define

\[
U_a=\sum_{h=0}^{a}\frac{\binom ah}{\binom Lh}\tau^h
[x^h]\left(\sum_{w=0}^{16}\binom{16}{w}x^wT_w(y,\tau)\right)^{L/16}.
\]

The resulting upper-bound expression is

\[
B_a=\binom La\tau^{-80a-1}y^{-e}\,e_0U_a^{256}\mathbf1.
\]

This expression concerns setup failure with exactly \(a\) active source rows. Its route averaging is without replacement. It is not a first moment over field words.

The search used its two existing starts and at most 25 optimizer iterations per start. The table records the selected witness. Its natural-log coordinates \((u,v)\) specify

\[
y=\exp(-e^u a/N),\qquad \tau=e^v.
\]

| Exact active rows \(a\) | Floating margin \(-\log_2 B_a\) | \(u\) | \(v\) | Selected iterations |
|---:|---:|---:|---:|---:|
| 16 | 1957.336892 | 7.183071336951522 | 2.947651636355612 | 11 |
| 32 | 1746.334103 | 7.120435500608540 | 2.258220251564710 | 9 |
| 64 | 1540.751175 | 7.062633053108466 | 1.991286947926587 | 9 |
| 128 | 1600.713707 | 7.024426560335106 | 1.907863435406201 | 9 |
| 256 | 2731.944273 | 7.028066385865181 | 1.923294185949148 | 9 |

All five selected optimizations reported success. The calculations were serial. Replay the bounded search with:

```powershell
python superregular_return_screen.py --active 16 32 64 128 256 --outer 256 176 --state 8 --k-log 20 --iterations 25
```

The working directory is `research/workstreams/char2_field_spin/theory`.

Positive margins at these five counts do not bound the intervening or larger classes. The calculation uses floating arithmetic, not outward-rounded intervals. It therefore does not establish a 40-bit whole-code guarantee.

## Independent rank-cap composition audit

The separate nine-state helper [state_rank_cap_transfer.py](theory/state_rank_cap_transfer.py) passed the conceptual composition audit. It is not an exact rank-transition matrix. Its contract applies to nonnegative continuation values \(V_0,\ldots,V_8\) satisfying

\[
V_1\le V_2\le\cdots\le V_8;
\]

there is no required comparison between \(V_0\) and \(V_1\).

For an incoming rank cap, the helper bounds every tilted upper tail after the state-root equations. Those tail bounds increase with the cap. Taking their successive differences produces a positive measure that dominates the actual law on this continuation cone. After the active-input zero equations, the next rank is at most \(\min(8,h+w-i)\); this is nondecreasing in \(h\). The return coefficient also increases with the incoming cap. Consequently each packet matrix preserves the cone.

Positive matrix products and nonnegative sums preserve that property. The normalized hypergeometric convolutions in [rank_cap_regional_screen.py](theory/rank_cap_regional_screen.py) therefore compose the packet contract without assuming an exact rank law. Starting from the all-ones terminal continuation is valid. Merely checking nonnegative matrix entries would not establish this argument.

The three new tests in [state_rank_cap_composition_audit_tests.py](theory/state_rank_cap_composition_audit_tests.py) pass. They compare normalized coefficients against an independent sequential Fraction polynomial, check regional scaling, and test all nine generating rays of the continuation cone under composition. The ten current local rank-cap tests also pass, for 13 passing tests in the final replay. These finite tests support the algebraic audit; they do not replace its uniform common-history argument.

The table above uses the two-state screen, not this nine-state helper. The next useful step is complete active-row coverage with directed replay, followed by a separate arithmetic-cost assessment of the larger outer. No whole-code or speed claim follows from this checkpoint.

## Complete floating cover of active counts 1 through 600

The follow-up screen covers every integer active-row count from \(1\) through \(600\). It uses 80 disjoint intervals and the unchanged coarse two-state packet matrix. After charging the number of classes in each interval and summing all interval bounds, the floating union margin is **119.379759 bits**. The weakest interval, \([108,113]\), has union margin 119.379837 bits.

This replaces isolated samples by complete coverage of that finite range. It remains a floating evaluation, not a directed certificate. Counts above 600 and the remaining density classes are not covered here. The larger-density proof gap for state eight and width sixteen is not resolved by this result.

For the interval argument, hold \(y,\tau\) fixed. The positive coefficient formula for \(U_a\) implies

\[
U_a\le U_b\quad\text{entrywise whenever }a\le b.
\]

Indeed, \(\binom ah\le\binom bh\) for every \(h\), and the additional coefficients are nonnegative. Positive products preserve this order. Thus an interval \([a_-,a_+]\) has the union bound

\[
(a_+-a_-+1)\,
\max_{a_-\le a\le a_+}
\left\{\binom La\tau^{-80a-1}\right\}
y^{-e}\,e_0U_{a_+}^{256}\mathbf1.
\]

The helper maximizes the source factor over every integer in the interval. It does not assume that the maximum is at the left endpoint. It evaluates the coefficient power once per fixed witness and reuses its coefficients for different interval endpoints.

The additional optimized anchors are \(1,2,4,8,512,600\). Each used the two existing starts and at most 25 iterations per start. Their selected floating margins are respectively 417.70, 835.49, 1673.12, 2572.79, 7671.98, and 9744.04 bits. Literal parameters are saved in [sparse_active_interval_screen.py](theory/sparse_active_interval_screen.py).

The helper also interpolates parameter pairs at rounded geometric midpoints of adjacent anchors. Interpolation only chooses candidate witnesses; it is not used as a proof of any intermediate bound. Each candidate defines its own fixed \(y,\tau\), and every accepted interval is evaluated explicitly. In particular, moving an interval endpoint never silently reuses an anchor's log-tilt with a different active count.

The module contains a compact literal fixture of all 80 endpoint triples. Its independent frozen replay reproduces the stated margin:

```powershell
python sparse_active_interval_screen.py --frozen
python -m unittest sparse_active_interval_screen_tests -v
```

All six tests pass. They check exact small-case regional and global monotonicity, coefficient reuse against the existing screen, fixed-mark conversion, interval-union accounting, source maximization, and complete fixture coverage. No optimization or native benchmark is needed for frozen replay.

## Directed sparse replay

The subsequent [Arb replay](theory/state8_sparse_certificate.py) retains
the same coarse packet cover and all 80 frozen active-row intervals.
Every accepted probability calculation uses rational marks and outward
balls; no floating value enters that arithmetic. It sums the full class
range with upper bound about `1.156409247462385e-36`, below `2^-119`.
The script exports its exact dyadic endpoint for the whole-code union.

Its five tests and an independent full replay pass. The audit additionally
checked 168 exact source-factor maxima, all 19 used anchors, directed
endpoints, precision restoration, truncation, and interval multiplicities.
The corresponding whole-code claim additionally uses the joint-return
middle bound and exact dense bound; sparse replay alone is not enough.

## Bounded four-state comparison

The optional comparison used the existing width-sixteen, four-state Cauchy diagonal kernel with the same GRS\([256,176]\) outer. Two starts and at most 20 iterations per start gave:

| Active rows | Floating margin | Log-tilt | Log-\(\tau\) |
|---:|---:|---:|---:|
| 64 | \(-3359.806012\) bits | 6.976259574863176 | 1.902042905410056 |
| 128 | \(-8029.733805\) bits | 7.025877073692352 | 1.955716799121857 |

Both optimizations reported success in eight iterations. This comparison used `joint_returns=False`; it does not test the stronger exact joint-return tables. These negative margins show that this particular upper bound does not close the sampled classes. They are not realized-route counterexamples or evidence of actual distance failure.

The next step is to resolve the remaining whole-code density gap before investing in a broad directed-certificate campaign. The sparse interval fixture is preserved for later reuse.
