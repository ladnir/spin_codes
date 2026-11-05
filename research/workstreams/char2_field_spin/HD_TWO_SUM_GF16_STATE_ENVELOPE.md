# A state-retaining all-cohort control for G16

This control keeps the two inner states while summing every group type.
It tests the frozen G16 two-sum construction with paired fair swaps and
the existing h4 inner. It changes neither the encoder nor its setup law.

The bounded FLOAT screen covers every integer active-group count
`18..32`, every constituent-row set, and every source support. **It does
not close this domain:** the best selected union is
`+16494.53668044997` log2 bits. This is a failed bound, not a construction
counterexample. No directed certificate or new native run is claimed.

## Regional envelopes and fixed group types

Fix the weight marker $0<z<1$ and an active-group count $r$. Let

\[
R_{ij}(m)=\max_{0\le f\le\lfloor m/2\rfloor}R_{ij}(m,f),
\qquad i,j\in\{0,1\},\quad 0\le m\le32,
\]

where $R(m,f)$ is the existing forced-nonzero-amplitude regional upper
matrix. The maximum applies uniformly to every fixed source subset;
it does not resample that subset. Select positive constants satisfying

\[
R_{ij}(m)\le A_{ij}B_{ij}^{m}
\]

at every occupancy. The implementation selects slopes and computes the
smallest intercepts over the complete stored curves.

A group type is its number of active constituent rows, from 1 through 16.
Types 1 and 2 use the exact singleton and pair laws. Types 3 through 16
use the fixed-row-set weights $W_a(w)$ from
[the arity reduction](HD_TWO_SUM_GF16_ARITY_DENSE.md). Those weights already
include the row-set multiplicity and the dense counting normalization.

Consider one path through the 16 regional transitions of an outer column.
Let $n_{ij}$ count its transitions from $i$ to $j$. Its 16 occupancy
marks are the multiset containing $n_{ij}$ copies of $B_{ij}$. Write

\[
P=\prod_{\ell=1}^{16}B_\ell,
\qquad e_w=e_w(B_1,\ldots,B_{16}),
\qquad A^n=\prod_{i,j=0}^{1}A_{ij}^{n_{ij}}.
\]

Fix floor markers $x_1,x_2,y_a\ge1$. Each marker remains fixed for its
type across every state path and every outer column. The type weights are

\[
\begin{aligned}
F_1(n)&=16^{1/16}255^{-1/2}x_1^{-9/16}
       (1+255x_1P),\\
F_2(n)&=120^{1/16}255^{-1}x_2^{-18/16}
       \left(1+510x_2P+255x_2^2(239P+e_{15})\right),\\
F_a(n)&=y_a^{-9/16}
       \left(W_a(0)+y_a\sum_{w=1}^{16}W_a(w)e_w\right),
       \qquad 3\le a\le16.
\end{aligned}
\]

The singleton floor is nine active row coordinates. For a pair, total
row weight at least eighteen relaxes the two individual floors of nine.
For a dense group, at least nine columns are nonempty. The factors
16 and 120 count the fixed sparse row choices once across the whole word.
No additional row-set or mask binomial belongs in these expressions.

## Summing fixed types without a common potential

Fix an ordered type vector \(\boldsymbol\tau\in[16]^r\). After the
fixed-row-set comparison, let \(M_{\boldsymbol\tau}\) be its column upper
matrix. Expanding the regional state paths gives

\[
M_{\boldsymbol\tau}(i,j)
\le \sum_{\pi:i\to j} A^{n(\pi)}
             \prod_{g=1}^{r}F_{\tau_g}(n(\pi)).
\]

For each path through the 16 column boundaries, Hölder sums the same
type vector across all columns. Define

\[
K_r(i,j)=\left(\sum_{\boldsymbol\tau}
                   M_{\boldsymbol\tau}(i,j)^{16}\right)^{1/16}.
\]

Then the type union is at most \(e_0K_r^{16}\mathbf1\). Minkowski over
the regional paths gives the computable entrywise upper matrix

\[
\widehat K_r(i,j)=
\sum_n \operatorname{mult}(i,j,n)A^n
       \left(\sum_{a=1}^{16}F_a(n)^{16}\right)^{r/16}.
\]

The ordered group-type sum factors as
\(\sum_{\boldsymbol\tau}\prod_gF_{\tau_g}^{16}
=(\sum_aF_a^{16})^r\). It already counts all cohort colorings.
The resulting failure-probability screen is

\[
\binom{8192}{r}\frac{z^{-734003}}{255}
 e_0\widehat K_r^{16}\mathbf1.
\]

There is one projective factor, no endpoint-potential factor, and no
extra witness-bank factor. A witness may be selected separately for each
disjoint weight class before summing those bounds.

Only 274 endpoint/count classes occur: 65 for each of endpoints 00 and
11, and 72 for each of 01 and 10. Grouping by transition counts is valid
because the four envelopes are the same in every region and the type laws
are symmetric in the 16 labelled positions. Position-dependent regional
envelopes would require a different grouping. Matrix chronology itself
is retained by the exact path multiplicities and the final matrix power.

## Bounded controls and the remaining loss

The first curves-only control used $z=.99$ and $z=.9995$. At `.99`,
the live-to-live log entry ranges only from about `-106.7500` at occupancy
zero to `-106.5693` at occupancy 32. At `.9995`, those endpoints are about
`-4.8025` and `-2.2822`. The latter mark therefore gave too little live
contraction to motivate a full gate. The tested marks were changed to
`.99`, `.98`, and `.96`, retaining the stronger marks from the earlier
scalar-envelope control. The original curves artifact preserves its
original helper hash.

The final gate used four shared floor-marker logarithms: `0,1,4,8`.
For each of its twelve schedules, it made one deterministic coordinate
pass over finite-point secant candidates. Every candidate intercept was
checked at all occupancies. The run built three regional banks and took
`7.45288520003669` seconds. It used 44 evaluations for each `.99` schedule
and 46 for each remaining schedule; no $z$ or floor-marker fit occurred.

| z | floor log 0 | floor log 1 | floor log 4 | floor log 8 |
| --- | ---: | ---: | ---: | ---: |
| .99 | +16494.5367 | +16718.5367 | +17390.5367 | +18286.5367 |
| .98 | +20488.0884 | +20326.0884 | +20777.9050 | +21673.9050 |
| .96 | +42322.8638 | +42160.8638 | +41674.8638 | +41026.8641 |

Every weight still fails. The smallest selected weight bound is
`+8311.382226380047` at $r=21$; the selected union is dominated by
$r=32$. Compared with the separate cut control's `+16580.99948205441`,
the union improves by only about 86.46 bits.

The best global schedule is `.99` with floor marker one. Its 00 envelope
touches occupancies 0 and 32; its off-diagonal envelopes touch 2 and 32.
The 11 curve is essentially affine throughout this range. Thus retaining
the states did not remove the principal loss at these marks.

Two relaxations remain distinct. Positive off-diagonal exponential
envelopes introduce births or returns at occupancy zero, although both
actual entries there vanish. Also, the 00 envelope bridges the exact
value one at occupancy zero and the rare-return values at positive
occupancy. The comparator's all-zero paths alone contribute

\[
\binom{8192}{r}\frac{z^{-734003}}{255}A_{00}^{256}
\left(\sum_a F_a(B_{00},\ldots,B_{00})^{16}\right)^r.
\]

This expression is a necessary diagnostic for the bound, not a lower
bound on the construction's failure probability. The additional Hölder
and Minkowski steps can also lose consistency between the dominant types
on different paths. No more alpha-one fitting is planned for this control.

## Fractional moments: the next interface

A fractional extension may address rare setup events that an alpha-one
counting bound cannot suppress. It requires a new forced-nonzero packet
interface; the existing uniform-amplitude fractional provider is not a
substitute. No fractional gate is implemented or claimed here.

The reviewed reduction first separates each mask mass from its normalized
nonzero-amplitude law. After subadditivity over mask words, Jensen averages
all inner keys at a fixed route and mask word. The existing averaged-refresh
matrix product then bounds that mean. State-path expansion puts fractional
powers on packet entries before averaging placements and swaps. Arrangement
counts and route normalizers keep exponent one.

Fixed row sets remain unioned outside the fractional power. The dense
weights therefore use the moment exponent $16\alpha$, not $W_a^\alpha$.
Directly powering point masses would instead require an additional
$255^{(1-\alpha)|T|}$ factor before normalized amplitude averaging.
The [fractional-mask note](HD_TWO_SUM_GF16_FRACTIONAL_MASKS.md) gives the
detailed chronology and the remaining implementation obligations.

The next proof step is to implement and test this fractional local interface,
not to promote the failed alpha-one controls or change the native code.

## Files, tests, and pins

The isolated helper is
[hd_two_sum_gf16_state_envelope.py](theory/hd_two_sum_gf16_state_envelope.py).
Its [eight tests](theory/hd_two_sum_gf16_state_envelope_tests.py) passed in
`0.008s`. They check exhaustive binary path counts, labelled symmetric
polynomials, sparse and dense type factors, full-curve domination, exact
zeros, fixed-type Hölder/Minkowski, one-type equality, matrix powers,
and resource guards. Two independent read-only audits found no material
issue; neither repeated a regional gate.

```text
state_envelope.py
5107b945ab18415191c48a3434a354c101d8fc5aeda66ec83a3fd2a5c3498af7
state_envelope_tests.py
886d3e979025610c27364acd4e0c2a12fffaac527e30446368a14a46a4780cf0
hd_two_sum_gf16_state_curves_r18_32_float.json
3e561e933be2d164af6c123786694eff45b861685370f0c546764cff48f1a908
hd_two_sum_gf16_state_r18_32_round1_float.json
112dcf3e2f363bf2ee3ef8e7938a5f6267576138b09432025dfd5aadbfbb3474
```

The completed gate is reproducible from the theory directory with
`python -B hd_two_sum_gf16_state_envelope.py --rounds 1 --output NEW.json`.
The output path must be new. The default computation guard is 30 seconds.
