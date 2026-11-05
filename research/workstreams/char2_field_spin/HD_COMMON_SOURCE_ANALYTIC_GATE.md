# A good-route regional return diagnostic

This diagnostic asks whether a route event common to all masks improves
one regional 0-to-0 bound. It fixes 18 potential source groups, nine actual
inputs, `z=691/700`, and `alpha=1/5`. It treats the two potential-source
histograms `(s,f)=(18,0)` and `(16,1)` separately. Here `s` counts singleton
source pairs and `f` counts full source pairs.

The helper reuses the saved unpowered regional curve from
`theory/hd_two_sum_late_alpha_quiet_gate_float.json`. No provider is rebuilt.
The curve's occupancy-nine entry is `log2(Rall)=-108.20440173434696`.
Its zero-intercept straight-line envelope is a distinct, looser control.
Both values are FLOAT proposals, not directed bounds.

## The common route event

Fix a potential source set `G` and an actual support `U` contained in `G`.
The region uniformly permutes its 4096 prescribed source pairs and swaps
each pair independently and fairly. Call the route good when every
16-coordinate emission block contains at most four positions from `G`.
This event is the same for every `U` contained in `G`.

The exact probability `p_G` is the positive coupon coefficient derived in
`HD_COMMON_SOURCE_ROUTE_PARTITION.md`:

\[
p_G=\frac{[X^sY^f]H_4(X,Y)^{256}}
 {2^s\binom{4096}{s}\binom{4096-s}{f}}.
\]

A routed pair spans two distinct emission blocks. Thus the half-block
polynomial, rather than an independent-position occupancy model, is needed.
The existing `good_probability` implementation computes this probability
as an exact rational. It does not resample `G` for different masks.

## A pointwise refresh bound and the return normalization

Write `q=256`, `D=q^16-1`, and `beta=D/255^16`. The live component of the
two-state transfer is a density bound relative to the uniform measure on
the `D` nonzero states. It is not a probability state with a stochastic row.

The fresh refresh is `D2 H16 D1`, with independent nonzero diagonal entries
and a superregular `H16`. Fix a nonzero input with support size `t` and a
nonzero target with support size `u`. The output-zero constraints leave at
most `t+u-16` free input coordinates; if `t+u<=16`, the event is impossible.
An information set among the input coordinates gives at most
`255^(t+u-16)` admissible full-support assignments. The input diagonal
assigns each one probability `255^-t`. The output diagonal maps a vector
with the required support to the target with probability `255^-u`.
Consequently every nonzero target atom is at most `255^-16`.
This holds for each fixed refresh input, not only an average flat input.

Now fix a current packet input `X`. Its nonzero amplitudes are independent
of the incoming state, conditional on the complete outer mask word.
The pre-refresh state is `Q+CX`, since `CA=0`. Under the flat reference
measure, a return to zero requires the single value `Q=-CX`.
Its mass is either zero or `1/D`. Therefore the live-to-zero entry is at
most `1/D`; there is no additional factor `beta` in this entry.
The factor already carried by the incoming density is propagated by matrix
multiplication. The formula `packet_ret/D` in `hd_xor_pair_route.py` uses
the same normalization, with a sharper marked return numerator.

On a good route, each actual block has at most four inputs. The corresponding
feedback columns are independent. A nonempty packet therefore cannot take
zero state to zero state. Dropping output marks only enlarges the other
entries, so the following matrix is a valid nonempty-packet cap:

\[
A=\begin{pmatrix}0&\beta\\1/D&\beta\end{pmatrix}.
\]

The live-output entries use the pointwise refresh bound above after bounding
the pre-refresh weighted mass by one. They need not sum to at most one.

## Empty packets and regional returns

An empty packet preserves zero state. A nonzero incoming state emits a
nonzero RS[16,4] word in at least one block, hence at least 13 nonzero bytes.
Its pre-refresh state remains nonzero. Its cap is therefore
`diag(1,beta*z^13)`. The helper verifies the exact rational inequality
`beta*z^13<=1` before deleting any empty packets.

Deleting these diagonal contractions preserves an entrywise upper bound
and packet chronology. There are at most nine nonempty packets. Hence a
good regional 0-to-0 cap is

\[
C_{\rm good}=\max_{2\le k\le9}(A^k)_{00}.
\]

One nonempty packet contributes zero. Two can contribute
`beta/D=255^-16`: the state can be born and later return. Thus neither
distinct packet placement nor the common good event zeros the regional
0-to-0 entry. The exact powers retain these paths.

## Conditional Jensen comparison

Let `F(Pi)` be the nonnegative regional marked 0-to-0 transfer after
averaging the independent amplitude and inner keys, with route `Pi` fixed.
The good-route argument bounds `F(Pi)` by `Cgood` on the common good event.
The saved unconditional max-over-`f` cap supplies `E F(Pi)<=Rall` for every
fixed actual support of size nine. This is deliberately weaker than
specializing the actual matching histogram to the potential source class.

Positivity gives `E[F(Pi) 1_bad]<=Rall`. Conditional Jensen therefore gives

\[
\mathbb E_\Pi F(\Pi)^\alpha
\le p_G C_{\rm good}^{\alpha}
 +(1-p_G)^{1-\alpha}R_{\rm all}^{\alpha}.
\]

The probabilities remain outside the conditional alpha powers.
The helper compares this sum with the unconditional Jensen cap
`Rall^alpha`, and records their minimum. It also reports the same comparison
against the saved straight-line envelope, with a separate label.

This calculation concerns one fixed actual mask and one regional entry.
It does not sum the masks in a complete column coefficient before alpha.
It does not cover other boundary states, outer syndromes, or whole words.
A successful regional diagnostic would motivate that next conditional
column interface; it would not itself provide a distance certificate.

## Implementation and review status

`theory/hd_shifted_common_source_analytic.py` pins the saved artifact with
SHA-256 `1f9bcb13d4a5120a130e3b23dfad264ffe25a3fc42f9b5c2dc4e8488732d4521`.
It also checks the reviewed provider and exact coupon-helper identities.
The output keeps exact rational probabilities, cap entries, and all nine
matrix powers. Only logarithms and the fractional comparison use FLOAT.

Seven light tests pass in 0.028 seconds. They include an exhaustive F5
pointwise refresh example, flat-state return normalization, exact path
powers, noncommuting empty-packet deletion, a tiny route probability,
conditioning exponents, and source/scope guards. The production scope test
uses mocked event probabilities. Independent source review found no issue.
The root task subsequently passed all 13 new and coupon-reference tests
and ran the frozen diagnostic once in 0.0045021 seconds.

Reproduce the single frozen diagnostic from `theory`, using a fresh output:

```powershell
C:/Python314/python.exe hd_shifted_common_source_analytic.py --output hd_shifted_common_source_analytic_gate_float.json
```

The helper refuses to overwrite its output. The saved artifact and inputs
have the following SHA-256 pins:

| File in `theory` | SHA-256 |
|---|---|
| `hd_shifted_common_source_analytic.py` | `6172f44ed3cea163f742d74bfac7038b6882b498eba79666024806fabad29093` |
| `hd_shifted_common_source_analytic_tests.py` | `d55cc40a48ebc4ef6234e4d4f9ad76bb459a9be4a3e48bd2a307d015a7ff4f6d` |
| `hd_shifted_common_source_analytic_gate_float.json` | `85642b69e127d3a02429019336c9546f14cb738b33230da2f11a3e42e6910780` |

The exact good cap has `log2(Cgood)=-127.27723991793381`.
The following are FLOAT regional comparisons, with no source union factor:

| Potential `(s,f)` | `log2 Pr[bad]` | Conditional good+bad | Unconditional raw curve | Improvement |
|---|---:|---:|---:|---:|
| `(18,0)` | -23.9531190674 | -25.4554134799 | -21.6408803469 | 3.8145331330 bits |
| `(16,1)` | -24.0499005321 | -25.4554152816 | -21.6408803469 | 3.8145349347 bits |

Against the looser straight line, the gains are 12.2065381107 and
12.2071454288 bits. These are not gains against the raw curve.
The good branch dominates both raw-curve comparisons; reducing only the
bad-event probability would now have little effect at occupancy nine.

## Feasibility of a conditioned column control

The next step must retain the good/bad conditioning through the whole
column coefficient. It cannot multiply the separately fractionalized
regional diagnostic across masks. The following reduction supplies a
possible small interface; it has not been implemented or numerically tested.

Put `delta=255^-16`, `Bgood=beta+delta/beta`, and
`pgood_atom=delta/Bgood^2`. The exact packet powers obey
`g_0=1`, `g_1=0`, and `g_m=beta*g_(m-1)+delta*g_(m-2)`.
Since `Bgood^2>=beta*Bgood+delta`, induction gives

\[
g_m\le\delta B_{\rm good}^{m-2}\quad(m\ge2).
\]

Here `beta>1`, so the maximum through `m` nonempty packets is `g_m`.
Thus the good conditional cap has a positive moment envelope: order zero
is one, and order `m>0` is `pgood_atom*Bgood^m`. This is the moment sequence
of an atom at `Bgood` and an atom at zero. It admits a synthetic positive
occupancy-one term, which is slack rather than a claim of a real return.

For the bad branch, suppose the saved unconditional curve is bounded by
`pold*Bold^m` for `m>0`. Divide `pold` by the bad-event probability and check
that the result is at most one. Together with order zero equal to one,
this gives another probability moment envelope for the conditional mean.
If the two regional parities have different bad probabilities, one may use
their smaller probability for this envelope, while retaining their exact
probabilities in the subsequent flag sum.

For either moment law `j`, define

\[
F_n^{(j)}(X)=\mathbb E_{L_j}(1+255L_jX)^n.
\]

In the positive quiet-mask comparison, let `h` count the selected source
groups. For any fixed arity-eight row sets, the regional incidences satisfy
`0<=n_i<=h` and `sum_i n_i=8h`. These row sets do not change across columns.
For independent variables with the two prescribed moment laws, write
`A=1+255L_iX` and `B=1+255L_jX`. Then, for `a,b>=1`,

\[
F_{a+1}^{(i)}F_{b-1}^{(j)}+F_{a-1}^{(i)}F_{b+1}^{(j)}
-2F_a^{(i)}F_b^{(j)}
=\mathbb E[A^{a-1}B^{b-1}(A-B)^2].
\]

Every coefficient on the right is nonnegative. Therefore a transfer of
incidence between two regions is coefficientwise convex, even when their
moment laws differ. Repeated endpoint interpolation bounds the product by
a convex combination of products with eight incidences equal to `h` and
the other eight equal to zero.

For a flag vector with `b` bad regions, each coefficient is consequently
bounded by the maximum coefficient among

\[
(F_h^{(\rm good)})^{8-k}(F_h^{(\rm bad)})^k,
\qquad \max(0,b-8)\le k\le\min(8,b).
\]

This is a maximum, not a sum over new row-set assignments. The bound holds
for each fixed assignment. All nine products can be cached for each
`h=0,...,18` and reused across the 17 possible bad-region counts.
After summing the positive group-mixture weights over `h`, take the alpha
power of each complete column coefficient. Then sum the exact probabilities
of the bad-region counts. For alternating parities, those probabilities
come from two independent eight-region binomials with their own parameters.

This requires about 171 cached products of degree at most 144, not a
six-type regional bank or an enumeration of `2^16` flag vectors. It still
covers only the quiet/zero-boundary component. The good moment envelope is
weaker than the raw curve at some high occupancies, so a column improvement
does not follow merely from the occupancy-nine result.
