# Goal 03: joint BCH spectrum under the alpha schedule

## Question

A one-data-symbol outer word at position \(i\) contains the field values

\[
(x,x,\alpha_i x),
\qquad
\alpha_i=\gamma^i.
\]

Goal 02 bounded the accumulator failure probability for BCH weights
\((22,22,22)\). This exploration asks whether multiplication by \(\alpha_i\)
makes the third weight 22 unlikely.

## Lag formulation

For \(t\in\mathbb Z/(2^{64}-1)\mathbb Z\), define

\[
w(t):=\operatorname{wt}(B(\gamma^t)).
\]

Multiplication by \(\alpha_i\) shifts the exponent by \(i\). The required
joint spectrum is therefore

\[
C_i(h,h')
:=|\{t:w(t)=h,\ w(t+i)=h'\}|.
\]

The complete sequence has \(2^{64}-1\) entries and cannot be materialized.
The exploration uses bulk sampling and an exact authenticated tail.

## Bulk measurement

The bulk experiment sampled one million independent pairs \((x,i)\). The
field value \(x\) was uniform and nonzero. The index \(i\) was uniform in
\([0,16{,}384)\).

The observed BCH weights had means 64.00799 and 64.004954. Their Pearson
correlation was \(-0.0004265\). The mean absolute difference was 6.353316.

The low-threshold counts were:

| Threshold | Source count | Target count | Joint count |
|---:|---:|---:|---:|
| 36 | 1 | 2 | 0 |
| 40 | 16 | 21 | 0 |
| 44 | 326 | 361 | 0 |
| 48 | 3,819 | 3,808 | 15 |
| 52 | 25,446 | 25,191 | 621 |

At thresholds 48 and 52, the joint counts were 1.03 and 0.97 times the
products of the sampled marginals. Thus, a uniform random lag behaves like an
independent permutation at the resolution of this experiment.

This conclusion does not apply uniformly to every lag.

## Exceptional initial lags

A second experiment fixed each lag and sampled 100,000 nonzero field values.
The first coefficients were strongly correlated:

| Lag | Weight correlation | Joint-to-independent ratio at weight at most 52 |
|---:|---:|---:|
| 0 | 1.0000 | 39.53 |
| 1 | 0.5278 | 23.90 |
| 2 | 0.2550 | 14.32 |
| 3 | 0.1225 | 8.91 |
| 4 | 0.0467 | 5.63 |
| 5 | 0.0155 | 3.49 |
| 7 | 0.0002 | 1.72 |
| 8 | 0.0010 | 1.36 |
| 13 | \(-0.00005\) | 1.02 |
| 15 | 0.0024 | 0.96 |

Ordinary weight correlation disappears quickly. Tail correlation persists
for several additional lags and is not captured by Pearson correlation.

## Exact authenticated tail

The file `bch_g4_le13.bin` contains every BCH word on at most 13 four-bit
packets. Each record includes its authenticated 64-bit message. The list has
38,560 words, including 1,420 words of binary weight 22.

For each weight-22 source in this list, the exact scan evaluated all 16,384
current coefficient lags. It then tested whether the target also occurred in
the authenticated list with binary weight 22. The scan found 3,622 pairs.

The largest contributions were:

| Lag | Weight-22 pairs |
|---:|---:|
| 0 | 1,420 |
| 4 | 718 |
| 8 | 410 |
| 12 | 225 |
| 16 | 122 |
| 1 | 104 |
| 3 | 81 |
| 20 | 66 |

The 3,622 pairs occupy 39 lags. The largest such lag is 43. These pairs are
actual one-data-symbol inputs with BCH profile \((22,22,22)\), not sampling
artifacts. The scan remains incomplete for weight-22 words whose packet
support exceeds 13.

Consequently, the current alpha schedule does not make the minimum profile
rare enough to omit. Applying the Goal 02 bound separately to these 3,622
authenticated inputs gives the union bound

\[
3{,}622\cdot2^{-36.4472}<2^{-24.6246}.
\]

This is only a partial-family bound. It does not cover other minimum-weight
words or other BCH profiles.

## Shifted geometric schedule

Consider the modified schedule

\[
\alpha_i':=\gamma^{64+i}.
\]

The 16,384 exponents now range from 64 through 16,447. The exact tail scan
found no authenticated weight-22 pair at exponents 64 through 16,383. A
separate exact scan found none at exponents 16,384 through 16,447. Thus the
shifted window contains no pair from the authenticated support-at-most-13
tail.

The shift preserves the streaming implementation. Compute

\[
q:=\sum_i\gamma^i x_i
\]

with the existing recurrence, then set \(p_1:=\gamma^{64}q\). The added work
is one fixed multiplication per complete outer word. It can be implemented
with 64 shift-and-reduce steps. Distinct nonzero coefficients and outer block
distance three remain unchanged.

This schedule defines a different construction. The current candidate still
uses exponents starting at zero.

## Correlation-free proof fallback

The exact joint spectrum may be unnecessary if the accumulator bound admits a
product envelope. Let

\[
a_x:=r(\operatorname{wt}(B(x)))\ge0.
\]

Multiplication by a nonzero \(\alpha\) permutes the nonzero field elements.
Hölder's inequality gives

\[
\sum_{x\ne0}a_x^2a_{\alpha x}
\le
\left(\sum_{x\ne0}a_x^3\right)^{2/3}
\left(\sum_{x\ne0}a_{\alpha x}^3\right)^{1/3}
=\sum_{x\ne0}a_x^3.
\]

Therefore, a bound of the form

\[
Q(h,h,h')\le r(h)^2r(h')
\]

would reduce the one-data-symbol shell to the ordinary BCH weight enumerator.
This reduction is proof-valid for every deterministic coefficient schedule.
Goal 02 has not yet established such a uniform product envelope.

## Conclusion

The bulk alpha schedule mixes ordinary BCH weights well, but its initial
coefficients preserve many dangerous low-weight words. The current
\((22,22,22)\) analysis is necessary. A shifted coefficient window is a
low-cost construction alternative, while the Hölder reduction is the cleanest
route that preserves the present construction.
