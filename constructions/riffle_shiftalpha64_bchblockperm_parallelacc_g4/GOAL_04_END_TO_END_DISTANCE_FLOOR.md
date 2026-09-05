# Goal 04: end-to-end finite distance floor

## Corrected result

Let (N=524{,}352) be the number of four-bit packets. The binary output length
is (L=4N=2{,}097{,}408). Every setup of **Riffle
ShiftAlpha64-BCHBlockPerm-ParallelAcc g=4** satisfies

\[
d_{\min}\ge33.
\tag{1}
\]

The statement is deterministic. Its failure probability is zero over the
within-block permutations and the global packet permutation. This finite
distance is about (1.57\mathbin{\cdot}10^{-5}L). Equation (1) does not prove
positive relative distance as the construction grows.

An earlier version of this goal reported a probabilistic distance of 26. That
bound was valid but strictly weaker than (1). It omitted the outer code's
global minimum-weight lower bound when interpreting the first-moment sum.

## Outer binary weight

The systematic field code has block distance three. To verify this statement,
fix a nonzero message and let (r) be its number of nonzero data symbols.

- If (r\ge3), the systematic part already has three nonzero field symbols.
- If (r=1), both parity symbols are nonzero.
- Suppose (r=2), at positions (i\ne j). If (p_0=0), then
  (m_i=m_j\ne0). Distinct coefficients imply
  (p_1=(\alpha_i+\alpha_j)m_i\ne0). Thus at least one parity symbol is
  nonzero.

Every nonzero field codeword therefore has at least three nonzero symbols.
Each nonzero symbol becomes a BCH codeword of weight at least 22. The binary
word entering the packet permutation consequently has weight

\[
w\ge3\cdot22=66.
\tag{2}
\]

Neither class of permutation changes (w).

## Accumulator inequality

Let (u_1,\ldots,u_N\in\mathbb F_2^4) be the permuted input packets. Define

\[
s_0:=0,
\qquad
s_t:=s_{t-1}+u_t,
\qquad
W:=\sum_{t=1}^{N}\operatorname{wt}(s_t).
\]

Because (u_t=s_{t-1}+s_t), the triangle inequality gives

\[
\begin{aligned}
w
&=\sum_{t=1}^{N}\operatorname{wt}(u_t)\\
&\le
  \sum_{t=1}^{N}\operatorname{wt}(s_{t-1})
 +\sum_{t=1}^{N}\operatorname{wt}(s_t)\\
&=2W-\operatorname{wt}(s_N)
\le2W.
\end{aligned}
\tag{3}
\]

Equations (2) and (3) imply (W\ge33), which proves (1).

## First open threshold

Distance 34 is not yet proved. The existing fixed-word contraction and
parity-dropped outer sum give, at bad-output threshold 33,

\[
\sum_{m\ne0}\Pr_{\Pi}[W(m)\le33]
\le43.401774584603615.
\tag{4}
\]

The right-hand side exceeds one, so the first-moment argument is trivial at
the first threshold not excluded by (1). Almost all of (4) comes from messages
with at least two nonzero data symbols after both parity blocks are discarded.

The earlier threshold-25 calculation remains a useful diagnostic. Its
parity-dropped first moment was (0.8845998202997373). The deterministic
argument makes the corresponding failure event impossible, so this number is
not the construction's setup-failure probability.

## Reproduction

Run:

```powershell
python scripts/analyze_riffle_shiftalpha64_end_to_end_floor.py --distance 25
python scripts/analyze_riffle_shiftalpha64_end_to_end_floor.py --distance 33
```

At threshold 25, the script now reports a zero failure probability and the
deterministic distance 33. At threshold 33, it reports the trivial first
moment in (4).

The analysis script has SHA-256
`D0C08473819ED061A52AE38B5C64EC01A985BAF471F32F127B34DADEEEEA84B6`.
The imported fixed-word script has SHA-256
`F7D68BB5E353BE3551067B15E5FFDD9A06265A2A4D768E5F4AD371F5445B01F6`.
The exact BCH spectrum has SHA-256
`26164B40CEE431CD24995887BA485563F02048D64FDEA00DC86B3D9A4897A94E`.

## Next goal

Retain at least one parity BCH block for messages with two or more nonzero
data symbols. The immediate target is a first moment below one at threshold
33, which would prove (d_{\min}\ge34) with positive setup probability.
