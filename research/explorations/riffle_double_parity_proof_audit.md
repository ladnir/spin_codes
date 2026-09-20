# Proof audit of Riffle Double-Parity g=8

## Conclusion

**Riffle Double-Parity g=8, as benchmarked with the sequential outer layout, cannot meet the current first-moment target.**
An exact family of five-block outer codewords gives

\[
\mathbb E[Z_d]\geq 2^{-23.016236846906736}.
\]

The target is $\mathbb E[Z_d]\leq2^{-40}$.
The lower bound therefore misses the target by about 16.98 bits.

Double-Parity repairs the two-block defect of Graph-24.
It does not repair the positive entropy slope caused by rate one-half blocks and eight-bit packets.

## Scope

The audit applies to the fast sequential layout used by the Double-Parity benchmark.
Each extended BCH $[128,64,22]$ word occupies 16 contiguous eight-bit packets before the global packet permutation.

The audit does not apply unchanged to a striped layout that places one BCH word across more physical groups.
Such a layout defines a different matrix and has a different cost.

Fix

\[
K=2^{20},\qquad
B=K/64=16{,}384,
\]

and let $n:=B+2=16{,}386$ be the number of outer BCH blocks.
The unpunctured binary length is

\[
N'=128n=2{,}097{,}408.
\]

Set $d:=\lfloor0.09N'\rfloor=188{,}766$.
The packet permutation has

\[
M:=N'/8=262{,}176
\]

packet positions.
The final $T:=\lfloor d/64\rfloor=2{,}949$ recursive nodes contain

\[
L:=8T=23{,}592
\]

packet positions.

## Terminal-placement event

Fix an outer word with $r$ nonzero packets.
Under a uniform packet permutation, the probability that all active packets occupy the final $L$ positions is

\[
P(r):=\frac{\binom{L}{r}}{\binom{M}{r}}.
\]

The recursive inner map processes these nodes last.
If the terminal-placement event occurs, earlier nodes receive the zero state and zero input.
The output is therefore confined to the final $T$ nodes.
Its binary weight is at most

\[
64T=188{,}736\leq d.
\]

Consequently, every outer codeword contributes at least $P(r)$ to $\mathbb E[Z_d]$.

## Double-Parity is an MDS outer code

Identify each 64-coordinate message block with

\[
\mathbb F:=\mathbb F_2[x]/(x^{64}+x^4+x^3+x+1).
\]

Fix distinct nonzero coefficients $\alpha_0,\ldots,\alpha_{B-1}\in\mathbb F$.
For data symbols $m_0,\ldots,m_{B-1}\in\mathbb F$, define

\[
p_0:=\sum_{i=0}^{B-1}m_i,
\qquad
p_1:=\sum_{i=0}^{B-1}\alpha_i m_i.
\]

The two parity equations have columns

\[
(1,\alpha_i)^{\mathsf T},\quad (1,0)^{\mathsf T},\quad (0,1)^{\mathsf T}.
\]

Every pair of columns is linearly independent.
The outer code is therefore an $[n,B,3]$ MDS code over $\mathbb F$.

Let $q:=2^{64}$.
For $3\leq w\leq n$, the number of outer codewords with block weight $w$ is

\[
A_w=
\binom{n}{w}(q-1)
\sum_{j=0}^{w-3}
(-1)^j\binom{w-1}{j}q^{w-3-j}.
\]

## Exact first-moment obstruction

Every active BCH block occupies at most 16 packets.
An outer codeword of block weight $w$ therefore has at most $16w$ active packets.
The function $P(r)$ decreases with $r$, so every such codeword has bad-event probability at least $P(16w)$.
Hence

\[
\mathbb E[Z_d]\geq A_w P(16w)
\]

for every $w\geq3$.

Exact integer and rational arithmetic gives:

| Outer block weight $w$ | $\log_2(A_wP(16w))$ |
|---:|---:|
| 3 | -63.407484138361 |
| 4 | -43.043708373439 |
| 5 | **-23.016236846907** |
| 6 | -3.266186470799 |
| 7 | 16.247074123757 |
| 8 | 35.553281649964 |

The row $w=5$ already rules out the $2^{-40}$ first-moment target.
The stronger rows show that the problem is not a narrow minimum-distance effect.

The argument uses the maximum possible packet support of each BCH block.
It therefore remains valid for every BCH weight spectrum and every choice of distinct Double-Parity coefficients.

## Comparison with Graph-24

The same terminal event exposes a separate defect in the current fast Graph-24 outer.
Choose two unpunctured data blocks and place the same nonzero message $v$ in both.
The XOR-parity block vanishes.
The 24-bit graph syndrome vanishes with probability $2^{-24}$ over the sampled graph.

Let $s(v)$ be the packet support of the BCH encoding of $v$.
For the current power-of-two length, the exact contribution of this family is

\[
\binom{16{,}128}{2}2^{-24}
\sum_{v\neq0}
\frac{\binom{23{,}592}{2s(v)}}
     {\binom{262{,}144}{2s(v)}}.
\]

Its base-two logarithm is

\[
-35.189306074384.
\]

Thus the fast sequential Graph-24 construction also cannot meet the current first-moment target.
The earlier one-block diagnostic gave approximately $-48.17$ because it omitted this larger equal-pair family.

## Low-support Double-Parity diagnostics

The MDS obstruction makes a complete low-support enumeration unnecessary for rejection.
The enumeration remains useful for understanding coefficient structure.

For the coefficients $\alpha_i=x^i$, three consecutive data positions admit parity-zero values proportional to

\[
(xt,(1+x)t,t).
\]

The BCH code contains exactly 12,587 nonzero messages with packet support at most eight.
This fact permits complete enumeration of consecutive triples through total packet support 26.

| Total packets | Codewords per consecutive triple |
|---:|---:|
| 21 | 3 |
| 22 | 0 |
| 23 | 36 |
| 24 | 688 |
| 25 | 117 |
| 26 | 5,674 |

Across all consecutive triples, this exact tail contributes approximately $2^{-57.0312}$.
Across every three-data support, exactly 49,146 codewords attain the minimum 21 packets.
Their contribution is approximately $2^{-57.3844}$.

These values initially looked favorable because they examine only block weight three.
The exact MDS spectrum shows that block weights five and above are decisive.

## Design implication

Double-Parity removes support-two outer codewords but adds only 128 global parity bits.
That fixed codimension cannot overcome the per-block entropy slope.

One active 64-bit data block provides roughly 64 bits of value entropy.
Its BCH word occupies at most 16 packets.
Placing those packets in a $9\%$ terminal region costs only about

\[
16\log_2(1/0.09)\approx55.6
\]

bits.
The choice of the active outer position adds further entropy.
As the outer block weight grows, the number of codewords increases faster than the terminal-placement probability decreases.

A larger rate-one-half BCH block does not change the asymptotic ratio.
For contiguous packets of width eight, a rate-one-half block still provides only two packets per eight message bits.

The next construction must change at least one of the following properties:

1. reduce the effective outer rate below one half;
2. place each local word across more independently permuted packets;
3. add a second global mixing stage that has no common terminal region; or
4. return to a smaller packet width.

Adding a constant number of structured parity blocks cannot repair the positive entropy slope.

## Reproduction

Run:

```powershell
python scripts\analyze_riffle_double_parity.py
```

The script authenticates the BCH packet-support mass, evaluates the exact MDS counts, and prints the terminal-event lower bounds.
Sampling options in the same script remain diagnostic and are not used by the rejection argument.
