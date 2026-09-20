# Goal 16: derive the packet-weight accumulator formula

## Result

The parallel accumulator has an exact four-count interface. For
(k\in\{1,2,3,4\}), let (h_k) be the number of input packets with Hamming
weight (k). Conditional on these four counts, the inner distribution does
not depend on the BCH-block weights or on the original block identities.

The exact input-output enumerator is a coefficient of a (5\times5) transfer
matrix. The five states are the possible Hamming weights of the four-bit
accumulator state. At packet width one, the formula reduces to the classical
accumulator enumerator.

This theorem replaces the BCH weight profile as the natural interface to the
inner code. It does not yet perform the outer-code sum.

## Experiment and notation

Fix a packet width (g\ge1) and a packet count (N\ge1). Let

\[
u_1,\ldots,u_N\in\mathbb F_2^g
\]

be the input packets. Define the accumulator states and outputs by

\[
q_0:=0,
\qquad
q_t:=q_{t-1}+u_t,
\qquad
y_t:=q_t.
\tag{1}
\]

For (k\in\{0,\ldots,g\}), define

\[
h_k:=\bigl|\{t:\operatorname{wt}(u_t)=k\}\bigr|.
\tag{2}
\]

The histogram satisfies (sum_k h_k=N). The total output weight is

\[
w:=\sum_{t=1}^N\operatorname{wt}(q_t).
\tag{3}
\]

For Riffle (g=4), the active interface is
((h_1,h_2,h_3,h_4)). The value (h_0) follows from (N).

## Five-state transition matrix

Suppose the old state has weight (a), and the new state has weight (b).
If their supports overlap in (i) coordinates, the input packet has weight

\[
k=a+b-2i.
\tag{4}
\]

For fixed (a,b,k), the number of new states that satisfy (4) is

\[
c_g(a,b,k)
=
\binom{a}{i}\binom{g-a}{b-i},
\qquad
i=\frac{a+b-k}{2}.
\tag{5}
\]

The value in (5) is zero when (i) is not an integer in the valid range.

Introduce variables (x_0,\ldots,x_g,z). Define the
((g+1)\times(g+1)) matrix

\[
M_g(\mathbf x,z)_{a,b}
:=
z^b\sum_{k=0}^g c_g(a,b,k)x_k.
\tag{6}
\]

The variable (x_k) records an input packet of weight (k). The factor
(z^b) records the output weight of the new state.

## Exact enumerator theorem

Let (e_0) be the unit vector indexed by state weight zero. Let
(\mathbf 1) be the all-one column vector.

**Theorem.** The number of input sequences with histogram
(\mathbf h=(h_0,\ldots,h_g)) and output weight (w) is

\[
A_{N,g}(\mathbf h,w)
=
[x_0^{h_0}\cdots x_g^{h_g}z^w]
e_0^{\mathsf T}M_g(\mathbf x,z)^N\mathbf 1.
\tag{7}
\]

**Proof.** Fix one transition from a state of weight (a) to a state of
weight (b). Equation (5) counts the new states for each possible input
weight (k). Thus the matrix entry in (6) counts every concrete transition
once.

Every input sequence determines one state path through (1). Conversely, the
state path determines each input packet as (u_t=q_{t-1}+q_t). This gives a
bijection between input sequences and state paths.

Multiplying the transition monomials records the packet histogram and total
output weight. Summing the final state proves (7). ∎

The number of input sequences with histogram (\mathbf h) is

\[
T_{N,g}(\mathbf h)
=
\frac{N!}{\prod_{k=0}^g h_k!}
\prod_{k=0}^g\binom{g}{k}^{h_k}.
\tag{8}
\]

The first factor chooses the packet-weight sequence. The second factor
chooses each packet value within its weight class. Therefore

\[
\Pr[\operatorname{wt}(y)\le D\mid\mathbf h]
=
\frac{\sum_{w=0}^D A_{N,g}(\mathbf h,w)}
{T_{N,g}(\mathbf h)}.
\tag{9}
\]

The probability in (9) is over a uniform input sequence with the prescribed
packet-weight histogram.

## Why the Riffle inner code has this conditional law

Fix the BCH output in each active block. The encoder independently permutes
the 128 bits of every block. Conditional on the weight of every resulting
packet, each packet support is uniform within its weight class. These packet
supports are conditionally independent.

The encoder then applies one uniform permutation to all packet positions.
Conditional on the global histogram (\mathbf h), the final packet-weight
sequence is uniform. Each packet value remains uniform within its weight
class. Hence the final packet sequence is uniform over the set counted in
(8). Equation (9) therefore gives the exact conditional law of the Riffle
parallel accumulator.

The common packet permutation still correlates the four accumulator lanes.
The transition matrix retains this correlation through the accumulator-state
weight. No independent-lane assumption is used.

## Reduction to the scalar formula

For (g=1), the transition matrix has two states. Let (h:=h_1). For
(h>0) and (w>0), equation (7) reduces to

\[
A_{N,1}(h,w)
=
\binom{N-w}{\lfloor h/2\rfloor}
\binom{w-1}{\lceil h/2\rceil-1}.
\tag{10}
\]

To obtain (10), split the output into runs of ones and zeros. The path starts
at zero and changes state exactly (h) times. It therefore has
(\lceil h/2\rceil) positive one-runs. The two binomial factors count the
positive one-runs and the compatible zero gaps. For (h=0), only (w=0)
has count one.

Thus (7) is a direct generalization of the standard accumulator formula.

## Packetization of one BCH block

One permuted BCH block contains 32 four-bit packets. For a packet histogram
(\mathbf h=(h_0,\ldots,h_4)), define

\[
|\mathbf h|:=\sum_{k=0}^4h_k,
\qquad
\|\mathbf h\|:=\sum_{k=0}^4kh_k.
\]

If (|\mathbf h|=32) and (\|\mathbf h\|=r), the number of weight-(r)
binary blocks with histogram (\mathbf h) is

\[
P_r(\mathbf h)
=
\frac{32!}{\prod_{k=0}^4h_k!}
\prod_{k=0}^4\binom4k^{h_k}.
\tag{11}
\]

Equivalently, these counts are the coefficients of

\[
\left(
t_0+4t_1x+6t_2x^2+4t_3x^3+t_4x^4
\right)^{32}.
\tag{12}
\]

Dividing (11) by (inom{128}{r}) gives the packet-histogram distribution
of a uniformly permuted BCH block of weight (r). Independent blocks compose
by multiplying their generating functions. The global packet permutation
then passes only the summed histogram to (9).

## Verification

The audit performs five exact checks.

1. It compares all 125 values of (c_4(a,b,k)) with direct four-bit state
   enumeration.
2. It compares the five-state and 16-state weighted transfer operators on
   eight integer-weight instances.
3. It compares every scalar coefficient through length 20 with (10).
4. It enumerates all `1,048,576` four-bit input sequences of length five.
   All 126 packet histograms and every output coefficient match (7).
5. It enumerates all 58,905 packet histograms of one 128-bit block. Their
   masses reproduce (inom{128}{r}) for every (r\in\{0,\ldots,128\}).

The fifth check includes the exact weight-22 and weight-24 BCH slices used in
Goals 11--15.

## Consequence for the proof strategy

The inner code should no longer receive an ordered BCH weight profile. It
should receive one packet histogram (\mathbf h). The outer and local
randomization stages must produce a bound on the distribution of
(\mathbf h).

A useful next interface is a positive rank-one envelope

\[
\Pr[\operatorname{wt}(y)\le D\mid\mathbf h]
\le
C\prod_{k=1}^4\tau_k^{h_k}.
\tag{13}
\]

Equation (12) makes the expectation of the right-hand side explicit for each
BCH weight. Hölder's inequality can then sum the block factors over the outer
field code before any BCH weight profile is selected.

This route would replace the current profile ladder by four optimized packet
tilts. Goal 15 supplies an exact occupation-three benchmark for testing the
loss in (13).

## Reproduction

```powershell
python scripts/verify_riffle_parallelacc_packet_weight_formula.py
```

The receipt is
`receipts/goal16_packet_weight_formula_audit.json`.
