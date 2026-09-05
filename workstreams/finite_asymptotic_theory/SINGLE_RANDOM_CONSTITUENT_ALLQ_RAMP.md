# All-occupation ramp for one random constituent

## Question

The exact occupation-two calculation classifies two local messages by rank.
That method does not scale directly to \(Q=L=8192\). A rank-\(r\) tuple can
have many repeated linear combinations, and a table indexed by all row
weights has dimension \(Q\).

The scalable route fixes one sampled constituent before it sums messages.
For a fixed code, the choices of nonzero local messages in different rows
factorize. Shared-code rank correlations then disappear from the message
sum; they are replaced by a deterministic condition on the constituent's
weight spectrum.

## Parameters

The target remains

\[
 B=256,\quad K=128,\quad L=8192,\quad N=2^{21},
 \quad D=\lceil0.109N\rceil=228590,
\]

with RandomStepConv-M22 and the region-permuted bit transpose.

For a binary \(K\)-by-\(B\) generator \(G\), define

\[
 A_w(G):=
 \left|\{u\in\mathbb F_2^K\setminus\{0\}:
             \operatorname{wt}(uG)=w\}\right|.
\]

## Uniform spectrum envelope

Let

\[
 \mu_w:=(2^K-1)\binom Bw2^{-B}.
\]

Suppose a full-rank generator satisfies

\[
 A_w(G)\le F\mu_w
 \tag{1}
\]

for every \(1\le w\le B\). After the row-coordinate permutation, the
counting mass at each ambient word is at most

\[
 F(2^K-1)2^{-B}.
\]

For a message with \(Q\) active rows, its complete outer counting measure is
therefore bounded by \(F^Q(2^K-1)^Q\) times \(Q\) independent uniform
ambient rows. This pointwise comparison remains valid inside the causal
inner transfer. A one-dimensional recurrence can then evaluate every
occupation \(Q\in\{1,\ldots,L\}\).

For \(F=8\), pairwise independence of the random generator's nonzero row
images gives the elementary one-sample bound

\[
 \Pr[(1)\text{ and full rank}]\ge0.7655276698957902.
\]

The all-occupation diagnostic does not close. Its largest term occurs at
\(Q=L\) and has logarithm

\[
 18102.25835146599.
\]

The failure is structural. The factor \(F^L\) costs

\[
 L\log_2 8=24576
\]

bits. Removing that factor makes the same endpoint pass by
6473.741648534 bits. Thus a uniform envelope at the dense endpoint requires

\[
 F<2^{6473.741648534/8192}=1.7293761115\ldots.
 \tag{2}
\]

The factor-8 computation is a rejected proof wrapper. It is not evidence
that the sampled code has poor distance.

## Shell-sensitive constituent event

A three-band event avoids the repeated worst-shell factor. Consider a
uniform binary \(128\)-by-\(256\) generator and require all four conditions.

1. The generator has full row rank.
2. No nonzero codeword has weight in \([1,24]\cup[232,256]\).
3. Each fringe band \([25,31]\) and \([225,231]\) contains at most 128
   nonzero codewords.
4. Every central shell \(32\le w\le224\) satisfies
   \(A_w\le1.5\mu_w\).

Pairwise independence gives variance at most the mean for every shell and
every union of shells. Markov's inequality handles the excluded extremes.
Cantelli's inequality handles the two fringe totals and the central shells.
The resulting diagnostic failure terms are

| Condition | Failure upper bound |
|---|---:|
| Excluded extreme shells | \(0.0000217745000474\) |
| Low fringe total | \(0.0027658845143881\) |
| High fringe total | \(0.0027658845143881\) |
| Central shell bounds | \(0.0537802623418237\) |
| Rank failure | \(2.94\times10^{-39}\) |

Their sum is \(0.05933380587064723\). Therefore one sample satisfies the
event with probability at least

\[
 0.9406661941293528.
 \tag{3}
\]

Ten independent attempts make the elementary abort bound smaller than
\(2^{-40.75}\). This observation does not supply an efficient acceptance
test for the event.

## Central compositions

For messages whose active rows all use central-shell codewords, the factor
is \(F=1.5\). Reweighting the all-occupation receipt gives

\[
 \max_{2\le Q\le L}\log_2 U_Q=-92.6514675456099.
\]

The maximum occurs at \(Q=2\). At \(Q=L\), the logarithm is

\[
 -1681.7288426263003.
\]

Thus every pure-central occupation closes in the diagnostic. The exact
shellwise occupation-one calculation remains stronger than the uniform
central envelope.

## Remaining finite proof

The remaining messages contain fringe rows. Their total constituent mass is
at most 128 codewords per fringe, compared with approximately \(2^{128}\)
central messages. This loss in message count should pay for their less
favorable input weights, but the proof must retain the composition.

The next calculation should use four row categories:

- inactive;
- low fringe, with weights \(25\) through \(31\);
- central, with weights \(32\) through \(224\);
- high fringe, with weights \(225\) through \(231\).

For each composition, a categorical reference distribution removes the
conditioning on category positions within each region. A convex simplex
cover can certify all integer compositions without enumerating them.

This shell-sensitive composition proof is the proposed route from
\(Q=3\) through \(Q=L\). It is not yet complete. In particular, the current
artifacts do not prove a 40-bit all-occupation statement for one sampled
\([256,128]\) constituent.

## Evidence

`evaluate_single_random_constituent_allq.py` implements the rejected uniform
factor wrapper. Its receipt is
`single_random_constituent_B256_allq_factor8_s22_d109.json`.

`evaluate_single_random_constituent_three_band_event.py` evaluates the
three-band event and reweights the pure-central occupation rows. Its receipt
is `single_random_constituent_B256_three_band_event.json`.

The exact occupation-one and occupation-two artifacts remain valid. They do
not depend on the rejected uniform factor-8 wrapper.
