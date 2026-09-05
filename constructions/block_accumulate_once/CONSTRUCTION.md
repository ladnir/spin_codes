# BlockAccumulateOnce construction

## Outer block code

Fix a length (N). Let

\[
B_N\subseteq\mathbb F_2^N
\]

be a nonzero binary linear code. The definition permits (B_N) to contain a
scaled double-parity code followed by scaled binary BCH constituents. All
such layers belong to (B_N); the inner map below sees only the resulting
binary word.

For (0\le h\le N), define the binary weight profile

\[
A_h(B_N):=
\left|\{x\in B_N\setminus\{0\}:\operatorname{wt}(x)=h\}\right|.
\]

## Setup

Setup samples one permutation

\[
\pi\gets S_N.
\]

The sampled permutation is public and remains fixed for every encoding.

## Accumulator

For (u=(u_1,\ldots,u_N)\in\mathbb F_2^N), define

\[
s_0:=0,
\qquad
s_i:=s_{i-1}+u_i,
\qquad
\operatorname{Acc}(u):=(s_1,\ldots,s_N).
\]

The construction does not append a termination bit. Addition is in
(\mathbb F_2).

## Code

The setup-dependent code is

\[
\operatorname{BAO}_\pi(B_N)
:=
\{\operatorname{Acc}(\pi(x)):x\in B_N\}.
\]

The acronym `BAO` is used only in formulas. The canonical construction name
is `BlockAccumulateOnce`.

## Variant boundary

BlockAccumulateOnce has exactly one global bit permutation and one binary
accumulator. A packet permutation, parallel accumulators, a second
accumulator, or another randomized layer defines a different construction.

The outer block code may scale with (N). A constant-size direct-sum block
code is not the intended asymptotic instance.
