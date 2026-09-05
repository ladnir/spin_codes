# Goal 06: amortized gap route

## Result

The whole-gap model avoids the failed 18-node recovery claim. It permits one
bounded startup loss per zero-input gap. A single 24-node distance statement
would close the current support-33 row:

\[
D(24)\ge144.
\tag{A24}
\]

Here

\[
D(n)
:=
\min_{s\in\mathbb F_2^{128}\setminus\{0\}}
\sum_{t=0}^{n-1}\operatorname{wt}(q_t(s)).
\]

No counterexample to (A24) was found. The state that refuted the 18-node
bound has 24-node weight 278.

Exact certificates also exclude weight at most 143 in three zero-anchor
slices. The anchor positions are 0, 12, and 23. These positions represent a
forward, centered, and terminal zero transient.

The evidence supports the amortized model. It does not prove (A24). The
remaining obstruction is an affine low-anchor problem with anchor weight one
through five.

## Why a 24-node bound suffices

The lifted transition (R) is invertible. Therefore, for all positive
integers (m,n),

\[
D(m+n)\ge D(m)+D(n).
\]

Indeed, fix a nonzero state (s). Then

\[
W_{m+n}(s)
=
W_m(s)+W_n(R^m s).
\]

The state (R^m s) is nonzero. Each summand is therefore bounded by the
corresponding value of (D).

The authenticated support-33 row contains at least 32,737 zero-input nodes
in at most 34 gaps. For gap lengths \(\ell_1,\ldots,\ell_k\), where
\(k\le34\),

\[
\sum_{i=1}^k
\left\lfloor\frac{\ell_i}{24}\right\rfloor
\ge
\left\lceil
\frac{32737-34\cdot23}{24}
\right\rceil
=1332.
\]

On the complement of the existing exact-return event, every complete block
starts in a nonzero lifted state. Condition (A24) would give total output
weight at least

\[
1332\cdot144
=191808
>188765.
\]

The strict margin is 3,043 bits. Thus each gap may lose its final 23 nodes
without reopening the row.

## Recovery of the 18-node counterexample

The Goal 05 counterexample has 18-node weight 97. Its next six node weights
are

\[
(23,33,33,27,40,25).
\]

Consequently, its 24-node weight is

\[
97+23+33+33+27+40+25=278.
\]

The counterexample is a startup transient. It is not evidence of a low
24-node orbit.

## Exact zero-anchor certificates

Fix an anchor position (j\in\{0,\ldots,23\}). The constraint

\[
q_j(s)=0
\]

defines a 64-dimensional linear shortening of the lifted state space. Delete
the forced zero output. The remaining observation code has length 1,472 and
dimension 64.

An exact linear-matroid union algorithm packed disjoint information sets in
three shortenings. A systematic enumeration then excluded every shortened
word of weight at most 143.

| Anchor (j) | Information sets | Candidates | Result |
|---:|---:|---:|:---|
| 0 | 22 | 1,082,372,342 | exhausted |
| 12 | 19 | 8,415,660,131 | exhausted |
| 23 | 20 | 4,150,424,788 | exhausted |

For example, the anchor-zero shortening has 22 disjoint information sets.
If a word had weight at most 143, one of the 22 restrictions would have
weight at most six. The certificate enumerates restriction weights zero
through five in all 22 sets. It then enumerates weight six in 12 fixed sets.
If every set had weight at least six, at most 11 sets could exceed six.
Therefore one of the 12 fixed sets would have weight exactly six.

The centered and terminal certificates use the same argument with 19 and 20
information sets. Their coverage profiles and candidate counts are recorded
in the audit receipt.

The all-even-zero subcode was checked separately. Its dimension is 31. An
exact enumeration checked 4,952,832 candidates and found no word of weight at
most 143. This result excludes the direct 24-node analogue of the alternating
zero mechanism from Goal 05.

## Remaining affine obstruction

Any counterexample to (A24) has some node of weight at most five because

\[
\left\lfloor\frac{143}{24}\right\rfloor=5.
\]

The zero-anchor certificates address anchor weight zero. For a fixed nonzero
anchor value (u), the constraint

\[
q_j(s)=u
\]

defines a 64-dimensional affine coset. There are

\[
\sum_{w=1}^{5}\binom{64}{w}
=8{,}303{,}632
\]

possible nonzero anchor values. Treating these cosets independently would
discard the useful structure and cause an unnecessary enumeration explosion.

The next proof step needs a joint decoder for two sparse objects:

1. the anchor value of weight at most five; and
2. a low restriction in one information set of the anchor kernel.

A meet-in-the-middle or information-set decoding certificate may handle this
joint problem. Repeating the zero-coset certificate for each anchor value is
not practical.

## Assessment

The amortized route is conceptually better than the paired nine-node route.
It charges one startup debt per actual gap and catches the observed recovery.
The exact zero-anchor results remove the known alternating-zero obstruction
at 24 nodes.

The route is not yet a proof. Its main risk is coding-theoretic rather than a
new two-lap phenomenon: the nonzero affine anchor slices remain unresolved.
If a joint affine decoder remains impractical, the fallback should restrict
the bad-state list to gap starts reachable from the existing wrap multiplier
and global packet structure. No construction change is justified by the
current evidence.
