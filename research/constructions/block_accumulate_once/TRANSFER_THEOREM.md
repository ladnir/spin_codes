# Exact transfer theorem for BlockAccumulateOnce

## Question

The accumulator can spread a permuted low-weight word across a large output.
The required suppression depends on both the outer minimum distance and the
number of outer words at each weight. We therefore seek a theorem whose only
outer-code input is a binary weight profile.

## Accumulator input-output count

Fix integers (1\le h\le N) and (1\le w\le N). Let

\[
T_N(h,w)
:=
\binom{w-1}{\lceil h/2\rceil-1}
\binom{N-w}{\lfloor h/2\rfloor},
\tag{1}
\]

where an infeasible binomial coefficient equals zero.

The value (T_N(h,w)) counts inputs (u\in\mathbb F_2^N) for which

\[
\operatorname{wt}(u)=h
\quad\text{and}\quad
\operatorname{wt}(\operatorname{Acc}(u))=w.
\]

To see this, view the accumulator output as runs of ones. If (h=2r), the
output has (r) one-runs and ends in zero. If (h=2r-1), it has (r)
one-runs and ends in one. Composing the one-run lengths and distributing the
zero positions gives (1).

For an integer (0\le D\le N), define

\[
P_{N,D}(h)
:=
\frac{1}{\binom Nh}
\sum_{w=1}^{D}T_N(h,w).
\tag{2}
\]

Fix a word (x\in\mathbb F_2^N) of weight (h). Under
(\pi\gets S_N), the word (pi(x)) is uniform among the weight-(h)
binary words. Thus (2) is exactly

\[
\Pr_{\pi\gets S_N}
\left[
\operatorname{wt}(\operatorname{Acc}(\pi(x)))\le D
\right].
\tag{3}
\]

## Transfer theorem

**Theorem.** Fix (B_N\subseteq\mathbb F_2^N). Suppose nonnegative numbers
(\overline A_1,\ldots,\overline A_N) satisfy

\[
A_h(B_N)\le\overline A_h
\qquad
\text{for every }1\le h\le N.
\tag{4}
\]

Then

\[
\Pr_{\pi\gets S_N}
\left[
d_{\min}(\operatorname{BAO}_\pi(B_N))\le D
\right]
\le
\sum_{h=1}^{N}\overline A_h P_{N,D}(h).
\tag{5}
\]

Consequently, if the right side of (5) is at most (2^{-\lambda}), then

\[
\Pr_{\pi\gets S_N}
\left[
d_{\min}(\operatorname{BAO}_\pi(B_N))\ge D+1
\right]
\ge 1-2^{-\lambda}.
\tag{6}
\]

If the right side of (5) is less than one, at least one permutation produces
minimum distance at least (D+1).

**Proof.** For each nonzero (x\in B_N), define the event

\[
E_x:=
\left\{
\operatorname{wt}(\operatorname{Acc}(\pi(x)))\le D
\right\}.
\]

Equation (3) gives (\Pr[E_x]=P_{N,D}(h)) when
(\operatorname{wt}(x)=h). A union bound over the nonzero words of (B_N)
gives

\[
\Pr\left[\bigcup_{x\in B_N\setminus\{0\}}E_x\right]
\le
\sum_{h=1}^{N}A_h(B_N)P_{N,D}(h).
\]

Apply (4) to obtain (5). The accumulator and the permutation are invertible,
so every nonzero outer word produces a nonzero codeword. Equations (5) and
(6) follow. ∎

## Target profile condition

For a target relative distance (\delta\in(0,1)), set

\[
D:=\lfloor\delta N\rfloor.
\]

The exact sufficient condition is

\[
\boxed{
\sum_{h=1}^{N}
\overline A_h
\frac{
\displaystyle
\sum_{w=1}^{\lfloor\delta N\rfloor}
\binom{w-1}{\lceil h/2\rceil-1}
\binom{N-w}{\lfloor h/2\rfloor}
}{\binom Nh}
<1.
}
\tag{7}
\]

Equation (7), not minimum distance alone, is the interface required from the
outer block code. A logarithmic minimum distance controls the first nonzero
term. The remaining low-weight spectrum controls the multiplicity paid by
the union bound.

## Intended scaled outer family

The intended specialization builds (B_N) from scaled BCH blocks and a
scaled double-parity layer. The scaling must make the complete profile in
(7) small enough. It is not sufficient merely to increase a parameter named
"block size."

A future specialization must establish:

1. the exact scaling of the BCH length, dimension, and binary distance;
2. the scaling of the double-parity symbols and coefficient field;
3. a rigorous envelope (\overline A_h) for all relevant binary weights;
4. the largest (\delta) for which (7) closes.

No such specialization is claimed in this record.
