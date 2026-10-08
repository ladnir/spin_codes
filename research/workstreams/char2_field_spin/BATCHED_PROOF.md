# A proof-oriented batched field inner

The two-scan construction in [SMALL_COEFFICIENTS.md](SMALL_COEFFICIENTS.md)
has complete finite certificates, but two global permutations are expensive.
This note studies one regional permutation followed by a batched inner.
The local statements below are proved or exhaustively checked as specified.
**There is no whole-code certificate for either batched candidate yet.**

## The implemented native-field candidate

Fix the consumer's field F = GF(2^m), in the declared polynomial basis.
The current prototype uses m=128 and repeated GRS[8,4] outer blocks.
There are L=K/4 outer blocks. Outer coordinate r goes to region r;
each of the eight regions independently permutes its L individual symbols.
The region-major stream is partitioned into steps of t=8 or t=16 symbols.
The state persists across region boundaries.

In a step, X is a t-symbol input and Q is a two-symbol incoming state.
Let A have rows A_j=(1,j), for polynomial-basis labels j=0,...,t-1,
and put C=A^T. The original prototype computes

    Y = X + A Q,          Q_next = M Q + C X.

Each M is the product of three independent shears described in
[the arithmetic note](lane_mixing/README.md).
Their coefficients are sampled from 128 distinct nonzero native-field
polynomials of degree at most seven. They are not sampled from GF(256)
when m>8. Optional independent coefficients from the same set multiply
each outer output immediately before routing and the inner.

The implementation also supports pre-output refresh:

    R = M Q,             Y = X + A R,          Q_next = R + C X.

Pre-output refresh lets fresh state randomness control the current output
constraints. It uses the same number of state multiplications.
Neither refresh placement, with or without the optional diagonal,
currently has a finite distance certificate.

For one fixed input word, fresh diagonal coefficients also randomize the
feedback sum. More generally, let the active input forms f_j be fixed
in a quotient by earlier zero equations. The vectors C_j tensor f_j
have rank at least min(2,w), where w is the number of nonzero forms.
Any two distinct C_j are independent. Independent coefficients from the
128-element set therefore give a specified feedback target probability
at most 128^(-min(2,w)). This local pivot bound does not yet account for
the output constraints or the complete code's source multiplicities.

## What the cheap maps imply

For either width, C A=0. Indeed, t, the sum of the labels, and the
sum of their squares all vanish in characteristic two. Consequently
C Y=C X in both recurrences above. Starting from state zero, a second
application of the same inner follows the same states and cancels its
first output contribution. Thus the inner, without its input diagonal,
is an involution. This identity alone does not violate a 10% target.

There is a general upper bound from inputs that never activate the state.
Consider any s-symbol state, t-symbol steps, and repeated [b,k] linear
outer blocks. Let T=N/t, assuming complete steps. Choose

    a = floor(s(T-1)/k) + 1

outer blocks, provided that many are available. Their message subspace
has dimension ak. Impose C X_i=0 for the first T-1 steps. These impose
at most s(T-1) linear equations, so a nonzero message remains.
The state stays zero through every emitted step, including the last.
The final state is discarded and need not be constrained.
The resulting output equals its outer word and has weight at most ba.
At rate one half and s=2, the resulting relative-distance ceilings are
approximately 50% for t=8 and 25% for t=16. Width t=64 would impose
an approximately 6.25% ceiling and cannot meet the current 10% target.

The intended route must permute individual symbols, not intact cohorts
of t outer blocks. Under an intact-cohort route, one cohort has 4t
message dimensions but at most 8s feedback equations. For t>=8 and s=2,
a nonzero transparent word then exists within only 8t output symbols.
Independent regional permutations of individual symbols avoid that
deterministic dimension-count obstruction; they do not by themselves
complete the distance proof.

## An exactly analyzable GF(256) alternative

The next proof-oriented candidate uses all outer, diagonal, and inner
coefficients in F_q, with q=256. This is a different ensemble from
the native-field small-coefficient prototype. Over GF(2^m), for m
divisible by eight, use an explicit F_q basis and scalar extension.
The entire generator, including the GRS evaluation points, must lie
in that embedded subfield. Bytewise operations in an unrelated native
polynomial basis are not an implicit substitute.

Scalar extension preserves minimum symbol distance exactly. Each
component in an F_q basis is a base-field codeword whose support lies
inside the extension-field word's support. Conversely, base-field
codewords are also extension-field codewords. This alternative can
therefore certify every listed m through its F_q code, but cannot
gain the larger field's asymptotic distance merely by increasing m.

For each input coordinate, sample an independent diagonal coefficient
uniformly from F_q*. For each step, independently sample M uniformly
from SL(2,F_q), and use the pre-output recurrence

    R = M Q,             Y = X + A R,          Q_next = R + C X.

Here X includes the fresh coordinatewise diagonal. All these setup
choices are independent of the route and then fixed for encoding.
The analysis below fixes one base-field source word before averaging over setup.
Its nonzero X coordinates are independent uniform members of F_q*.
If Q is nonzero, R is uniform over F_q^2 minus zero, independently
of the current X. If Q is zero, R is zero. Future refreshes make the
two state categories sufficient even after conditioning on earlier
output weights.

Uniform SL(2,F_q), rather than GL(2,F_q), suffices because it acts
transitively on nonzero vectors. It has an exact three-shear sampler.
With probability q/(q+1), sample a,c from F_q and b from F_q*, then
use U_c L_b U_a. These are precisely the determinant-one matrices
whose lower-left entry is nonzero. Otherwise sample d from F_q*
and e from F_q, and use [[d,e],[0,d^-1]]. The two disjoint cases
cover SL(2,F_q), with each matrix equally likely. This sampler has
not been substituted into the measured native-field prototype.

## Exact local transfer, including the joint return event

Fix one step's active lane set S, put w=|S|, and take 0<=z<=1.
The transfer entry from category a to category b is the expectation
of z raised to the output weight, restricted to the next category b.
Rows and columns are ordered zero state, nonzero state.

Write p=1/(q-1) and f(z)=p+(1-p)z. For a nonzero incoming state,
the exact output moment is

    G_w(z) = (1-t/(q+1)) z^(t-w) f(z)^w
           + w/(q+1) z^(t-w+1) f(z)^(w-1)
           + (t-w)/(q+1) z^(t-w-1) f(z)^w.

Omit the second term when w=0 and the third when w=t.
A nonzero R gives A R at most one zero coordinate. Each prescribed
lane is that zero with probability 1/(q+1); no lane is zero with
the remaining probability. At an active lane where A R is nonzero,
the diagonal-randomized input cancels it with probability p.

For a zero incoming state, Y=X. Let p0(w) denote Pr[C X=0]. Then
p0(0)=1, and for w>=1,

    p0(w) = ((q-1)^w + (-1)^w (q^2-1-w q(q-1)))
            / (q^2 (q-1)^w).

In particular, p0(1)=p0(2)=0. This formula counts nonzero solutions
to two linear equations with pairwise independent columns.

Define J_S(z)=E[1[Q_next=0] z^wt(Y)] for a nonzero incoming state.
The exact transfer is

    W_S(z) = [[z^w p0(w), z^w(1-p0(w))],
              [J_S(z),   G_w(z)-J_S(z)]].

The joint term J_S must not be replaced by a product of marginals.
It can depend on the lanes, not only w. For S={0},
J_S=z^(t-1)/(q^2-1); for S={1}, J_S=z^t/(q^2-1).
On return, R=C_j X_j. Thus output lane j equals j^2 X_j,
and another lane l equals (1+lj)X_j, proving these two examples.

There is a finite exact formula for every J_S. For a set U of lane
labels and a target h in F_q^2, let N_U(h) count nonzero assignments
x_j satisfying sum over j in U of C_j x_j=h. For n=|U|>=1,

    N_U(0) = ((q-1)^n + (-1)^n(q^2-1-nq(q-1))) / q^2,

    N_U(h) = ((q-1)^n + (-1)^n(-1+nq-e q^2)) / q^2,  h != 0,

where e is one if h lies on a line spanned by some C_j for j in U,
and zero otherwise. For U empty, N_U(h)=1[h=0].
To prove the formulas, Fourier-expand the two equations. Each
nonzero dual vector annihilates at most one C_j. The character sum
over one dual line is q-1 when h is on its perpendicular column
line, and -1 otherwise. Summing these disjoint lines gives the
displayed counts.

Take one representative r from each of the q+1 projective state
directions. Put v_j=A_j r, and let d_S(r) count nonzero v_j outside S.
For T contained in S, require v_j!=0 for every j in T and define

    h_T(r) = r + sum over j in T of C_j v_j.

Then

    J_S(z) = 1/((q+1)(q-1)^w)
             * sum over projective r
               sum over eligible T contained in S
                 z^(d_S(r)+w-|T|) (1-z)^|T| N_(S\T)(h_T(r)).

The subset expansion replaces each active-lane factor by
z+(1-z)1[X_j=v_j]. Its return constraint is C X=r. Every term is
nonnegative for 0<=z<=1. Homogeneous scaling accounts for the
q-1 nonzero states in each projective direction. At z=1 the formula
reduces to J_S(1)=(1-p0(w))/(q^2-1), as required.

## Regional routing and the remaining whole-code work

Although W_S can depend on S, uniform regional routing permits exact
averaging over masks. Fix a source word with h_r active coordinates
in region r. That region's active positions are a uniform h_r-subset
of its L positions. If t divides L, its transfer is therefore

    [u^h_r] (sum over S contained in {0,...,t-1}
                  u^|S| W_S(z))^(L/t) / binom(L,h_r).

Conditional on step occupancies, the lane masks are independent and
uniform within their prescribed sizes. Independent input diagonals
remove the original symbols' value correlations. The matrix product
must preserve state across regions. Its normalizing factor is the
product of the regional binomials, not binom(N,sum h_r).
Unaligned natural sizes need explicit handling of boundary steps;
the displayed factorization assumes aligned regions.

A whole-code certificate still needs the repeated-outer multiplicities
and complete coverage of the regional occupation vectors. An exact
two-state local transfer is only one component of that proof.

## One-active-outer-row screen

For one active outer row, each selected outer coordinate creates
exactly one active symbol in its region. The regional matrix has a
closed geometric form, so this test needs no large polynomial table.

Let H=L/t, let a=G_0(z), and let

    F_H(a) = (1-a^H)/(H(1-a)), with F_H(1)=1.

An empty region has matrix diag(1,a^H). An occupied region has matrix

    [[0,                 z F_H(a)],
     [J_bar_1(z) F_H(a), (G_1(z)-J_bar_1(z)) a^(H-1)]],

where J_bar_1 averages the t singleton masks. Multiplying these
matrices for the outer support and summing the exact MDS support
multiplicities gives the expected count of bad projective source
lines after a Chernoff bound. The calculation divides word counts
by q-1 because scalar multiples have identical output support.

The following results use floating-point arithmetic and a scalar
tilt search. They are **Q1 screens, not outward-rounded certificates**.
The target is K=2^16, rate one half, distance at least ceil(0.1 N).

| Outer | Width t | Q1 Chernoff margin | Minimum-weight contribution at the selected tilt |
|---|---:|---:|---:|
| GRS[8,4] | 8 | 19.880 bits | 19.970 bits, weight 5 |
| GRS[8,4] | 16 | 19.877 bits | 19.967 bits, weight 5 |
| GRS[16,8] | 16 | 49.791 bits | 49.939 bits, weight 9 |

There is also an exact obstruction to obtaining 40 bits from the
Q1 first moment for GRS[8,4]. Fix the minimum-weight outer support
{0,1,2,3,7}, which contains one projective outer codeword. Require
its region-0 symbol in the last ell steps, and its region-1 symbol
in the first ell steps. Impose the same conditions on regions 2
and 3. Put its region-7 symbol in the last ell steps. Require the
state to return to zero at the symbols in regions 1 and 3.

The output occupies at most 5 ell steps, so its weight is at most
5 ell t. Each return has probability 1/(q^2-1), independently of
the earlier history on this event. Its probability for each outer
row is exactly (ell/H)^5/(q^2-1)^2. At t=16, H=1024, ell=163,
the output weight is at most 13040, below 10% of N=131072.
Summing over L=16384 rows lower-bounds the Q1 expected bad-line
count by L*(163/1024)^5/(256^2-1)^2, approximately 2^-31.256.
Thus no refinement of this first moment can reach a 40-bit bound
for this ensemble. This expectation lower bound is not a lower
bound of the same size on the probability that a sampled code fails;
different source-line events can overlap.

The next useful experiment is therefore GRS[16,8] with the exact
GF(256) refresh distribution. Its Q1 screen leaves roughly ten bits
beyond the target. The unresolved task is to cover every higher
outer-row count without losing that allowance. This is a proposal,
not a promoted design or a replacement for the certified two-scan code.

## Reproduction and validation

The exact local implementation is `theory/batched_transfer.py`.
Three tests exhaust every GF4 support mask and all nonzero inputs,
check the complete joint output-weight distribution, and check C A=0.
They also verify the exact nonzero-solution counts for every GF4
syndrome. The arithmetic is rational; no floating-point comparison
is used in these tests.

An independent derivation and six additional checks are recorded in
[the joint-transfer audit](lane_mixing/GF256_JOINT_AUDIT.md).

```sh
python -B -m unittest discover -s research/workstreams/char2_field_spin/theory -p test_batched_transfer.py -v
python -B research/workstreams/char2_field_spin/theory/q1_batched_screen.py --width 8
python -B research/workstreams/char2_field_spin/theory/q1_batched_screen.py --width 16
python -B research/workstreams/char2_field_spin/theory/q1_batched_screen.py --outer-length 16 --outer-dimension 8
```

The existing two-scan certificate code and witnesses are unchanged.
No encoding benchmark was run for this note.
