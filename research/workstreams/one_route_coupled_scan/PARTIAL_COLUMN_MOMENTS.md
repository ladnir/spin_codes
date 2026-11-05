# Partial columns at the smallest outer support

The full-column gates do not cover a word that vanishes in some components
of an occupied RS column. This note bounds every such pattern when one
outer group occupies exactly 33 columns. It averages the actual Toeplitz
mixers and retains the paired route through the existing linear envelopes.

The selected calculation closes this class at 35% output distance with
ideal-setup failure below 2^-219, including all 8192 group choices. It does
not cover partial columns on 34 or more occupied columns, or multiple active
groups. The construction and the certified reference are unchanged.

## 1. Why the 33-column case has a small common parameter space

Let E=GF(q), with q=256, and let F be the degree-16 extension used in
[the construction](CONSTRUCTION.md). One outer group contains four copies
of RS[64,32] over E, extended linearly to F. Each column receives an
independent uniform invertible 4-by-4 Toeplitz matrix over E.

Fix a set H of 33 occupied column positions. For each j in H, fix a
nonempty coordinate mask S_j contained in {0,1,2,3}; let k_j=|S_j|.
Let S be the corresponding 256-coordinate support, empty outside H.
We will require an F-valued message to have exact outer support S.

First impose only the zero columns outside H. Invertibility of each mixer
forces the four RS rows to vanish there before mixing. Each shortened RS
row has dimension one. There are therefore fixed nonzero scalars f_j in E
such that its four-symbol column is

    f_j v,       v in E^4,

and the scalar extension uses the same formula with v in F^4. The scalars
f_j come from the degree-31 polynomial vanishing at the other 31 evaluation
points. They are nonzero at every j in H.

For a realized outer, impose the additional scalar zeros outside S_j in
each occupied column. Write d_S for the dimension of this shortened space.
Its dimension is the same over E and F because all constraint matrices
have entries in E. In particular, 0<=d_S<=4.

The full masks S_j={0,1,2,3} always give d_S=4. Thus outer randomness cannot
discount that full-column class. Partial masks impose additional random
linear constraints on the same four amplitudes; this is where averaging
the outer can help.

## 2. Count pairs to bound the larger rank potential

The inner continuation now uses a global rank marker

    tau=q^2=65536.

This is a proof parameter, not a change to the coefficient distribution.
Independent zero equations cost 1/tau in the rank potential, while
exceptional coefficient roots still have probability 1/q. The empty-gap
Bellman operator distinguishes these costs. An ordinary support-transfer
bound used at an active boundary must use the common factor
1/min(tau,q), not 1/tau for both cases.

For d>=1,

    tau^(d-1) <= (q^(2d)-1)/(q^2-1).                         (1)

The numerator counts ordered pairs of E-valued codewords supported within
S, excluding the pair of two zero words. It does not count F-projective
messages. We do not assume that an F-valued witness uses only two E lanes.
The pair count bounds the dimension potential for the entire shortened
F-space, so the witness may still be chosen after the complete setup.
Neither pair member must satisfy the inner zero equations or represent the
eventual bad message: the pair count replaces only the outer dimension
factor, and the inner continuation is uniform over the realized shortened
space with this source mask.

For an exact support that is not realizable over F, its contribution is
zero. If d_S=0, it also contributes zero; do not retain tau^(-1).
For a realizable support, (1) permits the upper bound

    E_outer[1[exact S is realizable] tau^(d_S-1)]
        <= E_outer[q^(2d_S)-1]/(q^2-1).                    (2)

An ordered pair of shortened RS inputs is represented by (v,w) in
E^4 times E^4. Its column pair is (f_j v,f_j w). Consequently its E-rank

    t=dim_E span{v,w}

is the same at all 33 occupied columns. The numbers of nonzero ordered
pairs of the two possible ranks are exactly

    N_1=(q^4-1)(q+1),
    N_2=(q^4-1)(q^4-q).                                  (3)

For rank one, count v!=0 and w in E v, together with v=0,w!=0.
For rank two, choose v!=0 and then w outside E v. These cases partition
all q^8-1 nonzero ordered pairs.

## 3. The actual Toeplitz probability for a column subspace

Fix a t-dimensional subspace U of E^4, with t=1 or 2. A column pair
satisfies the required zeros exactly when its mixer T sends U into the
coordinate subspace allowed by S_j. The following bounds are uniform in U.
Mixers in different columns remain independent for each fixed pair (v,w).

### A product-span bound for zero coordinates

First sample all seven Toeplitz diagonals independently from E, without
conditioning on invertibility. Suppose a output coordinates must vanish
on U, where a>=1. Their linear constraints on the seven diagonals have
rank at least

    a+t-1.                                                (4)

To see this, write T_(i,j)=d_(i-j), with indices from zero to three.
Associate an input vector v with the reversed polynomial

    p_v(X)=sum_(j=0)^3 v_j X^(3-j).

The equation (T v)_i=0 has coefficient polynomial X^i p_v(X) on the seven
diagonals, indexed from degree zero through six. The constraint span is
therefore the product of the a-dimensional span of the selected monomials
X^i and the t-dimensional polynomial space {p_v:v in U}.

Choose a basis of the latter space with distinct leading degrees. The
pairwise sums of its t leading degrees and the a selected monomial degrees
contain at least a+t-1 distinct integers. Products with distinct leading
degrees are linearly independent, proving (4). This is ordinary polynomial
multiplication; there is no cyclic reduction or wraparound.

The unconditioned probability of all these zero equations is at most
q^(-(a+t-1)). The probability that a square Toeplitz matrix is invertible
is 1-1/q; see the established
[Toeplitz column-mixer lemma](../char2_field_spin/HD_GROUPED_OUTER.md#invertible-toeplitz-column-mixers).
Conditioning on invertibility thus gives the upper bound

    q^(-(a+t-1))/(1-1/q).                                 (5)

If there are no zero coordinates, the probability is one; do not apply
(4) with a=0. Invertibility also makes the event impossible whenever the
allowed coordinate subspace has dimension less than t.

### Constants for one- and two-dimensional input subspaces

For t=1, the established atom bound is

    rho_T=1/((q-1)q^3).

An allowed coordinate subspace of dimension k contains q^k-1 nonzero
vectors. Since T is invertible, a uniform upper bound is

    eta_k=min(1,(q^k-1)rho_T).

For t=2, apply (5) with a=4-k, together with the dimension obstruction.
The resulting exact rational constants are:

| Allowed coordinates k | Rank-one bound eta_k | Rank-two bound theta_k |
|---|---|---|
| 1 | 1/q^3 | 0 |
| 2 | (q+1)/q^3 | 1/((q-1)q^2) |
| 3 | (q^2+q+1)/q^3 | 1/((q-1)q) |
| 4 | 1 | 1 |

The events allow further zeros inside S_j. This is intentional containment
counting, not a claim that either constant is an exact-support probability.
The bounds use the actual small-coefficient Toeplitz distribution, not a
uniform matrix over F or over all invertible 4-by-4 matrices.

For each fixed pair (v,w), multiply the appropriate column bounds using
independence. Equations (2) and (3) then give the source bound

    E_outer[1[exact S is realizable] tau^(d_S-1)]
        <= c_1 product_(j in H) eta_(k_j)
         + c_2 product_(j in H) theta_(k_j),               (6)

where the constants can be stored as integers:

    c_1=N_1/(q^2-1)=q^3+q^2+q+1,
    c_2=N_2/(q^2-1)=(q^2+1)(q^4-q).

As a check, if every mask is full, their sum is
q^6+q^4+q^2+1. This is slightly larger than the direct full-column source
factor tau^3=q^6, as expected from the projective-count bound (1).

## 4. Keep the scalar-region order in the finite sum

Use common positive envelopes from
[the empty-span argument](EMPTY_SPAN_RANK.md#6-a-two-state-envelope-for-the-full-support-one-group-gate),
verified for the same z and tau. A region contains m=2048 inner steps.
Let B_0 be the envelope when its source coordinate is zero and B_1 the
envelope when that coordinate is nonzero:

    B_0=diag(1,lambda^m),       B_1=B_env.

For one active group, a nonzero source coordinate is placed uniformly among
the region's 8192 positions by the actual pair permutation and fair swap.
B_1 averages these positions only after a common linear envelope has been
proved for every fixed route. Regional positions are independent.

For a four-coordinate mask M, retain the physical scalar-region order:

    C_M=B_(1[0 in M]) B_(1[1 in M])
        B_(1[2 in M]) B_(1[3 in M]).

The B matrices need not commute. Define three nonnegative matrices:

    C_0=B_0^4,
    C_1=sum_(M nonempty) eta_(|M|) C_M,
    C_2=sum_(M nonempty) theta_(|M|) C_M.                    (7)

Let ell=(1,0) select the initial quiet state. If v_min is the minimum
verified live potential, let t_free=(1,1/v_min)^T. This terminal vector
includes every final state; it does not require a reset.

For N=2^21 and e=floor(7N/20), the probability that some nonzero message
uses one group, occupies exactly 33 outer columns, and has output weight
at most e is bounded by

    U_33 = 8192 z^(-e) ell {
               c_1 [u^33](C_0+u C_1)^64
             + c_2 [u^33](C_0+u C_2)^64
           } t_free.                                    (8)

The coefficient extraction counts the choice of 33 columns once. The sums
in (7) count every nonempty scalar mask at each selected column once. No
additional binomial coefficient, factor of 15^33, or tau^3 belongs outside
(8). The initial rank factor has already been replaced by (6).

Equation (8) follows by first fixing an exact source support and the outer.
The common inner envelope is uniform over its actual input forms. Average
the route and inner coefficients, multiply by the source factor, and then
apply (6) to average the outer. Finally sum supports in chronological column
order. This order of argument preserves the setup-conditioned existential
message quantifier; it does not assume a fresh message-independent state.

Each coefficient extraction is a positive two-state dynamic program. For
t=1 or 2, initialize row vectors p_0=(1,0), p_d=0 for d>0. At each of the
64 columns update, for 0<=d<=33,

    p'_d=p_d C_0 + 1[d>0] p_(d-1) C_t.

The final contribution is p_33 t_free. All additions and multiplications
are nonnegative, permitting directed upward arithmetic.

## 5. Checked finite gate and its scope

The root agent evaluated [the partial-column checker](theory/check_partial_columns.py)
serially, using

    tau=65536,       z=999/1000,       lambda=4981/5000.

All 32 live-state envelope inequalities pass with exact rational arithmetic.
The positive-vector proposal took 48 bounded iterations; its floating-point
proposal is not the proof. The final calculation rounds upward at 80 and
120 decimal digits. Both precision runs give an upper bound approximately
1.08313275278292e-66 and verify the conservative statement

    U_33 < 2^(-219).

The independent finite-field checker also passes 750 GF(2) subspace
fixtures for the product-rank and conditional probability bounds, and 15
ordered-DP fixtures. These small-field checks test the implementation;
the general argument is given above.

This closes every scalar-mask pattern on exactly 33 occupied columns of
one active group, including all group and column choices. It uses the
original sampling distributions and a free terminal state.

The previous marker z=.998, lambda=.993 gave an uninformative bound about
6.45e230 for this larger support class. Changing the proof marker closed
the gate; no encoder change or inferred bad-code witness was involved.

There is still no whole-code certificate. In particular, the proof above
uses the constant column rank of a pair after shortening to 33 columns.
For larger occupied-column sets, a pair's rank can vary between columns;
the same product formulas must not be reused without a new bound.

## 6. A dimension cap for the remaining column classes

The following deterministic bound may reduce the cost of a later bridge.
Suppose an exact support occupies h>=33 columns, with allowed sizes k_j.
After imposing the other zero columns, each RS row has dimension

    b=h-32.

Any b occupied columns form an information set for the shortened RS rows.
Their invertible mixers preserve that property. Restricting their output
values to the selected scalar coordinates therefore gives

    d_S <= d_max := min_(J subset H, |J|=b) sum_(j in J) k_j.

Thus d_max is the sum of the b smallest allowed column sizes. Equivalently,
the zero-form rank is at least the sum of the b largest per-column hole
counts. This strengthens the crude cap 4(h-32) for partial columns and
holds for every allowed outer realization.

For tau>=q and d_S>=1, a safe scalar bridge is

    tau^(d_S-1)
        <= (tau/q)^(d_max-1) (q^(d_S)-1)/(q-1).

For 1<=tau<=q, omit the factor (tau/q)^(d_max-1) and retain the established
bridge. The scalar bridge may still be loose: the pair method above avoids
that dimension-cap loss for h=33. Additional random rank gains from holes
outside a chosen information set remain a separate proof opportunity.
