# Two adjacent groups under the actual pair route

Two adjacent active outer groups do not produce independent singleton
insertions. In alternate regions they form one complete physical pair.
In the other regions they occupy two distinct source pairs, which can
land in the same width-four batch. This note retains both effects.

The formulas below bound two specified source classes using a verified
common envelope. The first evaluated marker does not close them. No
whole-code certificate or low-distance counterexample follows.

## 1. Common envelopes before averaging routes

Use width four and the antipodal state cycle (0,2,1,3). Fix tau>=1,
0<z<=1, and a positive live-state vector v that has passed the empty-span
inequality T_tau(v)<=lambda v. The exceptional-root probability remains
1/256. Keep the quiet state separate, and write

    E = diag(1,lambda).

As in [EMPTY_SPAN_RANK.md](EMPTY_SPAN_RANK.md), a coefficient vector
(c_q,c_l)^T represents c_q on quiet states and c_l v(r,H) on live states.
Every envelope below acts on the same global potential tau^(-dim W).

For a nonempty input mask I of size k, let M_I be the ordinary fixed-mask
matrix evaluated with charge parameter psi=min(tau,256). The outgoing
state forms belong to the span of the incoming state and the k input
forms. After adding zero equations, their residual rank is therefore at
most r+k. Define the following bounds for each fixed mask:

    a_I = z^k max_{1<=s<=k} v(s,I),

    b_I = max_{(r,H)} M_I[H,empty]/v(r,H),

    d_I = max_{(r,H)} [sum_{H' nonempty} M_I[H,H']
              max_{1<=s<=min(r+k,|H'|)} v(s,H')] / v(r,H),

    A_I = [[0,a_I],[b_I,d_I]].

The maxima range over valid live states, 1<=r<=|H|. From a quiet state,
the output and next state equal the nonzero input. Its residual rank can
be any value from one to k: earlier common zero equations can correlate
inputs from different groups. In particular, two nonzero inputs must not
be assumed to have residual rank two. This proves the quiet row. The
ordinary rank transfer and the outgoing-rank bound prove the live row.

These inequalities hold for every fixed route, outer realization, and
retained common quotient. Thus we may average the nonnegative linear
matrices over independent route choices after establishing the inequalities.
We never average a nonlinear Bellman maximum over routes.

Define three coordinate averages:

    A_1 = (A_1bit + A_2bit + A_4bit + A_8bit)/4,
    A_P = (A_3 + A_12)/2,
    A_C = (A_5 + A_6 + A_9 + A_10)/4.

Here a numeric subscript denotes a four-bit input mask. The suffix `bit`
in the first line only distinguishes the singleton-mask matrices from
their average A_1. The P masks fill one physical pair. The C masks select
one coordinate from each physical pair. A common maximum over those masks
is also safe, but the averages preserve their exact probabilities.

## 2. One full pair versus two singleton pairs

There are g=8192 groups, M=g/2=4096 physical pair slots per region, and
m=M/2=2048 width-four batches. Each batch contains two physical pairs.

Fix the cyclic adjacent group pair (i,i+1 mod g). A region has matching
parity e=r mod 2. When e=i mod 2, the two active groups form one source
pair. Its destination pair slot is uniform among M slots. Its batch is
therefore uniform among m batches, and its input mask is 3 or 12 with
equal probability. Swapping the pair changes its values, not its support.
The full-pair regional envelope is

    B_P = R(A_P),
    R(A) = (1/m) sum_{ell=0}^{m-1} E^ell A E^(m-1-ell).    (1)

In the opposite matching parity, the two groups belong to distinct source
pairs. Each pair has one nonzero coordinate. Their two destination pair
slots are sampled uniformly without replacement from the M slots.

Exactly m of the binomial(M,2) unordered slot pairs occupy one common
batch. Thus their collision probability is

    p_col = m/binomial(2m,2) = 1/(2m-1).                  (2)

Conditional on collision, the batch is uniform. The two independent swap
bits select masks 5,6,9,10 with equal probability. This gives R(A_C).

Conditional on no collision, the two distinct batch indices are uniform
among binomial(m,2) choices. Each batch's pair slot and swap bit are
independent fair choices. Consequently its singleton coordinate is uniform
among four coordinates, independently for the two batches. The source
group labels may arrive in either order; A_I is uniform over their actual
input forms, so that order introduces no further factor.

Define

    R_2(A) = 1/binomial(m,2)
      * sum_{0<=i<j<m} E^i A E^(j-i-1) A E^(m-1-j).

The separate-pair regional envelope is exactly the following mixture of
the already established linear bounds:

    B_S = p_col R(A_C) + (1-p_col) R_2(A_1).             (3)

This is not R(A_1)^2, which would duplicate the region length, and it is
not an independent-with-replacement insertion law. At K20 the collision
probability is 1/4095; dropping it removes a permitted support pattern.

## 3. Positive formulas for the insertion averages

For any A=[[0,a],[b,d]], put

    S_m = (1/m) sum_{ell=0}^{m-1} lambda^ell.

Then

    R(A) = [[0,a S_m],[b S_m,d lambda^(m-1)]].            (4)

For two distinct insertions, let L=m-2 and C=binomial(m,2). Define the
nonnegative sums

    F_L = sum_{k=0}^L (L+1-k) lambda^k,
    G_L = sum_{k=0}^L (k+1) lambda^k.

Then

    R_2(A) = [[ab F_L/C, ad G_L/C],
              [bd G_L/C, ab G_L/C+d^2 lambda^(m-2)]].   (5)

To verify (5), write the three empty gaps as l_0,l_1,l_2>=0 with sum L.
The four entries of E^l_0 A E^l_1 A E^l_2 are

    ab lambda^l_1,
    ad lambda^(l_1+l_2),
    bd lambda^(l_0+l_1),
    ab lambda^(l_0+l_2)+d^2 lambda^L.

Summing these entries over the triples gives (5). All terms are positive,
including boundary insertions and zero-length gaps.

An alternative implementation uses degree-two matrix polynomials. The
coefficient of u^2 in (E+uA)^m is the numerator of R_2(A), with chronological
matrix order preserved. A sequential recurrence is

    S'_0 = S_0 E,
    S'_1 = S_1 E + S_0 A,
    S'_2 = S_2 E + S_1 A,

starting with S_0=I and S_1=S_2=0. Truncated positive polynomial powering
computes the same coefficient. Normalize once by binomial(m,2).

## 4. Two full groups and a free terminal

The regional routes are independent. For even i, the chronological product
for 256 occupied regions is (B_P B_S)^128. For odd i it is (B_S B_P)^128.
The products need not commute, so both orders must remain.

Let q_0=(1,0), v_min=min v(r,H), and

    t_free = (1,1/v_min)^T.

This terminal vector dominates the constant-one continuation on every
state. The ordinary envelopes need no special last-step correction and
impose no final reset.

There are 8192 cyclic adjacent group pairs, with 4096 of each parity. Two
full groups have shortened message dimension D=256 over F. Therefore,
for e=floor(0.35N), N=2^21, their union has the bound

    U_full,2 <= 4096 tau^255 z^(-e) q_0
        [(B_P B_S)^128 + (B_S B_P)^128] t_free.          (6)

It is conditional on any permitted outer realization. It covers messages
whose outer support is exactly full in these two groups and zero elsewhere,
not messages on arbitrary pairs of groups or partially populated groups.

## 5. Both groups on the same complete columns

Fix h occupied RS columns shared by the two groups. Require all four
components in every selected column to be nonzero in both groups. All
other columns and groups are zero. For h>=33, each group's shortened
dimension is 4(h-32), so the total dimension is D=8(h-32).

An empty column has envelope C_0=E^(4m). For each pair parity define

    C_even = (B_P B_S)^2,
    C_odd  = (B_S B_P)^2.

Summing all choices of h columns while retaining their order gives

    U_h,2 <= 4096 tau^(8(h-32)-1) z^(-e) q_0
      [ [u^h](C_0+u C_even)^64
        + [u^h](C_0+u C_odd)^64 ] t_free.                (7)

Coefficient extraction already counts every column set. Do not multiply
by binomial(64,h) again, or square that selection count because there are
two groups. The support set is shared in this gate. At h=33 the source
factor is tau^7; at h=64, (7) agrees with (6).

For later use, different complete-column sets can also be retained. Let
B_1=R(A_1) and C_single=B_1^4. For parity p, use

    [u^h_1 v^h_2]
       (C_0 + u C_single + v C_single + uv C_p)^64.       (8)

The source dimension is 4(h_1-32)+4(h_2-32), and the two parity terms each
have 4096 group-pair choices. Formula (8) preserves overlapping column
sets; it is not the same class as the shared-set coefficient in (7).
No calculation for (8) is asserted here.

## 6. What the first failed marker can expose

The first recorded evaluation used tau=65536, z=.999, lambda=.9962.
Root's directed calculation returned about 4.199e-4 for h=33 and
6.109e843 for h=64. These are failed upper bounds, not measured failure
probabilities. The large-h source factor alone is tau^255, so changing
tau has a much larger cost there than in the eight-dimensional h=33 gate.

The envelope also admits cheaper resets than the one-group singleton
case. With antipodal wiring, take a physical-pair input I={0,1} and current
state support H={0,2}. A reset Q_next=0 has output support {1,2} and requires

    X[0]=a_0 Q[0],
    X[1]=Q[2],
    Q[0]+a_2 Q[2]=0.

The ordinary support matrix charges two fresh coefficients. It drops the
possibly new deterministic rank of X[1]+Q[2]=0. A cross-pair collision can
reset a singleton state with only one fresh charge, but its route also
pays the factor 1/4095. These are local branch patterns, not constructed
outer messages or low-distance codewords.

The two contributions can be isolated directly in the regional product.
For instance, B_P[1,0] B_S[0,1] is a paired reset followed by a separate
region's activation. It need not use a collision. Conversely,
B_P[0,1] B_S[1,0] includes the rare cross-collision reset. If A_1 has
entries a_1,b_1,d_1, then

    B_S[0,1] = p_col a_C S_m
       + (1-p_col) a_1 d_1 G_L/binomial(m,2),

    B_S[1,0] = p_col b_C S_m
       + (1-p_col) b_1 d_1 G_L/binomial(m,2).

These formulas identify candidate losses without assuming which dominates
a particular numerical run. The live-to-live product must also be checked.
If boundary reset constants dominate, possible refinements include retaining
their deterministic rank or coupling the terminal reset function more
tightly to the preceding empty gap. Neither improvement is proved by (7).

## 7. Implementation audit and scope

[check_two_groups.py](theory/check_two_groups.py) implements (7). A read-only
audit checked its mask averages, collision factor, distinct-batch normalizer,
ordered polynomial products, source exponents, parity count, free terminal,
and positive upward-rounded arithmetic. No discrepancy was found. This
audit did not run the program or independently reproduce its numerical values.

The next gate is a verified useful bound for these actual two-group classes.
They still exclude partial columns, nonadjacent pairs, and larger active
group sets. No bound here may replace the full support coverage required
for a minimum-distance certificate.

## 8. Fresh source coordinates in the first occupied columns

The shared complete-column class admits a stronger active-step bound.
It retains some deterministic rank that the ordinary support matrix drops.
Fix h occupied columns and put b=h-32. The first b occupied columns,
in encoder order, provide independent source coordinates. They need not
be the first b physical columns: empty columns may occur between them.

To prove this fact, first fix both groups' outer setups and the full route.
Let V be the message space after imposing zeros at all absent columns.
The column mixers are invertible. Thus an absent mixed column is zero
if and only if all four constituent RS outputs there are zero.
An RS polynomial of degree less than 32 must therefore vanish at 64-h
specified evaluation points. Write it as the product of those linear
factors and a polynomial of degree less than b. Evaluation at any b
occupied columns is an invertible map on this remaining polynomial space.
Nonzero RS column multipliers do not change invertibility.

For one group, the four constituent evaluations in any b occupied columns
are consequently independent coordinates on its 4b-dimensional shortened
space. Applying an invertible four-by-four mixer separately in each
column preserves that property. The two groups have disjoint message
coordinates, so their first b occupied columns provide 8b independent
scalar forms on V. This conclusion holds over F, not just GF(256).

Consider a batch within these first b occupied columns. Let S_old be the
span of source forms consumed before this batch. Let W be the span of
the ordinary zero equations imposed so far. Both the incoming state
forms and W belong to S_old. Required nonzero forms are retained as
requirements; they are not additional equations in W. The active input
forms x_j, j in I, are mutually independent modulo S_old.

For an active coordinate j, the two possible zero equations are

    y_j = x_j + a_j Q_j = 0,
    Q'_j = x_j + a_j Q_j + (RQ)_j = 0.

If at least one equation is required, it has the unique fresh pivot x_j.
It increases the common rank by one for every value of a_j. If both
equations are required, their difference is (RQ)_j=0; there is still
only one fresh pivot. The support constraints already require the
corresponding rotated state coordinate to be zero in this case.

For a necessary support triple (I,H,J,H'), define

    c_fresh = |I intersect complement(J intersect H')|,
    c_old   = |(H minus I) intersect complement(J intersect H')|,
    psi     = min(tau,256).

The branch may therefore use the weight

    z^|J| tau^(-c_fresh) psi^(-c_old).                 (9)

This is a bound for the common rank potential, not a probability law.
It is valid for every tau>=1. To justify the product, process constraints
in the old-source span first. Their charged coefficients retain the
ordinary bound 1/psi. Then impose one pivot equation for each charged
fresh coordinate. These pivots are mutually independent and independent
of every old-source form. Each contributes exactly the factor 1/tau.
Further rejected branches or additional old equations can only decrease
the potential. In particular, no coefficient-root factor is also charged
for a fresh coordinate already accounted for by its pivot.

The outgoing-rank envelope can be tightened before summing over J.
If the incoming residual state rank is r and k=|I|, then

    rank((span(Q,X)+W)/W) = r+k.

Each of the c_fresh pivot equations lies in span(Q,X) and removes an
independent dimension there. Any other zero equations can only lower
the remaining dimension. Since Q' belongs to this span, every live branch
satisfies

    1 <= r' <= min(r+k-c_fresh, |H'|).                (10)

From a quiet state, Q'=X modulo W. Thus the fresh activation has residual
rank exactly k and its live-envelope entry is z^k v(k,I), without a
maximum over smaller ranks. Formula (10) is branch dependent; applying it
after aggregating different J values into one matrix entry would lose
the information needed for this improvement.

The exact route averages in Sections 2 and 3 are unchanged. Use fresh
active envelopes for the first b selected columns and ordinary active
envelopes thereafter. For a fixed target h, a column-selection recurrence
can track how many occupied columns have already been selected. The next
selected column is fresh precisely when that count is less than b.
An unselected column always uses C_0. This recurrence counts every h-set
once while retaining its chronological matrix product. Formula (7), with
one constant occupied-column matrix, describes only the earlier envelope.

The proof uses ordinary chronological zero equations. Lookahead equations
can involve source forms that have not yet been consumed and require a
separate argument. Nor does this proof apply to partially populated
columns: their initial shortening equations can couple the first occupied
columns to future columns. These local refinements do not by themselves
establish a multigroup certificate.

The prefix implementation in `check_two_groups.py` was audited read-only.
Its fresh-mask terms retain the disjoint charges in (9), and its live-rank
maximum applies (10) before summing output masks. Its column recurrence
advances the prefix only when a column is selected. No discrepancy was
found; the audit did not execute the program.

At the same marker tau=65536, z=.999, lambda=.9962, the root agent's
80- and 120-digit directed evaluations improve the h=33 bound from about
4.199e-4 to 2.566e-6. The h=64 bound drops from about 6.109e843 to order
10^766. The fresh-prefix bound still does not close the complete shared-
column two-group class, much less all two-group supports. These failed
upper bounds are not low-distance counterexamples.

The later [endpoint-tail refinement](ENDPOINT_TAIL.md) improves this same
class without changing the encoder. It closes the union over h=33,34,35
below 2^-117, including all adjacent pair and shared column choices.
Its broad h=33..64 bound still fails. Keep this earlier calculation as a
checkpoint for the specific loss repaired by the endpoint refinement.
