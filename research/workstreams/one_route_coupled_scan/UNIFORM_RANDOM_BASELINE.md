# A uniform-random baseline with one permutation

**Status: reviewed asymptotic proof draft.** The argument below establishes
every fixed distance strictly below the q-ary Gilbert--Varshamov bound. It gives no
concrete finite-length certificate or implementation timing. It is a new
baseline, not a certificate for the coupled scan or its small coefficients.

The ordinary expected-word calculation transfers to larger fields, but it
can be dominated by rare routes that make many messages bad together. The
proposed proof first excludes those routes through properties of shortened
outer codes. It then counts outputs using the inverse recursive encoder.

## 1. Ensemble and target

Fix a prime power q, a rational rate R in (0,1), and a distance delta with

    0 < delta < 1-1/q,
    H_q(delta) < 1-R,

where

    H_q(x) = [x ln(q-1) - x ln x - (1-x) ln(1-x)] / ln q.

The alphabet q is fixed as N grows. Binary extension fields q=2^m are
included, but this statement does not allow m to grow with N.

Take an admissible sequence N=MB with k=RB integral and B=Theta(ln N).
Independently sample uniform injective linear maps

    C_i : F_q^k -> F_q^B,       1<=i<=M,

and one uniform permutation Pi of the N coordinates. Let C be the image
of the concatenated block code after Pi. Thus C has dimension K=Mk=RN.
Its block images are independent uniform k-dimensional subspaces; their
random bases do not affect distance.

Choose memory L=Theta(ln N). Independently sample every a[i,j] uniformly
from F_q, for 1<=i<=N and 1<=j<=L. On routed input x define

    y_i = x_i + sum_(j=1)^L a[i,j] y_(i-j),
    y_i=0 for i<=0.                                         (1)

All setup is fixed once for every message. Equation (1) defines an
invertible linear map: its inverse recovers x_i from the displayed y's.
The code therefore has rate R and uses one global permutation.

### Asymptotic theorem

For every fixed q, R, and delta above, sufficiently large constants in
B=Theta(ln N) and L=Theta(ln N) give

    Pr_setup[d_min <= floor(delta N)] = o(1).                (2)

One sufficient memory schedule in the argument below is

    L >= k + ceil(8 log_q(N+1)) + 2.                        (3)

The constant in B depends on q, R, and delta. Increase that constant
first, then choose L by (3). No optimization of these constants is claimed.
For integral lengths, one may choose B among multiples of the denominator
of R and take the least such B above a fixed multiple of ln(MB).

Direct encoding uses O(NB+NL)=O(N log N) field operations. The outer has
M separate B-by-k generators, not one dense N-by-K generator. The inner
uses L field multiplications per output. This is a proof baseline; uniform
large-field multiplication is not asserted to meet a millisecond target.

## 2. Exact expected enumerators and their limitation

Let X mark input weight and S output weight. For the inner, retain the
number z of trailing output zeros, clipped at L. Put Q=q-1. The nonzero
entries of its transfer matrix are

    T[z,0]   = (1+QX) Q S/q,        z<L,
    T[z,z+1] = (1+QX)/q,            z<L,
    T[L,0]   = Q X S,
    T[L,L]   = 1.                                             (4)

When z<L, a fresh uniform coefficient vector makes the next output
uniform for every fixed input symbol. When z=L, the next output equals
the input. Consequently the exact expected input/output enumerator is

    F_N(X,S) = e_L^T T(X,S)^N 1
             = E_inner sum_x X^wt(x) S^wt(I(x)).              (5)

For the outer, put rho=(q^k-1)/(q^B-1). Its expected weight enumerator is

    A(X) = [1 + rho ((1+QX)^B-1)]^M.                         (6)

Write A_h=[X^h]A(X). A uniform permutation makes each fixed weight-h
support uniform. Although it does not randomize the nonzero values,
the inner's output-support law depends only on the input support: every
live step is uniform and every quiet step preserves zero/nonzero status.
Thus, with e=floor(delta N), the expected number of bad nonzero words is

    sum_(h=1)^N A_h / (binom(N,h) Q^h)
                   * sum_(d=0)^e [X^h S^d] F_N(X,S).        (7)

Equation (7) gives a valid whole-code first moment. Dividing it by q-1
also bounds failure because a nonzero codeword brings its scalar multiples.
However, that does not make (7) useful near q-ary GV.

For example, if all B coordinates of one block land in the last e positions,
all its q^k-1 nonzero words have output weight at most e, regardless of the
inner coefficients. The event has probability binom(e,B)/binom(N,B),
asymptotic to delta^B for B=Theta(log N). Its contribution to the expected
word count is roughly M(q^R delta)^B. This can grow while the probability
of any such route tends to zero. The following argument conditions on
good routed outers rather than paying that amplitude multiplicity.

## 3. Reverse the inner for a fixed output

For a coordinate set S, define the shortened dimension

    d_C(S) = dim_Fq {x in C : supp(x) subset S}.

Fix the routed outer C and a nonzero output vector y, before sampling
the inner. Define

    A(y) = {i : at least one of y_(i-1),...,y_(i-L) is nonzero},
    U(y) = A(y) union supp(y).

The inverse of (1) is independent uniform on coordinates in A(y), and
equals y on the remaining coordinates. Independence follows from the
independent coefficient rows a[i,*]. Its intersection with C is either
empty or an affine fiber of size q^d_C(A(y)). Hence

    Pr_inner[I^(-1)(y) in C] <= q^(d_C(A(y))-|A(y)|).         (8)

This uses d_C(A), not d_C(U); monotonicity later permits the upper bound
d_C(A)<=d_C(U). Every possible inverse input is supported in U. In
particular, if d_C(U)=0, no nonzero y can qualify. The inverse is always
injective, and its first active input equals the nonzero first output.

The proof sums (8) over all eligible output vectors. It does not fix a
message before setup. For each completed setup, every low-weight nonzero
message gives one of these outputs, so this sum bounds the expected
number of bad codewords conditional on C.

## 4. Activity intervals and the free terminal

List the nonzero output positions in increasing order. Start a new cluster
whenever the gap from the preceding nonzero position exceeds L. For a
cluster with first position s and last position t, use the interval

    [s, min(N,t+L)].                                       (9)

These c intervals are disjoint but may be adjacent. A gap of L+1 between
nonzero outputs produces adjacent intervals, not overlapping intervals.
Their union is U(y). Their c first positions are precisely the nonzero
outputs with an all-zero preceding memory. Therefore

    |A(y)| = |U(y)|-c.                                    (10)

Every interval except possibly the last ends with L forced zero output
symbols and has length at least L+1. The last interval is truncated at N
when necessary. We impose no final-state constraint, and may ignore even
the last interval's genuine trailing zeros. In particular,

    c <= floor(N/(L+1))+1.                                (11)

For a fixed interval representation of total size u, the number of output
assignments is at most

    q^(u-L(c-1)).                                         (12)

Requiring nonzero cluster endpoints and restricting internal zero runs
would only decrease this count. For a weight threshold e, a second bound
is the q-ary Hamming-ball volume inside those u positions.

There are at most (N+1)^(2c) representations with c intervals. For all
c<=r=floor(N/(L+1))+1, the sharper total bound

    sum_(c=1)^r binom(N+2c,2c)
        = exp(O(N log L/L)) = exp(o(N))                    (13)

is useful. Weakly ordered endpoint cuts give the binomial upper bound
and allow adjacent intervals. We count representations, not just connected
components of U.

## 5. A good routed outer on small interval unions

Put p=1-1/q. Strict concavity of H_q and H_q(p)=1 imply
H_q(x)>=x/p for 0<=x<=p. The entropy hypothesis therefore permits constants

    delta/p < a_0 < a_1 < 1-R.                            (14)

We first show that, with probability 1-O(N^-2), every set U of size at
most a_0 N represented by c intervals satisfies

    d_C(U) <= k(c-1).                                    (15)

Fix such a target U, of size u. Before routing, its coordinates form a
uniform u-subset of the N source coordinates. Generate that subset instead
by independent Bernoulli indicators of probability alpha=u/N, conditioned
on their total being u. That conditioning costs at most N+1: u is a mode
of Bin(N,u/N), so its mass is at least 1/(N+1). The case u=0 is immediate.

Call a block bad if its local code contains a nonzero word supported on
its selected coordinates. If a block selects t coordinates, projective
line counting gives

    Pr[block bad | t]
       <= (q^t-1)(q^k-1) / ((q-1)(q^B-1)).                (16)

Under the unconditioned Bernoulli experiment, t is binomial with mean
at most a_0 B. Split at t=a_1 B. Hoeffding's bound and (16) give

    Pr[block bad] <= p_B,
    p_B := exp(-2(a_1-a_0)^2 B)
                       + 2 q^(-(1-R-a_1)B).              (17)

All blocks are independent in this unconditioned experiment, including
their random subspaces. For any specified c blocks, the probability they
are all bad after conditioning is therefore at most (N+1)p_B^c.

Choose the constant in B=Theta(ln N) so that, for all sufficiently large N,

    p_B <= (N+1)^(-6).                                   (18)

For example, choose B>=C ln N with C strictly larger than
7/min{2(a_1-a_0)^2, (1-R-a_1)ln q}. Union over the interval representations
and choices of c bad blocks. The probability that (15) fails anywhere is
at most

    (N+1) sum_(c=1)^M [(N+1)^2 M p_B]^c = O(N^-2).        (19)

Outside this event, at most c-1 blocks contribute any shortened dimension,
and each contributes at most k. This proves (15). For c>M the conclusion
is automatic. Increasing the constant in B strengthens the polynomial
failure bound; (19) is not a concrete finite-length security claim.

## 6. A second good event for larger interval unions

We also need, uniformly on the interval representations in (11),

    d_C(U) <= [|U|-(1-R)N]_+ + epsilon_N N,
    epsilon_N -> 0.                                      (20)

Here [x]_+=max(x,0). Two separate events establish this statement.

First consider a local B-by-k random generator. For any restriction to
n rows, the probability that its rank falls d below min(n,k) is at most

    C_q q^(-d(|n-k|+d)),                                 (21)

for a constant C_q depending only on q. One way to obtain (21) is to count
rank-r matrices: their number is at most a constant times
q^(r(n+k-r)), then sum over r<=min(n,k)-d. Conditioning the full B-by-k
matrix to be injective costs only another constant. Its full-rank
probability is bounded below by the positive product
product_(j>=1)(1-q^-j).

Let d_B=ceil(B^(2/3)). Union (21) over every coordinate restriction and
every block. With failure at most

    M 2^B C_q q^(-d_B^2) = o(1),                         (22)

each restriction loses at most d_B ranks. This event depends only on the
outer generators and holds for every later route.

Let t_i be the number of coordinates of U routed from block i. On this
outer event,

    d_C(U) <= F_U + M d_B,
    F_U := sum_i [k+t_i-B]_+.                            (23)

For each fixed target U, use the Bernoulli experiment from Section 5.
The positive-part function is 1-Lipschitz, and
E|t_i-alpha B|<=sqrt(B)/2. Therefore

    E F_U <= [u-(1-R)N]_+ + N/(2 sqrt(B)).                (24)

More explicitly, [k+t_i-B]_+ is at most
[k+alpha B-B]_+ + |t_i-alpha B|; sum this inequality over the M blocks.

Changing one Bernoulli indicator changes F_U by at most one. Bounded
differences therefore yields

    Pr[F_U > E F_U + eta N] <= exp(-2 eta^2 N).            (25)

Conditioning the total adds the same factor N+1. Union (25) over (13).
Taking eta_N=(log L/L)^(1/4) makes its negative exponent dominate the
O(N log L/L) count. Thus, with probability 1-o(1), (20) holds with

    epsilon_N = eta_N + 1/(2 sqrt(B)) + d_B/B = o(1).      (26)

No independence between the two good events is required. Fix any routed
outer satisfying both (15) and (20) for the remaining inner calculation.

## 7. Small activity unions

Consider nonzero outputs whose interval union U has size at most a_0 N.
If c=1, (15) gives d_C(U)=0, so no such output can belong to the code.

For c>=2, multiply the assignment count (12) by (8), using (10) and (15).
For each interval representation the conditional expected bad-word count
is at most

    q^(u-L(c-1)) q^(k(c-1)-u+c)
        = q^((k-L)(c-1)+c).                              (27)

The memory schedule (3) bounds the sum of (27) over all representations by

    sum_(c>=2) (N+1)^(2c) q^((k-L)(c-1)+c)
       <= sum_(c>=2) (N+1)^(8-6c) q^(2-c)
       = O(N^-4).                                       (28)

This count even permits arbitrary output weight. The L forced zeros at
each completed activity interval pay for its source dimension and endpoint
choices. The final interval remains free throughout the calculation.

## 8. Large activity unions and the GV exponent

Now take u=aN>=a_0 N and output weight at most e=floor(delta N). By (14),
delta/a<p. The q-ary ball-volume bound gives at most

    q^(u H_q(delta/a))                                   (29)

assignments inside U. For completeness, evaluate
(1+(q-1)S)^u at S=(delta/a)/((q-1)(1-delta/a))<1 and multiply by
S^(-delta N); its nonnegative coefficients give (29).

Equations (8), (10), and (20) make the contribution of each representation
at most

    q^(N G(a) + epsilon_N N + c),
    G(a) := a H_q(delta/a) - min(a,1-R).                  (30)

This exponent is uniformly negative for a in [a_0,1]. On
[a_0,1-R], H_q(delta/a)<1, so G(a)<0. On [1-R,1],
the function a H_q(delta/a) is increasing: its derivative is
-log_q(1-delta/a)>0. Consequently

    G(a) <= H_q(delta)-(1-R) < 0

there. Continuity on the compact interval gives a fixed negative upper
bound. Since epsilon_N=o(1), c=o(N), and the total number of representations
is exp(o(N)), their complete contribution is q^(-Omega(N)).

## 9. Conclusion and scope

For every good routed outer, the conditional first moment over the actual
uniform inner is O(N^-4)+q^(-Omega(N)). Markov's inequality bounds the
existence of any bad nonzero codeword. Add the outer-good-event failures
from (19), (22), and (25). Their sum tends to zero, yielding
statement (2).

The proof separates rare bad outer/routes from the number of messages they
make bad. It does not assume that an adversary chooses its message before
setup. It uses one permutation, independent local random outer maps, the
actual uniform recursive coefficients, and a free terminal state.

Root and two independent analytic reviewers checked the inverse fibers,
interval endpoints, fixed-size conditioning, rank-deficit union, concentration,
and entropy counts. They found no gap. A bounded exhaustive check additionally
verified the interval identities and assignment bound over 262140 nonzero
length-eight words over four symbols, for L=1,2,3,4, covering 322 interval
representations. This finite check is not a distance certificate.

A finite certificate would additionally require explicit choices of B, L,
all error thresholds, and numerical upper bounds. The asymptotic constants
can be impractical, especially near GV. No finite certificate, optimized
implementation, or transfer to small coefficients is asserted.

The next step is to state explicit finite versions of the two outer properties
and the inverse-inner bound. Then replace one component at a time, beginning
with a block-dispersing structured route. Keeping these interfaces explicit
separates a construction defect from slack introduced by the analysis.
