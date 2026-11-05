# Fresh Toeplitz batches for the uniform baseline

**Status: analytic extension; finite instantiations are certified
[separately](STRATIFIED_CLOSED_CERTIFICATE.md). No measured implementation.**
This note replaces the scalar recursive inner in
[the uniform baseline](UNIFORM_RANDOM_BASELINE.md) by a batched inner. It
preserves the exact inverse distribution needed by that proof. A batch can
occupy one quarter of the memory window without losing the asymptotic
sub-GV result. No transfer to the earlier small-coefficient coupled scan is
claimed.

The finite statement below conditions on an arbitrary routed outer code.
It therefore also specifies what a future cheaper outer and structured
route must provide. The inner's coefficients remain uniform over the full
field; replacing them by small coefficients requires a separate argument.

## 1. Construction and its inverse

Fix a finite field F_q and integers b,m,n>=1. Put

    N=nb,                 L=mb.

Fix a linear routed outer code C subset F_q^N, independently of the inner
setup. Partition its input x into n consecutive batches X_t in F_q^b.
The state S_t contains the preceding m output batches, in a fixed order,
and has L scalar coordinates. Initially S_1=0.

For each t, independently sample b+L-1 coefficients d_(t,r) from F_q,
with -(L-1)<=r<=b-1. Define a b-by-L Toeplitz matrix by

    T_t[a,j] = d_(t,a-j),       0<=a<b, 0<=j<L.

The encoder computes

    Y_t = X_t + T_t S_t,                                  (1)

then removes the oldest batch from the state and appends Y_t. There is
no feedback within a batch. Thus (1) is a different recurrence from
performing b consecutive scalar updates. Every T_t is sampled once and
fixed for all messages. The terminal state is unconstrained.

For any output y, its preceding-output states S_t(y) are determined. The
inverse is

    X_t = Y_t - T_t S_t(y).                               (2)

Consequently the inner is an invertible linear map, for every setup. It
preserves the outer's dimension and rate and adds no global permutation.

### Exact joint uniformity within an active batch

Fix a nonzero S in F_q^L. Let ell be its first nonzero coordinate. Then

    (TS)_a = d_(a-ell) S_ell
             + sum_(j>ell) d_(a-j) S_j.                  (3)

The b coefficients d_(a-ell), for a=0,...,b-1, are distinct. After fixing
all other diagonals, process a in increasing order. The second term only
uses diagonals with smaller indices than a-ell. Equation (3) therefore
defines a triangular map with nonzero diagonal S_ell from these b
coefficients to TS. This map is bijective. Hence

    TS is uniform on F_q^b.                              (4)

This is joint uniformity of the entire batch, not only uniform marginals.
Fresh independent matrices at different batches give independence between
their inverse-input batches for a fixed y. Sharing diagonals between
batches would remove this justification.

## 2. Activity intervals and assignment counts

Call an output batch nonzero when Y_t is not the zero vector. For a fixed
nonzero y, let A_batch contain the batches whose preceding state is nonzero.
Let U_batch be their union with the nonzero output batches. Expand both
sets to scalar-coordinate sets A(y) and U(y).

List the nonzero output batches. Start a new cluster when consecutive
listed indices differ by more than m. If a cluster starts at s and its
last nonzero batch is t, represent it by

    [s,min(n,t+m)].                                       (5)

These c batch intervals are ordered and disjoint, but can abut. For
example, a gap of m+1 between nonzero batches gives adjacent intervals.
Their union is U_batch. The first batch of each interval is nonzero and
has preceding state zero. It is the only batch in that interval outside
A_batch.

Write v for the total number of batches in the intervals and u=bv. Then

    |U(y)|=u,                  |A(y)|=u-bc.               (6)

The scalar set A is itself a union of at most c intervals, after removing
the c activation batches. The first c-1 intervals have length at least
m+1 batches and end in m zero output batches. Thus they force L(c-1)
zero scalar outputs. The last interval remains free, including when it
is clipped at N. In particular,

    c <= floor((n-1)/(m+1))+1.                           (7)

For a fixed interval representation, put

    f=bv-L(c-1),                 P(z)=1+(q-1)z.

There are f unforced scalar output positions. Each of the c activation
batches must be nonzero, and f>=bc. For an integer threshold 0<=D<=N,
define

    W_(q,b)(f,c,D)
      = sum_(d=0)^D [z^d] (P(z)^b-1)^c P(z)^(f-bc).      (8)

All coefficients in (8) are nonnegative. The expression bounds the number
of actual outputs of weight at most D with that representation. It ignores
constraints on internal zero runs and on the last nonzero batch of each
cluster, so it can overcount. It simultaneously retains the forced zero
tails and the weight threshold. Useful weaker bounds are

    W_(q,b)(f,c,D) <= q^f,
    W_(q,b)(f,c,D) <= V_q(f,D),                           (9)

where V_q(f,D)=sum_(d=0)^min(D,f) binom(f,d)(q-1)^d. Also
W_(q,b)(f,c,D)=0 when c>D.

## 3. A finite whole-code bound for any fixed routed outer

For a scalar-coordinate set S, define its shortened outer dimension by

    d_C(S)=dim_Fq {x in C : supp(x) subset S}.

Fix y before sampling the inner. Equations (2)--(4) make its inverse input
uniform on an affine translate of F_q^A. The inverse equals y outside A.
Its intersection with C is empty or has q^d_C(A) elements. Therefore

    Pr_inner[I^(-1)(y) in C] <= q^(d_C(A)-|A|).           (10)

Every possible inverse is supported in U. Its first activation batch
equals the nonzero first output batch, so the inverse is nonzero. If
d_C(U)=0, the probability in (10) is zero.

Let J_(n,m) be the following convenient superset of the actual interval
representations: ordered disjoint nonempty batch intervals inside [1,n],
with the first c-1 lengths at least m+1 and the last length at least one.
Intervals may abut. For J in this family, define A_J by deleting the first
batch of each interval from its scalar union U_J. Write c_J and v_J for
its interval count and total batch length, and f_J=bv_J-L(c_J-1).

**Finite conditional bound.** For every fixed C and every 0<=D<=N,

    Pr_inner[d_min(I(C))<=D]
      <= min(1, (1/(q-1)) sum_(J in J_(n,m))
            1[d_C(U_J)>0]
            q^(d_C(A_J)-b(v_J-c_J))
            W_(q,b)(f_J,c_J,D)).                         (11)

For C={0}, take d_min=+infinity. To prove (11), sum (10) over all nonzero
outputs of weight at most D, classified by their actual representation.
Equation (8) bounds each class. A bad nonzero codeword brings q-1 distinct
nonzero scalar multiples, all of the same weight. Dividing the expected
count by q-1 bounds the probability of any bad word.

This proof permits the bad message to depend on the completed setup. It
does not choose a witness before sampling the inner. If an outer/route
event supplies the required dimension bounds except with probability
epsilon_outer, add epsilon_outer to the conditional bound once.

### Grouping the finite bound by interval count and size

For fixed c and v, the number of representations in J_(n,m) is

    K_(n,m)(v,c)
      = binom(v-m(c-1)-1,c-1) binom(n-v+c,c),             (12)

when m(c-1)+c<=v<=n, and is zero otherwise. Subtracting the minimum lengths
leaves a weak composition among c intervals. Their c+1 gaps are nonnegative
and sum to n-v; zero gaps allow abutment. These two compositions give (12).
The family includes impossible final clusters, which only increases (11).

Suppose integer caps D_U(v,c),D_A(v,c) hold for every such representation:

    d_C(U_J)<=D_U(v,c),        d_C(A_J)<=D_A(v,c).

Then (11) is at most

    (1/(q-1)) sum_(c,v) K_(n,m)(v,c) 1[D_U(v,c)>0]
      * q^(min(D_A(v,c),b(v-c))-b(v-c))
      * W_(q,b)(bv-L(c-1),c,D),                           (13)

clipped at one. If only D_U is available, monotonicity permits D_A=D_U.
Equation (13) is an evaluable interface, not an evaluated finite certificate.

## 4. Small activity unions: the batch-size constraint

Fix an integer k>=1 and a_0 in (0,1). Assume that every union U of c
intervals with |U|<=a_0 N satisfies

    d_C(U)<=k(c-1),                                      (14)

as in the uniform baseline. An output with one such interval is impossible.
For c>=2, use d_C(A)<=d_C(U) and the first bound in (9). The expected
bad-word contribution of one representation is at most

    q^(bv-L(c-1)) q^(k(c-1)-bv+bc)
      = q^((k-L)(c-1)+bc).                              (15)

There are at most (n+1)^(2c) representations. Put g=L-k and M_0=n+1.
If M_0^2 q^(b-g)<1, summing the resulting geometric series gives

    small-union failure
      <= M_0^4 q^(2b-g)
          / ((q-1)(1-M_0^2 q^(b-g))).                    (16)

One sufficient schedule is

    L-k >= 2b+ceil(8 log_q(n+1))+2.                     (17)

Its right side makes (16) O((n+1)^-4). The c=2 term already contains
q^(k-L+2b), so this argument does not permit b=L/2 when k>0. Increasing
L while preserving that ratio does not repair the term.

A simple admissible choice is

    L=4b,
    2b >= k+ceil(8 log_q(N+1))+2.                        (18)

More generally, any admissible fixed ratio b/L<1/2 permits (17) with sufficiently
large logarithmic parameters. Equation (18) is a clean proved ratio,
not a claim that one quarter is optimal.

## 5. Large activity unions: retain both the random coordinates and zero tails

Fix a rate R in (0,1). The dense outer condition is needed on A, not only
on U. Assume it holds uniformly for the interval sets under consideration:

    d_C(A) <= [|A|-(1-R)N]_+ + epsilon_N N,
    epsilon_N -> 0.                                    (19)

The uniform baseline establishes this property for all unions with the
number of intervals in (7). Removing each activation batch leaves another
union in the same class.

Fix q,R,delta with 0<delta<1-1/q and H_q(delta)<1-R. Put p=1-1/q, and
choose a_0>delta/p as in the baseline. Define the normalized ball exponent

    V_delta(r) = r,                         0<=r<=delta/p,
    V_delta(r) = r H_q(delta/r),             delta/p<r<=1.

The value at r=0 is zero. For f=rN and D=floor(delta N), (9) is at most
q^(N V_delta(r)). The first branch covers thresholds above the typical
weight of an unrestricted word; it avoids incorrectly extending the
entropy expression beyond its ball-volume range.

For a representation with u=aN>=a_0 N, put

    alpha=b/L,             t=L(c-1)/N.

Then f/N=a-t, while |A|/N=a-alpha t-b/N. Combining (9), (10), and (19),
the normalized exponent is at most

    V_delta(a-t)-min(a-alpha t,1-R)+epsilon_N+b/N.        (20)

For every alpha bounded above by a fixed alpha_max<1, the first two terms
are uniformly negative over a>=a_0 and 0<=t<=a. To see this, first suppose
a-alpha t>=1-R. Monotonicity of V_delta gives

    V_delta(a-t)-(1-R) <= H_q(delta)-(1-R)<0.

Otherwise put r=a-t. The expression becomes

    V_delta(r)-r-(1-alpha)t.                             (21)

Here V_delta(r)<=r. If t>0, the final term makes (21) negative. At t=0,
r=a>=a_0>delta/p, so V_delta(a)<a. Continuity on the compact parameter
set gives a uniform negative gap, including transitions between the cases.

Thus a constant-ratio batch does not lose the asymptotic GV exponent.
Counting only output weight and then bounding d_C(A) by d_C(U) would lose
the linear activation penalty bc. Equation (20) avoids that loss by keeping
the exact random-coordinate set A and the L(c-1) forced zeros together.

For b,L=Theta(log N), the total number of representations is exp(o(N)).
For example, weakly ordered endpoints give the bound

    sum_(c<=floor((n-1)/(m+1))+1) binom(n+2c,2c).

Its logarithm is O(N log L/L) when 1<=b<=L. In the fixed-ratio case,
n=Theta(N/L), so the simpler bound 2^O(n)=exp(o(N)) suffices.
Consequently (20), epsilon_N=o(1), and b=o(N) give total large-union
failure q^(-Omega(N)).

### Asymptotic corollary

Use the uniform baseline's independent local [B,k] random outers and one
uniform permutation, with B=Theta(log N), k=RB. Choose admissible lengths
with b=Theta(log N), L=4b, and (18). The baseline's two good outer events,
followed by (16) and (20), give

    Pr_setup[d_min<=floor(delta N)] = o(1)

for every fixed strict sub-GV distance. The alphabet q is fixed as N grows.
The small-union calculation, rather than the large-union exponent, imposes
the present restriction b/L<1/2. This remains an asymptotic statement;
no concrete K20 margin follows from it.

## 6. Cost: less setup immediately, less arithmetic only with a suitable kernel

A fresh dense b-by-L matrix uses bL random field coefficients per batch.
The Toeplitz matrix uses b+L-1. Over n=N/b batches, the stored coefficient
counts are

    dense:       NL,
    Toeplitz:    N(b+L-1)/b.

The working output window contains L field symbols. Storing only diagonals
preserves the setup reduction; materializing each dense matrix would not.
With b=L/4, the number of sampled field coefficients is O(N), rather than
O(NL). Fresh diagonals for every batch are part of the proof.

A direct Toeplitz matrix-vector product still uses bL field multiplications
and bL field additions per batch, including addition of X_t. The structure
alone does not make that implementation faster.

The product can instead be computed by polynomial multiplication. Define

    D(Z)=sum_(r=-(L-1))^(b-1) d_r Z^(r+L-1),
    S(Z)=sum_(j=0)^(L-1) S_j Z^j.

Then (TS)_a is the coefficient of Z^(L-1+a) in D(Z)S(Z). Let M_F(s)
denote the cost, in actual F_q operations, of multiplying polynomials of
length O(s). Computing the whole product gives a valid upper bound

    inner cost: O((N/b) M_F(b+2L)).                      (22)

At b=L/4 this is O((N/L) M_F(L)). If a suitable implementation gives
M_F(L)=O(L log L), the inner costs O(N log L), or O(N log log N) with
logarithmic memory. Even Karatsuba's O(L^(log_2 3)) field-operation bound
is asymptotically smaller than direct O(NL) encoding at this batch ratio.
Neither statement supplies a practical crossover or a measured runtime.

These are full-field operations. For F_(2^128), a multiplication is not
a GF256 multiplication or an XOR. Padding, scratch space, transform setup,
and the actual polynomial kernel must be included in any implementation
comparison. The dense random outer still costs O(NB) field operations;
this note only changes the inner.

The linked finite checkpoint combines the exact inverse-coordinate set A
with stratified-route moments. It certifies selected memory/batch sizes,
including one finite b=L/2 instance not implied by the asymptotic corollary.
The next gate is to benchmark direct and polynomial Toeplitz products at
those sizes. Preserve the scalar baseline as a control. No benchmark,
production edit, or paper change accompanies this note.
