# A forward bound for short excursions

Status: a reviewed analytic complete-class bound and an exact finite
certificate for the sparse-activity wedge recorded below, not a full-code
certificate. This argument uses the actual one-diagonal inner, fixed MDS
outer, and packet route. It retains outer support information discarded
by the coarse [active-union gate](ACTIVE_UNION_GATE.md). Root owns numerical
evaluation; the recorded checker was independently audited read-only.

The coarse gate can fail even when its route moment is exact. An activity
union can support an outer word without making its required resets likely.
For a two-batch excursion, a sparse outer activation forces a dense state
before reset. Its reset probability therefore gives a forward bound that
depends on the number of active outer rows, not only on whether their
shortened space is nonzero.

## Fix the route and the source before exposing inner coins

Use the notation of the active-union gate: N=bn=MB, k=B/2, p rows per
packet, G=M/p groups, and b dividing M. Fix a marked excursion layout
occupying v batches in c intervals, and let U contain all its coordinates.
Fix the packet route. The source space is

    C_U = C_out intersect F^U,       d=dim C_U.

For group g, let T_g(U) be the number of regions whose routed packet lies
in U. Because U contains whole packets, every row of this group has the
same occupancy. Define

    S=#{g:T_g(U)>k},       R=pS.

Only these R rows can have nonzero shortened codewords. The exact MDS
formula gives

    d=p sum_g(T_g(U)-k)_+ <= kR.                             (1)

Every outer row contributes one coordinate per region, and each batch
lies within one region. Thus every x in C_U satisfies

    wt(x restricted to any one batch) <= min(b,R).           (2)

This holds for every source word in the space, not only for a word fixed
before the route. We now fix one nonzero projective source representative
x after the route and before all independent inner coefficients.

## A two-batch reset

The coordinate-form recurrence is

    V_t=D_t Z_t,       Y_t=X_t+r_t .* V_t,
    Z_(t+1)=V_t+Y_t,

where each D_t is fixed and superregular and all coordinates of r_t are
independent uniform field elements. A two-batch excursion starts with
Z_a=0 and X_a!=0. Therefore

    Y_a=X_a,       Z_(a+1)=X_a,       V_(a+1)=D_(a+1)X_a.

Superregularity gives wt(x)+wt(D_t x)>=b+1 for every nonzero x. By (2),

    wt(V_(a+1)) >= b+1-min(b,R).                             (3)

The required reset is Z_(a+2)=0. It is equivalent to

    X_(a+1)=(1+r_(a+1)) .* V_(a+1).

Its probability is zero if X_(a+1) is nonzero outside supp V_(a+1).
Otherwise each nonzero coordinate of V_(a+1) fixes one independent
coefficient, so the probability is exactly q^(-wt(V_(a+1))). Hence it is
at most q^(-(b+1-min(b,R))).

For several prescribed two-batch resets, expose the setup chronologically.
If a required earlier activation or state condition has failed, the source
cannot realize the layout. Otherwise its current reset uses fresh
independent coefficients and obeys (3). Multiplying these conditional
bounds is valid even when other, longer excursions occur between them.
No future nonzero event is used to condition and renormalize the coins.

## Count source lines and short intervals

Let m_2 count the guaranteed closed excursions of length exactly two.
The first c-1 intervals have minimum length two and the last has minimum
length one. Their excess length is v-(2c-1). Every longer closed interval
uses at least one excess batch. Thus every layout in the class has

    m_2 >= max(0,3c-v-2) =: m.                              (4)

The last interval receives no reset charge; the terminal remains free.
There are (q^d-1)/(q-1) projective source lines in C_U. Fixing their
representatives before the inner setup and summing the preceding reset
probabilities bounds the conditional existence event by

    min(1, (q^d-1)/(q-1)
                 q^(-m(b+1-min(b,R)))).                    (5)

No output-weight or other support restrictions were used in (5), so it
bounds every output realizing the layout. Its setup quantifiers permit a
bad source to be chosen after all setup randomness is revealed.

Choose an integer row threshold 0<R0<=b. On R<=R0, (1) and (5) give

    Pr_inner[an eligible source exists | route]
      <= min(1, q^((k+m)R0-m(b+1))/(q-1)).                  (6)

If R0=0, this part is zero because C_U is the zero space. Small-R and
large-R routes form a partition; their bounds are added, not multiplied.

## Average the active-group tail under the actual route

For the fixed layout let upsilon_j=|U intersect region j|/M. A single
group's occupancy is a sum of independent Bernoulli variables with these
probabilities. Write

    mu_exact=Pr[sum_j Bernoulli(upsilon_j)>k].

For any r distinct source groups, the product of their bad-occupancy
indicators is increasing in every regional membership indicator. The
monotone without-replacement comparison from
[Partial-packet moments](SHUFFLED_PACKET_ROUTE.md), applied one region at
a time, gives

    Pr[these r groups all have occupancy>k] <= mu_exact^r.

The comparison does not assert independence of the actual groups.
A factorial-moment bound consequently gives

    Pr_route[S>=r] <= binom(G,r) mu_exact^r.                 (7)

Indeed binom(S,r)>=1 on S>=r, and its expectation is the sum over the
r-element group subsets. A useful placement-independent bound is

    mu_exact <= E binom(sum_j Bernoulli(upsilon_j),k+1)
              = e_(k+1)(upsilon_1,...,upsilon_B)
              <= binom(B,k+1)(v/n)^(k+1).

The last inequality is the elementary-symmetric-mean inequality, using
sum_j upsilon_j=Bv/n. Define

    mu=min(1,binom(B,k+1)(v/n)^(k+1)),
    r=floor(R0/p)+1,
    Tail(v;R0)=min(1,binom(G,r)mu^r).                        (8)

If r>G, set the tail to zero. Since R=pS, (8) bounds Pr[R>R0]. This is
an averaged bound for each fixed layout; the layout is not chosen after
the route without subsequently paying its multiplicity.

## Complete class bound and scope

With the exact layout count

    J(v,c)=binom(v-c,c-1)binom(n-v+c,c),

the existence probability for the entire (v,c) class is at most

    J(v,c) min(1,
        q^((k+m)R0-m(b+1))/(q-1) + Tail(v;R0)),              (9)

where m is (4) and 0<R0<=b. Replace the first term by zero at R0=0.
The feasible conditions and deterministic touched-region exclusions in
the active-union gate still apply. One may take the minimum of (9) and
any independently justified bound for the same complete class.

For the planned K20 parameters b=128, B=64, p=4, the class (v,c)=(128,64)
has m>=62. Choosing R0=80 makes the small-R exponent

    (32+62)*80-62*129 = -478.

The remaining routes have at least 21 active groups and are covered by
their factorial tail. The following check includes both terms, every
layout count, and every (v,c) pair in its claimed range.

## Exact finite sparse-activity wedge

For K=2^20, N=2^21, q=2^32, B=64, k=32, p=4, b=128, the integer checker
[check_short_excursion_forward.py](theory/check_short_excursion_forward.py)
covers precisely the feasible classes

    1 <= v <= 512,       1 <= c,       2c-1 <= v,
    3c-v-2 >= ceil(v/4).                                  (10)

Thus every counted layout has at least ceil(v/4) guaranteed closed
two-batch excursions. The guarantee counts only the first c-1 intervals;
no final-state constraint is added. The check covers every placement,
support mask, and source amplitude in these classes. It does not impose
an output-weight cutoff and therefore also covers all their outputs of
weight at most floor(0.35N).

There are 10,954 checked classes and 54 additional classes excluded by the
deterministic touched-region test. For every checked class the program
chooses one row threshold R0, using only integer arithmetic. It verifies
separately

    J(v,c) q^((k+m)R0-m(b+1)) < 2^-80,
    J(v,c) binom(G,r)
       [binom(B,k+1)(v/n)^(k+1)]^r < 2^-80,                (11)

with the first term replaced by zero when R0=0. The first inequality
deliberately drops the favorable projective divisor q-1. The second is
verified by exact integer cross multiplication, not logarithms. Hence
the probability that any nonzero source realizes any class in (10) is

    < 2*10954 / 2^80 = 21908 / 2^80 < 2^-65.               (12)

The weakest recorded large-row bound still exceeds 570 bits, at
(v,c,R0,r,m)=(512,216,100,26,134). This diagnostic is not used in (12).
In particular, (10) includes the complete (66,33) class that defeats the
coarse active-union proof; the repair comes from retaining the sparse
outer activation constraint, not from improving its route moment alone.

Reproduction command from this workstream directory:

```sh
python3 -B theory/check_short_excursion_forward.py
```

The root's recorded serial run used a 60-second timeout and a 256 MiB
memory limit, passed in 0.373 seconds, and wrote no artifacts. Independent
source review checked the class coverage, zero-source exclusions, strict
integer comparisons, and union accounting. The checker itself also has
a 60-second elapsed-time guard; the memory limit is an external runner
constraint.

This gate becomes weak when few closed excursions have length two. It
does not yet bound longer excursions using their intermediate state
supports. The next extension should preserve the activation constraint
while deriving the actual one-diagonal transition probabilities; the
two-diagonal uniform-reset bound is unavailable here.
