# Counting support histories and output cancellations together

2026-10-27. The field-quartet construction is unchanged. This note sharpens
the conditional field count. It proves a sufficient routing condition at
one difficult row count, but not the required ensemble routing theorem.
There is no whole-code 45% certificate yet.

The subsequent [whole-run analysis](QUARTET_BURST_CUT.md) proves an ensemble
bound for full-width bursts, uniformly over all row sets of size2100, and
extends it to histories with a bounded deterministic deviation charge.
General partially active transitions remain open.

## Why combine the counts?

The previous activity bound counts state supports. The turnover bound also
counts choices of zero outputs, but charges every history the same global
allowance. A nearly constant support and a rapidly changing support should
not receive the same charge. Moreover, the required number of cancellations
depends on that same support history.

Retain the setup and notation of [the coverage bound](QUARTET_HIGH_TURNOVER.md).
There are N=2^21 state positions, 4096 updates of width 512, and M=8192 outer
RS[256,128] rows over J=GF(2^32). Put q=2^32 and W=943718. For a nonzero
message, let n be its number of active outer rows. Let A be total state
activity, D the number of support deaths, C the covered inactive positions,
and w the terminal state weight. The initial state is zero; 0<=w<=512.

Fix the route, quartet partitions, an active-row set R, an exact state-support
history S, and a selected subset Z of zero outputs in consecutive-support
intersections. Write z=|Z|. The established field count is

    Pr_field[some message realizes (R,S,Z)]
        <= 2 q^(A-128n-C-z).                                  (1)

The probability uses the actual once-sampled accepted maps. The factor two
bounds their global acceptance density. Messages are not fixed independently
of the setup: (1) already counts candidate messages and state amplitudes.

For a low-output message, consecutive supports have total union size A+D.
Every position outside their intersection is forced nonzero in the output.
Thus some admissible Z has

    z=max(0,A+D-W),    0<=z<=A-D-w.                            (2)

The second inequality is necessary; if it fails, the history cannot have
output weight at most W.

## One normalized count for decorated histories

A decorated history records the support and which live-to-live transitions
carry selected zero outputs. Assign such histories a probability mass solely
as a counting device. It does not describe encoder states.

Given parameters beta,d,e>0 with beta<1 and d+e<1, use these transitions:

| Incoming support | Outgoing outcome | Probability |
| --- | --- | ---: |
| zero | zero | a=1-beta |
| zero | live | beta |
| live | zero | d |
| live | live, unselected | c=1-d-e |
| live | live, selected | e |

Begin each coordinate at zero. Conditional on these parameters, the 512
coordinate chains are independent and use the same transition probabilities.
This defines a normalized mass mu_(beta,d,e)(S,Z).

Now integrate the zero-state transition row against Dirichlet(1,1), and the
live-state transition row independently against Dirichlet(1,1,1). Call the
result mu. It remains normalized. Define

    n0=N-A+w,    n1=A-w,    I=A-D-w.

The transition counts are N-A-D, D+w, D, I-z, z, in table order. Consequently,

    mu(S,Z) = 1 / [(n0+1) binom(n0,D+w)
                    * binom(n1+2,2) multinomial(n1;D,I-z,z)]. (3)

For every allowed parameter choice, binomial and multinomial probabilities
are at most one. Hence

    mu(S,Z) >= mu_(beta,d,e)(S,Z)
                / [(n0+1) binom(n1+2,2)]
            >= mu_(beta,d,e)(S,Z) / 2^64,                    (4)

where (N+2)^3<2^64. These inequalities hold pointwise, simultaneously for
all parameters and decorated histories.

Use the normalized row-set prior

    nu(R)=1/[(M+1)binom(M,|R|)].

Sum (1) over triples satisfying

    2 q^(A-128n-C-z) <= 2^-74 nu(R) mu(S,Z).                 (5)

Normalization bounds the resulting probability by 2^-74, conditional on
any fixed route and quartet partitions. It includes every row count,
history, selected zero set, and free terminal state. Averaging geometry
preserves the bound. Outside this single exceptional event, no compatible
triple satisfies (5).

## A linear necessary condition for low output

Let L=-log2 mu_(beta,d,e)(S,Z). Its expansion is

    L = -N log2 a + A log2(a/c)
          + D log2(ac/(beta*d)) + w log2(c/beta) + z log2(c/e).

Combining (4) with the negation of (5) gives

    A-128n-C-z + [L+log2((M+1)binom(M,n))+139]/32 > 0.       (6)

Here 139=1+64+74. Assume c/e<=q. The coefficient of z in (6) is then
nonpositive. Replacing z by the lower bound A+D-W from (2) can only
increase its left side, even when A+D-W is negative. Every remaining
low-output message therefore satisfies

    u_A A - C + u_D D + u_w w + u_0 > 0,                   (7)

with

    u_A = log2(a/e)/32,
    u_D = -1 + log2(ac^2/(beta*d*e))/32,
    u_w = log2(c/beta)/32,
    u_0 = W-128n + [-N log2 a - W log2(c/e)
                    + log2((M+1)binom(M,n)) +139]/32.

The helper `selected_prior_weights` derives these combined coefficients
first, then rounds each upward using rational log intervals. This order
matters because W multiplies the eliminated coefficient with a negative
sign. Nonnegative A,D,w make the rounded score a valid upper bound.

No union cost is needed when choosing beta,d,e for a history class: (4)
puts all these consequences on the same good event. The current replay uses
beta=d=1/10 and e=1/6. These are proof parameters, not encoder parameters.

## A sufficient routing condition at r=2100

At r=M-n=2100, the certified activity good event gives A>=743225. Directed
evaluation of (7) gives approximately

    C + 0.7443 D < 111080.924 + 0.07603 A + 0.08983 w.

These decimals are explanatory only; replay uses rational upper coefficients.
Suppose the geometry additionally satisfies, for every support history and
every inactive-row set of size 2100,

    C + 0.74 D >= 0.24 A - B.                              (8)

Substituting (8) in the rounded score, using A>=743225 and w<=512, gives

    score <= u_0 + (u_A-0.24)*743225 + B + 512 max(0,u_w).

The discarded D coefficient is nonpositive. At B=5000 the upper score is

    -50486951550610057 / 8796093022208 < 0.

This contradicts (7). The same directed check permits every integer
B<=10739. Thus an ensemble bound for (8) with a sufficiently small offset
would finish this row-count class, without the full transfer calculation.
The geometry tail and the other row-count classes remain unproved.

The coefficient 0.74 is important: the older cut bound charged D at weight
one. Its excess cannot simply be dropped when comparing against (7).

## Fixed-geometry diagnostic and numerical scope

The integer min-cut helper now accepts a separate death-edge capacity.
For slope numerator 24, coverage denominator 100, and death capacity 74,
it minimizes 100C+74D-24A. It independently recounts the resulting support
and checks flow conservation. The short-run erasure threshold is now
24*length>=74 for every nonterminal run of a minimizer.

For the existing full-size diagnostic seed 20261027 and fixed inactive
set {0,...,2099}, it verifies (8) for all histories with exact B=3301.30.
The explicit chronological feasible flow proves the weaker B=9918.38 on
that same geometry. Both are below the sufficient threshold 10739.
Neither result is a probability bound or a certificate for all row sets.
No bad-code witness was found.

Two neighboring slopes use the same graph, fixed row set, and death weight.
The sufficient offset is recomputed from the same rounded joint score:

| Slope | Largest sufficient integer B | Exact min-cut B | Chronological-flow B |
| --- | ---: | ---: | ---: |
| 0.24 | 10739 | 3301.30 | 9918.38 |
| 0.25 | 18171 | 3988.44 | 14971.01 |
| 0.26 | 25604 | 7832.88 | 23175.07 |

All three sampled inequalities suffice for this row-count class. At slope
0.26 the optimizing support is almost full: A=2080885, C=533010, D=253.
Its geometry is different from the small-cut optimizer at slope0.24,
which has A=47026, C=636, D=9931. The slope0.26 run took about56 seconds;
the other two took about9 and10 seconds.

Inserting the joint score into the old FLOAT transfer does not close it.
At theta0.13, ell0.03, tau0.05, the refined 17-point p grid, rank marker1.1,
and joint marker1.5, the screen is +306806.85. Its artificial regional flow
alternates near weights2 and460. A preceding rank/turnover test gave
+296584.30. Both are worse than the best unmarked screen +271673.74.
The useful advance is the direct sufficient condition (8), not those screens.

## The ensemble count still overcounts geometry witnesses

A first-moment geometry test at slope0.26, death weight0.74, and b=256 gives
a FLOAT log2 union moment +1500177.33. With the complete row-set union already
included, its 60-bit offset screen is187529.67, far above the sufficient
25604. The sampled cuts therefore do not justify an ensemble claim.

Its dominant artificial regional bridges are approximately509->486->6 and
6->485->509. One spends several updates almost full and then almost empty;
the other grows in reverse. These are weighted counting trajectories, not
observed encoder behavior or a lower bound on setup failure. This diagnostic
motivates counting minimizing cuts, rather than every support history.

### Canonical minimizing cuts

For fixed geometry and an inactive-row set, maximize

    F(S)=lambda*A-C-tau*D.

Among maximizers, choose one with the fewest selected state positions. The
geometry has finitely many supports, so this choice exists. Any violation
of the proposed expansion inequality has such a maximizing witness.

Assume lambda+tau<=1. If an inactive input has exactly one selected state
neighbor, delete that neighbor. Let k>=1 be the number of inactive input
positions that become uncovered. This deletes one selected state and changes
D by at most one, including at the free terminal boundary. Therefore

    F(after)-F(before) >= -lambda+k-tau >= 0.

A strict increase contradicts maximality; equality contradicts the chosen
minimum cardinality. Thus no covered inactive input has incidence one.

A stronger local rule follows by deleting every incoming selected state
in one quartet. Let a>0 be their number, and h their overlap with the next
support. Let f be the number of inactive input positions outside that next
support. All f become uncovered. For each deleted coordinate, the change
in deaths is p+v-1, where p and v are its preceding and following support
bits. Since p<=1, the total change is at most h. Hence

    F(after)-F(before) >= -lambda*a+f-tau*h.

The same maximality and tie-breaking argument gives

    f < lambda*a+tau*h.                                   (9)

This rule does not require lambda+tau<=1. We use that assumption for the
additional inactive-birth exclusion at empty incoming quartets.

The moment helper's `--canonical` option imposes (9), and excludes inactive
births when a=0. Every geometry retains at least one maximizing cut, so the
filtered sum still bounds the probability of a geometry violation. These
rules restrict witnesses to a bad geometric inequality; they must NOT be
imposed on arbitrary message histories in the field-distance count.

The simple incidence-one filter alone changes the offset screen only from
187529.67 to187422.95. The stronger quartet-erasure rule gives187364.56,
still far above25604. These small improvements are evidence against relying
on the local rules alone, not evidence that the encoder has bad distance.

## Checks and next discriminating question

`selected_prior_certificate()` exactly checks three small decorated-history
normalizations, 268 pointwise comparisons, and 1680 selected-zero elimination
inequalities. Sixteen exact weighted paths verify the free-terminal factor
in the full diagnostic. `conditional_weighted_cut_closure()` verifies the sufficient
condition with directed rational coefficients. The min-cut helper checks
393216 exhaustive toy cuts, six feasible flows, and 61440 run erasures.
Two additional exhaustive toy cases check the canonical condition with
strict and tied objectives. The slope0.25 full-size replay also checks the
canonical rules on its exact minimizing support. The geometry-moment
helper checks50 rational local coefficients and9 canonical boundary cases.

All runs are serial and remote. The two FLOAT screens retain their 120s and
128MiB numeric-array limits. The graph run uses a 120s timeout and 768MiB
process limit; its temporary build is removed. No raw data is retained.

Next determine whether the weighted-cut offset admits a usable ensemble
tail, including the complete row-set union. Analyze whole runs or loss in
the chronological flow, rather than repeating the failed one-step tilts.
The sampled chronological-flow offset leaves little room for that union;
its exact min-cut leaves more.
Do not infer either tail from a sampled graph or silently replace the
once-sampled minimum-distance claim by a fixed-message estimate.
