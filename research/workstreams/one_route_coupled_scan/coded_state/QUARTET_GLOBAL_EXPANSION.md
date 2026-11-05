# Global quartet expansion and the low-rank remainder

The construction and the 45%/40-bit target remain fixed. This checkpoint
proves another all-history class theorem. It does **not** close minimum
distance for the whole code. The earlier rank, high-coverage, and
[full/empty burst](QUARTET_GLOBAL_BURSTS.md) theorems remain available;
local seven-type refinements remain parked.

## What the expansion bound measures

Work over J=GF(2^32), with N=2097152, 4096 updates of width512,
8192 RS[256,128] rows, and W=943718. Use the actual once-sampled accepted
maps, uniform independent quartet partitions, and regional permutations
specified in [the global reset](QUARTET_GLOBAL_RANK_RESET.md).
The terminal state is free.

Fix a sequence S of state supports, starting empty. Let A be its total
activity, D its total deaths, and w its terminal weight. Write

    m = sum_t min(|S_t|,512-|S_t|),
    I = A-D-w,                 z = max(0,A+D-W).

Here I counts positions where the state survives an update. Discard z>I:
such a history cannot have at most W nonzero outputs. Given the quartet
partitions, let Gamma_t be the union of quartets meeting S_(t-1), and B0
the births outside those live quartets. The covered input footprint is
V=union_t(Gamma_t union B0_t). Put

    F=|V|,    E=F-A,
    R=sum_i min(128,|V intersect row i|),    d=F-R.

Thus d is the dimension left by shortening the outer code to V. It depends
on the route as well as the support history and partitions. Define the
total quartet expansion G=sum_t(|Gamma_t|-|S_(t-1)|). Exact counting gives

    G+|B0|=E+w,    E>=D,    |B0|<=D+w,    D<=E<=D+3m.       (1)

The last inequality uses |Gamma_t|-u<=3 min(u,512-u). A history close to
empty/full states therefore has coverage close to A+D, rather than an
arbitrary footprint of the same total activity.

## An all-history theorem

Except with probability less than 2^-76 over setup, no nonzero J-message
has output weight at most W and a realized support history satisfying

    m+32z-ceil(log2 binom(I,z)) >= 32d+7216.                (2)

A simpler sufficient condition is m+11z>=32d+7216, since I<=N and
binom(I,z)<=N^z. In particular m>=32d+7216 is sufficient without using
the output restriction.

The quantifier in (2) includes every message, support history, outer-row
set, and free terminal state. This is a restricted event over J, not a
claim that its support-class labels are preserved under projection of a
128-bit alphabet. A completed unrestricted J-distance theorem would
extend by the usual basis-component argument.

### Quartet expansion with exact integer verification

For a fixed incoming support of weight u, let G_t be the extra slots
covered by a uniform quartet partition. Its moment is exactly

    E[4^(-G_t)] = [r^u](256+16r+96r^2+256r^3+256r^4)^128
                         / (256^128 binom(512,u)).

All513 integer coefficient comparisons establish

    E[4^(-G_t)] <= 2^(-min(u,512-u)).                      (3)

This includes every incoming support; it is not a sampled spectrum or a
typical-state estimate. The partitions at distinct updates are independent.

By (1), q^(-E)<=4^w 4^(-G)(q/4)^(-D). The incoming minority mass is
m-min(w,512-w), and 2w+min(w,512-w)<=1024. Consequently

    E_partitions[q^(-E)] <= 2^1024 2^(-m) 2^(-30D).        (4)

### Sum all histories, including their multiplicity

For one binary coordinate over T=4096 updates, the number of paths
starting at zero and having r deaths is binom(T+1,2r+1), with the
terminal bit unrestricted. Set delta=2^-15. The sum of death weights is

    Z1=sum_r binom(4097,2r+1) delta^(2r)
       <=4097 exp(4096^2 delta^2/6).

Indeed, (2r+1)!>=6^r r!, and binom(T+1,2r+1)
<= (T+1)T^(2r)/(2r+1)!. Coordinates factor when summing these weights,
so the sum Z over all state-support histories satisfies

    Z=Z1^512 <2^6144 exp(1/8+4/3) <2^6147.                (5)

The last comparison uses ln2>2/3, verified by a rational Taylor bound.

### Apply the joint-rank bound before averaging

For fixed S, partitions, and route, the projectively normalized field
count from the earlier global proof is at most

    (gamma/q) binom(I,z) q^(d-E-z),    gamma<2.

On (2), the factor binom(I,z)q^(d-z) is at most 2^(m-7216).
Apply that inequality on the event indicator **before** dropping the
indicator and averaging partitions or routing. It is valid even though
d depends on the route. Equations (4) and (5) then bound the total by

    2^-31 * 2^(1024-7216) * Z < 2^(7140-7216) = 2^-76.

The chosen z is determined by S; there is no uncharged union over z.
An independent audit checked the setup conditioning, terminal factor,
history sum, and this order of inequalities.

## What remains outside the certified classes

A remaining bad history must violate (2), have F<53N/64, fail the earlier
rank thresholds, and contain a partial state. In particular, its shortened
outer dimension must satisfy

    32d > m+32z-ceil(log2 binom(I,z))-7216.

This is a concrete requirement for the remaining routing argument. It
does not imply that those geometries exist with significant probability.

Two additional counts are available. For fixed A,D, the number of state
histories is at most

    min{binom(N,A),
        min(binom(A+511,511),binom(N-A+511,511))
          *binom(A,D)*binom(N-A,D)}.                     (6)

For one coordinate of weight a and d deaths the exact count is
binom(a,d)binom(4096-a,d). Summing diagonal terms is bounded by the two
Vandermonde sums separately, proving (6). A simpler minority count is
2^4096 binom(N,m), obtained by rounding each state to empty/full and
recording its changed positions. These are upper bounds, not equivalent
descriptions of feasible trajectories.

Floating point checks combining these counts with the existing uniform
route moment and an expansion tail pass several dense partial-state
profiles. A wider3392-point diagnostic still fails. Its worst sampled
bounds range from about+33000 to+892000 bits; the failing points often
have E=D. Such point checks neither cover the full domain nor prove a
construction obstruction. They show that independently bounding the
history count, expansion tail, and equalized routing moment loses too
much of their joint geometry.

## Accepted maps and regional quotas

The companion diagnostic `quartet_accepted_rank_probe.py` retains both
superregularity tests before summing optional output supports. Sum over
fixed candidate row sets H, with h=|H|, and retain the event that H is
exactly the rows of covered degree at least129. On this event every
input outside H is zero, and

    R=128h+|L outside H|,

where L is the live-quartet part of V. In a quartet let U be the incoming
support, T its next support, Jout its output support, and Hq its H-slots.
A live quartet must satisfy

    |U|+|T union Hq|>=5,    |U|+|Jout union Hq|>=5.

Empty incoming quartets can have births only inside Hq. The same H is
used in every region; it is not chosen anew at each update.

For a Bernoulli surrogate with probability p', let P be its16-update
kernel and pi=Pr[Bin(8192,p')=h]. The quota-conditioned regional kernel
has the valid upper bound

    K_h <= E00 + pi^(-1)(P-E00).                         (7)

The empty regional path has weight exactly one for every label mask,
so it pays no conditioning penalty. Positive regional paths do. The
existing renewal identity sums paths visiting at least129 positive
regions. Entrywise minima of (7) over different p' remain upper bounds.

Additional necessary conditions are

    |V intersect H|>=129h,    |V outside H|<=128(8192-h).

Their valid scalar markers are implemented, but did not improve the
tested one-step potential. They are weaker than the individual degree
conditions defining H. For a fixed candidate H, condition on its binary
slot masks and V, before imposing the degree-classification event.
The remaining row assignments are uniform within each population.
Writing c_j,d_j for their covered counts in region j gives the stronger
classification filter

    Pr[Z_H>=129]^h * Pr[Z_outside<=128]^(8192-h),

where the two Poisson-binomial parameters are c_j/h and d_j/(8192-h).
The row negative-association coupling from the burst proof justifies the
powers. The two populations are independent. Endpoint populations are
handled separately. Evaluating this filter within the full history sum
remains open; it is not part of the proved theorem above.

## Verification and next discriminator

The exact513-comparison expansion check is available through
`quartet_global_coverage_probe.py --expansion-cert`; it also enumerates
short binary paths to check the activity/death and free-terminal counts.
The accepted kernel
passes its p=0 and p=1 controls, including the forced four-output support
from a singleton with zero input. Floating regional probes retain all513
boundary weights. They still fail: for h=1, the regional bound improves
from about+4258 to+224 bits at the tested markers. Four-tilt entrywise
minima give about+164k,+586k,+874k bits at h=512,2048,4096 respectively.
These are upper-bound failures, not distances or probabilities of bad codes.

All numerical work was serial, file-free on Peach, single-threaded,
memory-limited, and watchdog-bounded. No benchmark or construction change
was made. The full goal remains active.

Next retain temporal concentration together with the actual row-degree
classification in the global count. Before another large numerical cover,
test a dense burst with a persistent sparse fringe: its forced input rows
must acquire129 covered coordinates in the same global route. Repeating
scalar potential tuning or extending the parked local DP is not the next
step.
