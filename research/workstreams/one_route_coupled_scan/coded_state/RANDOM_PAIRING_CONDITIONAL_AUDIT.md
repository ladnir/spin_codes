# Conditional audit of the random-pairing proof

2026-10-26. Construction unchanged. No coverage search or native benchmark.

## Conclusion

Granting the current proposed flow bound does not close the distance proof.
Its threshold removes a positive dimension gap in an obstruction test.
At that threshold, the field exponent in the corresponding upper bound
is zero. Counting the support histories and output-zero choices still
costs exponentially many terms.

For one explicit half-full family, the direct union bound needs another
57,395 field-symbol units of coverage at the selected target-row count.
The causal-flow loss allowance would fall from 17.738% to 6.788%. The
existing mean-field diagnostic predicts about 12.34% loss. This comparison
is not a lower bound on actual loss or a counterexample to the code.
It does show why completing the present numerical target is not the next
priority: that target does not supply the missing end-to-end estimate.

Preserve the local lemmas, restricted certificates, and partial flow
coverage. First analyze the joint weighted count, instead of counting
every history against the worst permitted geometry. Do not change the
construction on the strength of this counting failure alone.

## Fixed construction and scope

Use the [random-pairing construction](RANDOM_PAIRING.md) over
q=2^32, with K=2^20, N=2^21, M=8192 RS[256,128] outer rows, n=512 state
coordinates, and L=4096 updates. The 256 regional permutations and all
perfect matchings form the geometry. Independently sampled accepted SL2
pair maps provide the remaining setup randomness. Setup is fixed for
every message. With update indices 1 through L,

    Q_0=0,   Q_t=A_t Q_(t-1)+X_t,   Y_t=Q_t+delta Q_(t-1).

The terminal state is free. The bad event is that some nonzero message
has output weight at most W=floor(0.46 N)=964689. This audit does not
establish that event's probability, a 49% certificate, or a new timing.

## A general conditional count

Fix the geometry and a proposed support history H_0,...,H_L, with H_0
empty. This is an arbitrary candidate history, not a claim that a
realized history is independent of setup. We will union over candidates.
Define

    A = sum_t |H_t|,
    D = sum_t |H_(t-1) minus H_t|,
    B = sum_t |H_t minus H_(t-1)|,
    I = sum_t |H_(t-1) intersect H_t|.

Then B=D+|H_L| and I=A-B. The union of consecutive supports has total
size A+D. Every birth and death forces a nonzero output. Therefore the
class is impossible if B+D>W. Otherwise a low-weight output requires
at least

    z=max(0,A+D-W)

zero-output positions among the I intersections. Choose z such positions;
there are binom(I,z) choices. The zero-output equations are independent
edges on coordinate paths, leaving A-z state-amplitude variables.

At update t, let F_t be all coordinates in pairs meeting H_(t-1). Let
U_t=H_(t-1) union H_t, J_t=F_t minus U_t, and V_t=F_t union H_t.
Write J=sum_t |J_t| and R=sum_t |F_t|. Then

    |V|=A+D+J.

Let p count the active next-state coordinates outside F_t, summed over
updates. These are births from empty pairs, including the initial state.
Their inverse inputs fix the first values of distinct coordinate-path
components. Thus the pin rank is exactly p, and

    R+p=|V|=A+D+J.

This is the [pin argument](UNIFORM_PAIR_PINS.md) with arbitrary matchings.
It allows changing weights, resets, and a nonzero terminal state.

Fix a nonempty exact set S of active outer rows. Let T be its complement,
r=|T|, and c_i the number of row i's positions in V. Two necessary
compatibility conditions are:

- Every J_t position routes to S. Such a position is an inactive endpoint
  paired with an active one, with zero next state. Superregularity forces
  its input to be nonzero.
- Every active row has c_i>=129, the RS minimum weight.

If either condition fails, this class has no witness. Otherwise shortening
the outer rows to V gives the upper bound q^d on candidate messages, where

    d=sum_(i in S)(c_i-128).

This counts messages with exact active set S by a linear-space relaxation;
it need not count that exact set sharply.

For each fixed candidate input, the pins and the selected zero outputs
leave at most q^(A-z-p) state assignments. For each valid assignment, the
independent accepted pair maps have joint point probability at most
2 q^-R. The factor two bounds their density loss over the entire code.
Relaxing the amplitude count does not apply a live-pair probability bound
to zero incoming vectors: only assignments with the stipulated supports
are used in that probability step.

Consequently, conditioned on the geometry, the probability that this
history and active-row set have a bad message is at most

    1[compatibility] min(1, 2 binom(I,z) q^E),
    E=d-J-D-z.                                               (1)

The cap follows the union over messages and zero selections. It precedes
the expectation over geometry. The all-zero message and S empty are
excluded; impossible support classes are discarded before relaxation.

Let C_T=sum_(i in T)c_i. The shortening identity simplifies (1) to

    E=A-K+128r-C_T-z.                                       (2)

In the regime A+D>W, this becomes

    E=-(K-W)+128r-(C_T+D).                                  (3)

Thus the bridge from coverage to field cancellation is available. The
unresolved issue is its sum over all compatible histories and row sets.

## What the granted flow bound would prove

The current selected target is r=1974. Since K-W=83887, its threshold is

    C_T+D >= 128*1974-83887 = 168785.

Suppose, hypothetically, that this holds simultaneously for every T of
that size and every history in the stated activity range, except with
setup probability epsilon. Equation (3) then gives only E<=0. Substituting
that uniform estimate into (1) gives a cap of one, not a small failure
probability. A stronger conclusion requires a positive reserve or a joint
count retaining more than this cutoff.

More generally, consider a family with at most H support histories and
at most C zero selections per history. Suppose a geometry event G,
depending only on route and matchings, guarantees compatibility implies

    C_T+D >= 128r-(K-W)+g

for every history in the family and every T of size r. Suppose also
Pr[G fails]<=epsilon. Independence of the remaining pair maps and the
conditional bound give

    Pr[bad message in this family]
      <= epsilon + min(1, 2 H C binom(M,r) q^-g).             (4)

This is a conditional theorem, not a proved premise about current setup.
The geometry event is simultaneous over histories, so its probability
does not need a second history union. The field-count term still does.
Conflating these two unions would incorrectly remove H from (4).

For a whole-code bound, sum disjoint family contributions or upper bounds
that cover them. A common geometry failure event can be paid once.
If separate geometry events are used, their failure budgets must be added.

## One explicit counting budget

Consider full-length histories with |H_t|=256 for every t>=1. At each
of the 4095 subsequent transitions, remove exactly 32 coordinates and
add exactly 32. This family has exactly

    H=binom(512,256) binom(256,32)^(2*4095)

support histories. Its parameters are

    A=1048576, D=131040, I=917280, z=214927.

The following intervals follow from exact integer counts and bit lengths.
They are not floating-point tail estimates.

| Factor | Interval containing its base-two logarithm |
|---|---:|
| Support histories H | [1109590,1109591) |
| Zero selections binom(I,z) | [720465,720466) |
| Target sets binom(8192,1974) | [6519,6520) |

For the direct bound (4), the smallest integer g satisfying

    2 H binom(I,z) binom(8192,1974) 2^(-32g) <= 2^-41

is 57395. This leaves a separate 2^-41 geometry budget for a total
2^-40 bound on this one family, not the whole code. The corresponding
coverage target is 168785+57395=226180. Integer rounding gives a certified
good-geometry log2 upper bound of -63 for this chosen g.

These are all combinatorial candidate histories, not histories shown to
be compatible with the outer or the sampled geometry. Their entropy can
therefore overestimate the cost of a sharper proof. The 57395 reserve is
the requirement of this particular direct bound, not a necessary condition
for distance and not a claim of low-weight codewords.

## Why completing the current numerical target is insufficient

The [causal flow](RANDOM_PAIRING_FLOW.md) gives, for every history,

    C_T+D >= (A+2(|H_1|-|H_L|)-Loss-3n)/4.

The explicit family has equal endpoint weights. To reach 226180 using
this inequality, it suffices to prove Loss<=142320. This is 6.788% of
the 4095*512 random-service coordinate updates. The current target is
Loss<=371900, or 17.738%.

The existing independent-partner approximation predicts 12.34% loss.
It would give coverage about 197313, above the old obstruction threshold
but below the 226180 target of the direct count. This approximation cannot
rule out a proof. It diagnoses a mismatch between the proposed certificate
and the intended composition. Finishing more histogram cells at the old
allowance does not remove that mismatch.

There are two distinct losses to investigate: the causal flow can
underestimate actual coverage, and (4) can overcount compatible histories.
An exact cut or a successful geometric tail alone does not address the
second loss. Conversely, weak causal allocation is not evidence that the
encoder itself needs stronger randomness.

## Remaining whole-code obligations

| Class or step | Audit status |
|---|---|
| Zero state history or empty active-row set | Only the zero message; exclude |
| Forced output weight B+D>W | Impossible; exclude |
| A+D<=W | z=0; exponent (2) still applies, but the selected high-activity flow premise does not cover this regime |
| A+D>W with activity below K | Equation (3) applies; the geometric guarantee is weaker and no complete substitute is certified |
| High activity, r=1974 | Current flow premise remains unproved globally; even granting it leaves the quantified counting gap |
| Other r, from 0 through 8191 | Some may admit simpler bounds; no complete partition or transfer from r=1974 is established |
| Arbitrary activation times and resets | Included in the general conditional lemma, not in all existing numerical certificates |
| Setup failure budgets | Must cover all geometry events used, without conditioning away dependencies |

The directed restricted results for R=192 and R=256, including up to
16384 swaps in half-full states, remain recorded in
[the active-row continuation](RANDOM_PAIRING_ACTIVE_ROWS.md). They do not
cover the example above or all the missing rows of this table. The cyclic
fixed-pair certificate is for a different construction and cannot fill them.

## Next discriminating step

Keep the code fixed. Start with the explicit 32-swap family and derive a
joint weighted sum retaining matching crossings, outer-row compatibility,
pins, and zero-output choices. Use (1) before replacing the geometry by
one uniform cutoff. First ask whether even an optimistic version of that
joint bound pays the roughly 1.84 million bits of the direct union count.
Do not launch an all-profile or all-histogram certificate until this
end-to-end test has a plausible negative exponent.

If the joint count still fails, identify whether coverage or admissible
trajectory multiplicity is responsible. Only then compare a small inner
change, such as occasional stronger mixing, with the preserved control.
This audit gives no basis for declaring random pairing insecure or for
claiming that simply increasing pair-coefficient randomness fixes the gap.

## Checks

`audit_random_pairing_composition.py` independently checks 2304 small
support/matching geometries, 21168 compatible exponent identities, and
169 zero-edge/pin rank instances. The general statements above follow
from the set and path arguments, not from extrapolating these toy tests.
It also evaluates the displayed full-size combinatorial counts exactly
and checks the minimal integer reserve on both sides of its threshold.

The serial run took 2.254 seconds, with a 60-second watchdog and fixed
problem sizes keeping working storage below 256 MiB. More than 20 GiB
was free before the run. It wrote no data files. Existing numerical
coverage, frozen helpers, production code, and the paper were unchanged.
