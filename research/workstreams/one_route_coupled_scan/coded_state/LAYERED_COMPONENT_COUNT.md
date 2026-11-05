# Joint component and exceptional-update count

The construction is unchanged: RS[256,128], the regional route, and the
layered pair-shear inner over J=GF(2^32). The GF32 sandwich remains the
whole-code control. This note closes the previously hardest retained
mixed-history class at the exploratory 40% threshold. It does not close
all histories or establish a near-half distance theorem.

The improvement is not a stronger worst-case bound on the component count.
We sum components and exceptional updates together, then interpolate two
valid amplitude bounds. A short lookahead in the existing pairing schedule
prevents incompatible support patterns from inflating that sum.

## Class and probability space

Put M=8192, b=512, N=2^21, and q=2^32. Consider the suffix of R=192 regions,
with u=RM, L=u/b batches, and T=L-1=3071 noninitial incoming states.
The state before the suffix is zero. The terminal state is free.
After the first batch, every pair is either type10 or type11; there are
no00 or01 pairs. Thus the first coordinate of every pair is active.

The class fixes these aggregate counts:

    P=(u-b)/2=786176 pairs;
    f=589632 full pairs;
    e=368520 exceptional full pairs;
    B=49136 newly active coordinates at internal transitions.

Output weight is at most W=838860=floor(.40N). All support placements,
exceptional placements, amplitudes, and internal output-zero masks with
these counts are included. This is not a fixed-message experiment.

Fix any route satisfying the existing simultaneous head/edge envelopes.
For this class, they imply d<=d_bar=163923 for every inverse free set F,
where d is the outer shortening dimension on F. The checker recomputes
this bound rather than accepting it as an input. Probabilities below are
over the independent inner setup, conditional on that route. The existing
good-routing event has its own failure budget; it is not resampled for
different messages or histories.

## Count both amplitude constraints on the same history

Write V=P+f=1375808, D=P+f-2e=638768, and K_R=u-N/2=524288.
Fix a support/exception/output-zero history H. Let Z(H) count internal
zero-output edges whose two adjacent states are active. These equations
are independent along temporal coordinate paths. The state-amplitude
space therefore has dimension a=V+b-Z(H), before nonzero inequalities.

Let kappa(H) be the number of unpinned temporal components from the
[fixed-coordinate consistency lemma](LAYERED_FIXED_CONSISTENCY.md).
That lemma gives the following probability contribution for H:

    eta * q^min(a+d_bar-D-e, kappa(H)+K_R-D-e),    eta<2.

For any theta in [0,1], min(x,y)<=theta*x+(1-theta)*y. Choose theta=27/32.
After moving history-independent terms outside the sum, the class bound is

    2^[1+27*(V+b+d_bar)+5*K_R-32*(D+e)]
       * sum_H 2^[5*kappa(H)-27*Z(H)].

The sum still obeys the fixed f,e,B counts and the output-weight bound.
Both amplitude arguments are applied to the same history. This does not
assume their failure events are independent. It also avoids substituting
an upper bound on d into a claimed lower bound on an extra rank saving.

## Three states suffice for each coordinate path

At a temporal vertex, retain actual support and whether its auxiliary
kernel component is still unpinned. There are three possibilities:

| Index | Actual support | Current component |
| --- | --- | --- |
| 0 | inactive | pinned |
| 1 | active | pinned |
| 2 | active | unpinned |

An unpinned component contributes a factor 2^5 only when it ends. A later
pin can otherwise destroy it, so charging that factor when a component
starts would overcount its surviving dimension.

For each fixed input coordinate, X=Q_next+gamma*Q with gamma in {0,1}.
For a free coordinate use the tag gamma=-1. Given an actual output-zero
flag, the update is deterministic:

- A zero output or gamma=1 links adjacent temporal vertices.
- An inactive next vertex or gamma=0 pins the next vertex.
- A zero output together with gamma=1 pins the linked component, since
  delta!=1.
- An unlinked edge finalizes the old component and starts the next one.

The two endpoints must share a branch choice. A type10 pair has tags
(1,-1). A regular full pair has (-1,-1); an exceptional full pair has
(0,1). The count includes one exceptional marker per pair, not independent
markers that allow the endpoints to disagree.

The first state has all auxiliary components pinned because X_0=Q_1 is
fixed. At the final update, omit output constraints and births. A free
terminal coordinate contributes one new component; each old unpinned
component is finalized. This bounds the free terminal state without an
unjustified return-to-zero requirement.

## Positive weights and pairing lookahead

Use positive markers m_f,m_e,m_B,m_W, with 0<m_W<=1, for full pairs,
exceptional pairs, births, and counted output weight. Counted weight is
the first output plus internal outputs; the final output is omitted.
It is therefore at most W for every target history. Removing the markers
costs m_f^-f m_e^-e m_B^-B m_W^-W.

At time t, pairing uses bit t mod9 of the coordinate label. Attach a
positive potential h_(a,c)(s) to coordinate state s, where a and c are
the label's current and next pairing bits. If a=0, actual support must
be active because this class contains neither00 nor01 pairs.

For one current pair, both endpoints share the next two label bits (c,d).
Its input potentials are h_(0,c) and h_(1,c); both output potentials are
h_(c,d). Form the exact local weighted sum, including the shared branch,
component completion, zero rank, births, and output weight. Divide by the
input potentials. Let Lambda_(c,d) be the maximum of this ratio over the
six allowed incoming pair states. This gives a pointwise inequality,
not a claim that neighboring pairs or successive states are independent.

Within each batch, pairs are disjoint. Multiply their local inequalities;
then compose batches. Each (c,d) occurs in exactly b/8 pairs because the
three pairing bits are distinct, including at the nine-step wraparound.
Thus every Lambda_(c,d) occurs (P-b/2)/4=196480 times.

For c in {0,1}, define the initial factor

    I_c = h_(0,c)(1)*m_W * (1+h_(1,c)(1)*m_W).

Each occurs b/4=128 times. The exact terminal local sum, divided by
h_(0,c)h_(1,c), has maximum E_c; each E_c also occurs128 times.
The resulting complete upper bound on the class sum is

    (I_0 I_1)^128
      * product_(c,d) Lambda_(c,d)^196480
      * (E_0 E_1)^128
      * m_f^-f m_e^-e m_B^-B m_W^-W.

This counts all histories in the stated class. It retains enough of the
pairing schedule to anticipate future forced-active coordinates without
enumerating512-coordinate support masks.

## Exact replay and checks

The proposal code uses floating-point optimization only to choose markers.
The verifier stores twelve explicit positive rational markers/potentials.
It evaluates every local row with fractions, takes exact maxima, and uses
integer powers and bit lengths to bound logarithms. For x=n/d>0,

    log2(x) < [bit_length(n^p)-bit_length(d^p)+1]/p.

Apply this to every factor with a positive exponent, using reciprocal
markers for the coefficient extraction. No floating logarithm or optimizer
success flag enters the verified result.

| Replay | Conditional log2 failure upper bound |
| --- | ---: |
| Integer-power precision p=256 | -38885 |
| Integer-power precision p=512 | -44123 |

The earlier bound for the same class was +734997. The intermediate
coordinate-independent component sum did not close it. Sharing pair
branches, interpolating the amplitude bounds, and retaining only the next
pairing orientation gave floating upper bounds about+58394 and+9334.
Retaining two future pairing bits yielded the successful proposal above.
These intermediate numbers are diagnostic, not certified margins.

Exhaustive checks enumerate650 width-four histories with two noninitial
incoming states and a free terminal state. The transfer's component count
agrees with the independent path-component oracle on every history.
An exact weighted sum over these histories satisfies the positive-potential
bound, including the initial and terminal factors. Previous GF4 matrix-rank
checks remain applicable to the component oracle.

The same fixed rational markers also pass27 selected neighboring profiles:
f offsets{-1024,0,1024}, e offsets{-1024,0,1024}, and B offsets{-512,0,512}.
Their route envelopes are recomputed. The worst p=256 upper bound is-14603.
These checks are not interval coverage and do not cover profiles with01
or00 pairs.

Run serially from the repository root:

    python -B research/workstreams/one_route_coupled_scan/coded_state/check_layered_component_count.py

The bounded proposal, if needed, is:

    python -B research/workstreams/one_route_coupled_scan/coded_state/layered_component_transfer.py --lookahead 2

No encoder operation, permutation distribution, native timing, or production
code changed. Next apply this count to the other retained failing profiles,
then generalize the local tables to01/00 pairs before claiming full coverage.
The now-closed class is at40%, not a certificate for the near-half target.
