# Stronger pair randomness: what the proof actually gains

The next research model keeps the outer, route, state size, and reset. It
replaces only the two-shear distribution with the rejected-uniform SL2
distribution from [the bounded diagnosis](LAYERED_BOUNDED_DIAGNOSIS.md).
The old candidate and the whole-code GF32 control remain unchanged.

This change does more than improve a probability constant. It removes
fixed-input relations inside live pairs. The resulting compatibility rank
has an exact formula, so a two-state support count replaces the old
three-state component count and its exceptional labels.

One previously failing empty-heavy profile now has a complete restricted
certificate, including every active-row count. Other profiles remain open.
No new whole-code distance or encoding-time claim is made.

## Construction and probability space

Use K=2^20, N=2^21, 8192 RS[256,128] rows, and 256 independently permuted
regions. The state has b=512 coordinates over J=GF(2^32). At batch t,
pair coordinates whose indices differ in bit t mod 9. For each pair,
independently sample A from SL(2,J), rejecting unless both A and A+delta I
are superregular; delta is a fixed GF256 element other than zero or one.
All maps and routing choices are sampled once and fixed for every message.

The pair update is

    Y = X + (A+delta I)Q,
    Q_next = delta Q + Y = A Q + X.

The earlier local proof gives point probability at most
1/(q^2-7q+2) for the image of every fixed nonzero pair, where q=2^32.
Over the entire encoder, its inflation relative to q^-2 per live pair is
less than two. Thus inverse counting gives two random coordinates for
every nonzero pair and no random coordinates for a zero pair.

An accepted map uses three shears, rather than the previous two. This adds
one variable product per pair. It does not change the number of global
routes. The operation-count increase is not an end-to-end timing estimate.

## Why the old step was difficult

In the old distribution, some live pairs leave a fixed input coordinate.
Exceptional pairs leave two fixed coordinates. Those constraints can link
state values across time. A component can remain free, become pinned later,
or disappear. Counting such components separately from exceptional updates
loses their common dependence.

The stronger distribution has a simpler dichotomy:

- A live pair supplies two random inverse-input coordinates.
- A zero pair fixes X=Q_next on both coordinates.

The second case still contains useful information. If the next coordinate
is active, its input is nonzero and fixes the new state value. If it stays
inactive, its input is zero and constrains the outer word. Neither case
should be treated merely as a missing random coefficient.

## Exact pin identity

Fix a complete support history, including the terminal state, and a set
of zero outputs between active states. State values at active positions
are variables. Each zero output gives Q_next=delta Q along one coordinate
path. Let A_count be the number of active state positions and Z the number
of these edges. Their equations are independent, so the relaxed amplitude
space has dimension

    a = A_count - Z.

Nonzero inequalities are omitted when taking this dimension. They cannot
increase any subsequent count.

Let L map this space to the inverse-input coordinates that are fixed,
rather than random. Let H count initial active coordinates and subsequent
active coordinates born from an incoming zero pair. Then

    rank(L) = H,             dim ker(L) = a-H.

To prove this, consider the components of the temporal paths after keeping
the zero-output edges. Each fixed nonzero input names the value at a
component's first vertex. An initial active vertex has no predecessor.
A birth from a zero pair has an inactive predecessor, so its output cannot
be zero. Distinct such inputs therefore name distinct components. Every
other fixed input is identically zero. No component can receive a second
pin later while it remains active: its incoming pair would then be live.

This is precisely the complication removed by the distribution change.
There are no remaining linear relations between distinct nonzero fixed
input coordinates in this relaxed state space.

The test independently constructs GF4 constraint matrices for all 28,561
three-vertex support/zero-output histories on four coordinate paths.
Their ranks agree with both identities. These finite checks support the
implementation; the path argument establishes the general statement.

## Combining the pins with the outer

Consider a suffix of R regions, with zero incoming state. Put k=R-128.
Fix the exact set S of active outer rows, of size s. Let J_count count
fixed-zero inputs routed into S. A shortened active RS row has dimension
k; these zeros give the valid count bound q^(ks-J_count). If a row has too
many zeros to remain active, its actual contribution is zero.

Let D be the number of random inverse-input coordinates. On the retained
good-route event, let d_s bound outer shortening to those coordinates in
S. Counting states first gives q^(a+d_s). Counting compatible messages
first gives q^(a-H+ks-J_count). Both bounds concern the same history.
For any theta in [0,1], their minimum is at most

    q^[a + theta*d_s + (1-theta)*(ks-J_count-H)].

Multiply by the inverse point probability, at most 2*q^-D, and require
every fixed-nonzero input to route into S. This uses neither independence
between constraints nor a uniform-state assumption.

For the implemented upper bound, omit terminal support and fixed-input
constraints. If V counts the noninitial incoming active coordinates,
then a<=V+b-Z, where Z now counts only internal zero-output edges. Retain
only the initial pins and internal births from zero pairs in H. These
omissions enlarge the bound and leave the terminal state free.

For fixed S, its routed positions form an independent uniform s-subset in
each region. As in the retained proof, use Bernoulli colour probability p
and divide by the exact per-region conditioning probability

    pi_s(p) = binom(8192,s) p^s (1-p)^(8192-s).

Write theta=m/32 and beta=1-p+p*2^(-(32-m)). Local factors are now:

| Event | Factor in the weighted support sum |
| --- | --- |
| Internal zero output between active states | 2^-32 |
| Birth from an empty pair | p*2^(-(32-m)) |
| Coordinate staying zero after an empty pair | beta |
| Initial active coordinate | p*2^(-(32-m)) |
| Initial inactive coordinate | beta |

Markers additionally record empty pairs, full pairs, births, and output
weight. The zero-output factor is always 2^-32, not 2^-m: the exact pin
identity retains the same zero-output rank in both amplitude counts.

The resulting local matrix has four pair-support states: 00,01,10,11.
Positive coordinate potentials give pointwise row bounds which compose
along the actual cycling pair schedule. A two-update variant retains
crossing pairs on four coordinates; it changes the analysis, not the code.
An independent literal sum checks 64 exact two-update contractions.

## Restricted certificate and remaining failures

For the stronger-pair ensemble, the following aggregate profile closes:

    R=192, P=(R*8192-512)/2=786176,
    empty pairs z=393088, full pairs f=294816,
    internal births B=98272.

All orientations of its remaining 98,272 one-sided pairs are included.
The actual matching-moment envelope gives d_s<=min(64000,64s). Four frozen
positive-weight proposals are replayed as rational inputs. Directed integer
arithmetic checks every s=1,...,8192 and includes all choices of S.
Adding 13 bits for the union over s gives:

| Output threshold | Log scale 2^24 | Log scale 2^28 |
| --- | ---: | ---: |
| floor(0.40 N) | -5716 | -5717 |
| floor(0.49 N) | -5716 | -5717 |

These bound the bad-word contribution intersected with the good-route
event, not an unnormalized conditional probability. The route failure
budget must be added once in a future whole-code proof. The equality of
these two results occurs because the retained weights barely use the
output-weight restriction for this particular profile. It does not imply
that other profiles reach 49%.

The earlier coarse ideal-pair bound for this empty-heavy example was
positive (+478536). The new certificate concerns the stronger-pair model,
not the old two-shear construction. Partial base-field history certificates
do not yet yield a whole-code certificate over GF(2^128).

Two other retained families remain unresolved at 40%:

- R=160, z=f=327552, B=81888. At s=4096,m=28, the one-update bound was
  about +168436 before explicitly excluding absent pair types. The
  two-update bound with no one-sided pairs is +23047. At m=27 it improves
  to a directed bound of -218. This closes that source-count slice, not
  the entire profile; all-source coverage remains to be checked.
- R=256, z=f=524160, B=0. At s=4096,m=28, two-update counting still gives
  about +3798111. Short lookahead does not resolve this example.

## A structural explanation of the second failure

The last profile has no births and no one-sided pairs. Its supports must
decrease with time, and each incoming support is a union of pairs in the
current matching. A proper nonempty support cannot obey both requirements
through nine distinct pairing phases.

For a proof, choose a coordinate active at the last of those nine phases.
Closure under the last matching puts its partner there too. Both were
active earlier because supports only decrease. Closure under the preceding
matching doubles this set again. Continuing backward forces all 512
coordinates into the initial support, contradicting its being proper.

Consequently this history consists of a full-state prefix, at most eight
partly active incoming states, and an inactive suffix. The terminal update,
whose constraints the current count omits, is handled separately. Exact toy
enumeration checks all decreasing paths through the distinct phases at
widths 2,4,8,16; the same backward argument applies at width 512.

The coarse count hides this temporal concentration. It permits the outer
zeros and the partial-state freedom to be arranged much more favorably
than this history allows. Adding another random map cannot fix that
particular counting loss. On the full-state prefix, every pair is live,
so the stronger-pair model already supplies essentially full-batch inverse
randomness, with the established total density loss below two.

## Next discriminating step

Keep this three-shear research model fixed for one episode-aware pass.
Retain the full-state prefix, short partial tail, and inactive suffix in
the actual regional routing count. Start with the no-birth/no-singleton
family above, then charge births and one-sided pairs as explicit defects.
The uniform-feedback control is a reference for that argument, not a
certificate that can be imported without paying for the partial tail.

Only if partly live, defect-rich histories remain difficult should we
compare a second crossing layer or a periodic stronger refresh. That test
must show which new rank or point-probability fact it buys. Do not increase
randomness just to compensate for a proof that discarded where zeros lie.

No encoder, outer, global route, paper, frozen helper, or previous proof was
modified. No native timing was run. This round used serial local Python
checks, small fixed arrays, bounded proposals, and no generated data files.
The C: volume had over 20 GiB free before numerical work.

Reproduce from this directory, serially:

    python -B diagnose_uniform_pairs.py
    python -B diagnose_pair_pins.py --test
    python -B diagnose_pair_pins.py --tensor-test
    python -B diagnose_pair_pins.py --geometry
    python -B diagnose_pair_pins.py --replay
    python -B diagnose_pair_pins.py --replay --distance 490

The proposal command is optional and does not replace the frozen replay:

    python -B diagnose_pair_pins.py --case 0 --cover

Proposals have per-call time caps; the cover uses at most 14 proposals
and a 90-second loop budget. The deterministic selected-profile replay has
a 45-second cap. Only compact source inputs and these results are retained.
