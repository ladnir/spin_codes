# Grouping state histories inside the covered footprint

The earlier global proof takes a fractional moment separately for each
state-support history. Histories with the same covered input footprint
share the same outer parity rank, so that separation pays too much for
their multiplicity. This checkpoint sums those histories first. It gives
a new global bound for the unchanged construction; its full numerical
evaluation and the 45%/40-bit certificate remain open.

Use J=GF(2^32), q=2^32, width512, T=4096, N=2097152, and W=943718.
The probability space and accepted setup are those of
[the global rank reset](QUARTET_GLOBAL_RANK_RESET.md).
The terminal state is free. Here V_t denotes covered positions at update
t, rather than state support. Each candidate sequence V is fixed before
exposing setup. It is tested for being the actual footprint later.

## Summing the state histories first

Fix V, quartet partitions, and the route. Its outer parity rank is
R(V)=sum_i min(128,|V intersect row i|). Let A_t consist of positions in
V_t whose quartet at update t+1 is wholly contained in V_(t+1).
For the terminal update set A_T=V_T, and put L=sum_t|A_t|.

Every compatible state support S_t is contained in A_t. Otherwise an
active state position would force coverage of a quartet containing an
uncovered position at the next update.

Choose 0<x<=1. For each genuine compatible state history, the joint-rank
proof bounds its weighted projective field count by

    (gamma/q) q^(-R(V))
      (qx)^births x^deaths (1+qx)^survivors,    gamma<2.   (1)

All genuine histories producing V have the same R(V). Sum (1) over them
before taking a fractional moment. Only then enlarge the sum to arbitrary
states S_t contained in A_t. The added histories need not produce V;
they are nonnegative terms in an upper bound, not additional applications
of the rank lemma to the wrong footprint.

For one coordinate the allowed-state transfer is

    H = [[1, qx], [x, 1+qx]],

and a step forbidding a nonzero state uses H0=[[1,0],[x,0]]. The positive
vector v=(1,1+1/q) satisfies

    H v <= lambda v,    H0 v <= v,    lambda=1+(q+1)x.

The initial state is zero and the terminal vector (1,1) is at most v.
Multiplication of these inequalities, separately on each coordinate,
bounds the whole state-history sum by lambda^L. No terminal-state reset
or independence between histories has been introduced. In particular,

    Pr_field[bad output with footprint V | partitions,route]
      <= min(1,(gamma/q) x^(-W) q^(-R(V)) lambda^L).      (2)

For 0<theta<=1, put a=q^theta, z=x^theta, and b=lambda^theta. Average
the theta power of (2), then sum over nonempty candidate footprints.
For each fixed V, rank depends only on routing and L only on partitions;
their geometry expectations factor. This gives

    (gamma/q)^theta z^(-W)
      sum_(V nonempty) E_route[a^(-R(V))] E_parts[b^L].  (3)

This order moves the state-history entropy inside the fractional moment.
The bound includes an overcount by the zero state inside a nonempty
candidate footprint. Removing it before the moment is possible, but
(X-1)^theta must not be replaced by X^theta-1 when theta<1.

## A footprint transfer with one boundary count

In a quartet containing u incoming covered positions, an outgoing covered
mask of size v contributes b^u if v=4, and one otherwise. Thus the exact
quartet polynomial is

    F(r,y)=(1+r)^4((1+y)^4-y^4)+(1+br)^4 y^4.

The width512 transition is [r^u y^v]F(r,y)^128 / binom(512,u).
Start at the empty footprint, compose T transitions, and use terminal
weight b^v. The terminal weight supplies A_T; it must not be dropped.
A total-coverage marker replaces y by xy, with a separately named marker
when the output marker x is also in use.

The scalar exponential-potential diagnostic still fails at all eight
tested coverage fractions between1/8 and4/5. It averages away temporal
coverage information and permits unrealizable footprints. The failed
envelope is not a lower bound on code failure.

## Retain forced births before composing

There is a simple necessary realizability condition. A partially covered
quartet cannot have an active incoming state: that would cover the whole
quartet. Every covered position in such a quartet is therefore a forced
birth. Let B_t be this set. At the first update, B_1=V_1, including fully
covered quartets, since the initial state is zero.

Compatibility requires B_t subset A_t. Set f=sum_t|B_t|. At a forced
birth the transfer is Hb=[[0,qx],[0,0]], which satisfies
Hb v<=c v for c=(q+1)x. The stronger state-history bound is consequently

    c^f lambda^(L-f).                                  (4)

After taking the moment, write m=c^theta and b=lambda^theta. A covered
position is now either forced or optional. At the next update, every
quartet meeting a forced position must be wholly covered. This gives the
quartet polynomial

    y_o^4(1+m r_f+b r_o)^4
      +(1+r_o)^4(1+4y_f+6y_f^2+4y_f^3).                 (5)

Here r_f,r_o mark incoming forced/optional positions, and y_f,y_o mark
outgoing forced/optional covered positions. Uncovered positions have no
marker. A full outgoing quartet contributes four optional positions;
a partial one contributes only forced positions. The first footprint is
entirely forced, and the final vector is m^f b^o. Take the128th power and
divide by the incoming multinomial coefficient to average the actual
uniform quartet partition.

Equation (5) preserves forced propagation, but still relaxes the
requirement that a full covered quartet with zero incoming state must
activate all four outgoing state positions. It is an upper counting
interface, not an exact characterization of realized footprints.

There is a further safe improvement: when a quartet has no incoming
covered positions, a full outgoing footprint is also entirely forced.
Replace the full-quartet term in (5) by

    y_f^4+y_o^4((1+m r_f+b r_o)^4-1).                   (7)

This handles every first activation, including one after a zero prefix.
With (7), start at the empty footprint and use the same kernel at all T
updates. No special first-update distribution is needed. In contrast,
the weaker labeling in (5) can have f=0 for a nonzero history; it cannot
be used to assert the existence of a labeled forced birth.

The helper evaluates (5), not (7), through a coarse bound using the largest quartet row under
a product potential. Those eight probes also fail. The full averaged
two-count transfer in (5) has not been composed through all regions.
Do not interpret the largest-row probes as that composition.

## Routing variance should use coverage, not state holes

Let v_t=|V_t|, and define

    Q=sum_t v_t(512-v_t)/512,
    D_V=sum_(t<T) (v_t-v_(t+1))_+.

For actual histories, D_V<=E, where E is covered footprint minus state
activity. If p_j is coverage fraction in region j, then

    sum_j p_j(1-p_j) <= D_V/1024+1/2+Q/8192.             (6)

To prove (6), round every v_t/512 at a common uniform threshold. The
average number of closed runs in the resulting binary sequence is
D_V/512. The burst bound gives expected row variance at most
D_V/1024+1/2. The gap between this expected variance and the variance of
the mean regional coverage is a sum of threshold variances. Within each
16-update region, variance of an average is at most the average of the
individual variances; their sum is Q/8192. Finally,
v_(t+1)>=|S_t|=v_t-E_t implies D_V<=sum_t E_t.

A dense covered plateau has Q=0 even if the state itself has many holes.
Those holes still pay field suppression. A bound based only on state
minority incorrectly charges them as routing variance. This is a
structural distinction, not a change of construction.

## Status and next step

Two independent audits checked (1)--(3), the grouping order, and the
footprint kernel. The strengthened positive-vector bound is exact.
The forced-birth bound, three-label kernel, and improvement (7) also
passed independent mathematical audit.
The helper `quartet_grouped_coverage_probe.py --checks` passes504 exact
allowed-mask chain checks,106 independently enumerated footprint-kernel
checks at widths4/8, and1024 exact two-update forced-birth checks.
Numerical envelopes remain floating diagnostics and do not close.

Next compose the forced-birth footprint count with actual regional
coverage and the outer row-degree constraints. Test a dense burst plus
sparse fringe before attempting a global numerical cover. This is the
new global reduction; the earlier state-history expansion theorem remains
proved and the old local seven-type refinements remain parked.
