# Fast high-distance SPIN: first geometry screen

No new lower-expansion candidate passed the tested 40% bounds. The existing
GRS[256,64], width-four/state-four, keep-two construction remains the positive
proof control, not a solution to the 3--5 ms forward-encoding target. The new
screens rule out spending this iteration on a native state-three or wider
3K prototype on the strength of these bounds alone. They do not prove that
those codes have low distance.

## Construction and scope

The requested source dimension is K=2^20 and the final nominal rate is 1/2.
Every candidate here uses a fixed GRS[b,k] outer over the embedded GF256,
independent uniform row permutations within its b coordinate regions, and
the existing normalized full-Cauchy inner. All input and state diagonals are
independent uniform nonzero GF(2^16) scalars. There is one regional routing
stage. A 128-bit record contains eight parallel scalar lanes; the claimed
target is symbol/record distance, not bit distance.

The notation `keep(5,5,6)` means that consecutive width-eight packets retain
five, five, and six fixed output coordinates, periodically. Packet order is
preserved in the transfer product. It is not an average of three transfers.
Rows are rounded to complete packet periods, so nominal padded rate 1/2 can
differ slightly from requested K divided by output length. For example,
GRS[255,85] has 12,480 rows, padded source dimension 1,060,800, and output
length 2,121,600. Padding constraints can only improve distance relative to
the full padded source space; they must still be stated by the implementation.

The values below are floating fixed-witness upper-bound margins for one
active-row count A or one raw-input weight near density .03. Positive means
that this bound is small for that class. Negative means that this bound is
uninformative. Neither sign is a whole-code certificate or an observed
codeword. Numerical searches used at most 55 L-BFGS-B iterations per start;
some stopped before convergence. The saved witness remains reproducible
without relying on convergence or rerunning the optimizer.

| Raw expansion; outer | Width / state; retention | A=16 | A=64 | Density .03 |
|---|---|---:|---:|---:|
| 2K; [256,128] | 4 / 4; keep4 | -2515.7 | -- | -220206.8 |
| 2K; [256,128] | 4 / 8; keep4 | -955.2 | -- | -160191.4 |
| 4K; [256,64] | 4 / 2; keep2 | -264.4 | -- | -23583.6 |
| 4K; [256,64] | 4 / 3; keep2 | 580.5 | -98.4 | -5969.7 |
| 3K; [192,64] | 4 / 8; keep(3,3,2) | 1349.3 | -- | -37379.6 |
| 3K; [192,64] | 8 / 8; keep(5,5,6) | 1350.6 | -- | -35164.8 |
| 3K; [255,85] | 8 / 8; keep(5,5,6) | -- | -2009.6 | -30259.4 |

The 3K width-eight design was a concrete traffic-reduction candidate:
25% fewer raw records than the 4K control, larger outer blocks that still
fit GF256, and XOR-circulant fixed maps. It failed the tested bound in both
an intermediate active class and the low-middle class. Larger circulants
also have more arithmetic per raw record; raw expansion alone would not
establish a speedup even if its proof screen passed.

## What is an obstruction, and what is not

The [regional window argument](HIGH_DISTANCE_WINDOW_OBSTRUCTION.md) is a
real coefficient-independent obstruction. If fixed windows leave

    D = sum_i max(0, hits_i - (b-k)) > b*s

shortened source dimensions, at most b*s exit equations cannot annihilate
all those sources. A surviving nonzero source has output supported inside
the windows, for every fixed assignment of local linear maps. This remains
valid with periodic puncturing when the windows contain complete periods:
their retained-output fraction equals their raw-input fraction. A zero
punctured image also fails the requested injective encoding contract.

The new helper applies this exact fixed-window lower-bound gate as well as
replaying numerical witnesses. A zero lower bound excludes nothing: variable
window lengths, multiple windows, and connected boundary windows still need
their own treatment. The [realized-window fixtures](WINDOW_DIMENSION_SEARCH.md)
show why arbitrary local state matrices cannot repair a bad routed geometry.
Replacing regional routing by one full uniform permutation is not a repair;
the [global suffix obstruction](HIGH_DISTANCE_GLOBAL.md) is stronger there.

The zero-feedback ceiling in [HIGH_DISTANCE_GATE.md](HIGH_DISTANCE_GATE.md)
applies to its stated identity-feedthrough inner. It should not be silently
transferred to a differently normalized or punctured full-Cauchy map.

In contrast, all negative entries in the table are upper-bound failures.
They do not establish a quiet-gap witness, a setup-failure lower bound, or
an impossibility result for those geometries.

## Precise remaining proof direction

Increasing the coefficient field alone does not address these witnesses.
For example, the 2K A=16 witness uses rank tilt exp(2.0115), approximately
7.47, far below even 255. Adding the existing minimum-source-weight mark
to that optimization returned its neutral value exactly; its margin stayed
-2515.6615 bits. A total-weight lower bound is therefore not the missing
information in this tested relaxation.

The present positive transfer retains only actual-zero versus nonzero state.
It separately bounds support shortening and then takes a worst incoming
state rank at each packet. A new lower-expansion proof should retain more of
the joint source/state history rather than count every bad message: a rare
routed subspace can contain many bad messages, so that first moment need not
bound the relevant setup event sharply.

One concrete next proof object is a full-Cauchy rank-cap transfer. After
selecting z output zeros at a packet with w active inputs, let p=min(w,z)
be the input pivots. Eliminate their fixed Cauchy columns, then analyze the
remaining z-p state equations before exposing the corresponding fresh input
diagonals. If the surviving incoming state rank is H, the proposed successor
cap is min(s,H+w-p). Births start with rank at most min(s,w), rather than
the unrestricted nonzero state. Returns should retain the incoming cap in
their zero-flat count. The existing
[rank-history audit](GF8_RANK_HISTORY_AUDIT.md) establishes this style of
composition for a different identity-feedthrough map, not automatically for
this full-Cauchy packet.

The exact gap is to prove the Schur-complement exposure order, the joint
return flat bound, and preservation of the nonnegative continuation cone
for the full-Cauchy map and each retained-output count. Independent maxima
for successor-rank masses are not valid; simultaneous upper-tail domination
is required. Until those facts are proved and the low-middle deficit closes,
this is a research direction rather than a candidate certificate. It may
still be too weak: earlier rank-history refinements in another geometry
gave hundreds, not tens of thousands, of bits.

The sharper source alternative is a support/cut count that tracks rowwise
shortened dimensions jointly with state ranks. It must control the rare
window subspaces explicitly and account for unequal regional occupancies.
Replacing the existing feasible-support polynomial by a narrower polynomial
without a valid homogenization argument would not suffice.

## Replay and next decision

The new helper changes no existing proof or native implementation. It stores
the exact integer classes and floating marks used above, and offers a bounded
replay with numerical regression checks:

```text
python -B research/workstreams/char2_field_spin/theory/hd_fast_screens.py --check --windows
```

Only lightweight local proof screens were run. No native build, benchmark,
remote command, full parameter sweep, or directed certificate was run.
The next iteration should either prove and screen the rank-cap/cut refinement
on the 3K [255,85] geometry, or explicitly choose a different mixing topology.
It should not promote a 2K/3K implementation or the 4K state-three variant
as a 40%-distance design on the present evidence.

## Expanding inners: a deterministic finite window ceiling

The later [global-prefix strengthening](#global-prefix-strengthening)
removes the regional-routing assumption and reduces the exit cost to one
state. The original regional calculation below is retained for comparison.

A rectangular inner can emit more symbols than it consumes. This reduces
raw traffic at fixed final rate, but it does not remove the causal window
argument. The following finite version needs no random-routing assumption.

Fix R source rows of a linear MDS[b,k] outer. Each of its b coordinate
regions contains one symbol from each row, in an arbitrary permutation.
The inner is causal, linear, zero-preserving, starts at zero, and carries
s field symbols of state. It consumes w inputs per packet. Its output count
has a fixed period of t packets, with positive counts m_1,...,m_t. Assume
w*t divides R. Choose the first W positions in every region, with w*t
dividing W. These windows retain exactly

    b * W * (m_1+...+m_t)/(w*t)

output coordinates. Their output fraction is W/R, even when some m_j>w.
Write p=b-k. Every route has exactly bW total window hits. Therefore its
shortened source dimension satisfies

    D = sum_i (hits_i-p)_+ >= bW-pR.

Impose zero state after each window; these are at most bs homogeneous
equations. Any independent source restriction of codimension c, including
padding constraints, costs at most c further dimensions. Consequently

    bW - (b-k)R > bs+c

guarantees a surviving nonzero source with output weight at most the
displayed window count, for every route and every local coefficient choice.
If its output is zero, encoding loses injectivity instead. Thus choosing
W/R below the requested relative distance gives a deterministic failure.
Charging all b exits is conservative, including at the global terminal end.

In particular, let W be the largest multiple of w*t strictly below delta*R.
Since W >= delta*R-w*t, a necessary condition for escaping this test is

    delta <= 1-k/b + s/R + c/(bR) + w*t/R.

This is only a necessary condition, with alignment slack. For fixed packet
width and state it approaches the outer redundancy fraction 1-k/b.
At final rate 1/2, targeting 40% therefore requires raw expansion at least
approximately 5/3, even if the inner expands its input. Larger local
coefficient fields do not change the argument.

Two proposed examples fail deterministically at K20. The table first uses
the full padded source space. Subtract the actual padding codimension from
the final column for the requested source restriction.

| Outer; packet output counts | R | W | bW-(b-k)R | bs | Surviving dimension at least |
|---|---:|---:|---:|---:|---:|
| [192,128]; w=4,s=4, keep(5,5,6) | 8256 | 3300 | 105216 | 768 | 104448 |
| [256,192]; w=4,s=4, keep6 | 5504 | 2200 | 210944 | 1024 | 209920 |

Both have padding codimension 8192, so substantial kernels survive every
padding convention imposing that many independent linear restrictions.
Their window fractions are strictly below .40.

### A non-deterministic obstruction just below the ceiling

GRS[256,144], w=8,s=8, keep9 has raw expansion 16/9 and nominal rate 1/2.
Its redundancy fraction 7/16 is above .40, so the deterministic inequality
does not reject it. Independent uniform regional permutations nevertheless
give a substantial failure lower bound for the same fixed-window event.

Here R=7296 and W=2912. Exclude e=18 complete source rows, conservatively
covering the four-row-interleaved padding layout. This restriction is valid
when the remaining complete rows belong to the requested source coordinates;
unrelated padding layouts must be checked separately. Every retained row
has hit count B distributed as Binomial(256,W/R). Thus the retained
shortened dimension S satisfies

    E[S] = (R-e) E[(B-112)_+] = 2901.4074876777...,
    S <= kW = 419328.

No independence between different row hit counts is asserted. Since the
exit capacity is 2048, the elementary bounded-variable inequality gives

    Pr[S>2048] >= (E[S]-2048)/(419328-2048)
                > 0.00204516.

The helper evaluates this as an exact rational number; decimals are only
display values. This is a genuine setup-failure lower bound, at most about
8.94 bits of reliability, not a failed upper-bound screen. Dropping eighteen
rows reduces the expectation, so this conclusion does not rely on the padded
source dimensions. Increasing state to sixteen removes this particular test
but would require separate cost and whole-code analyses.

### The next rectangular point also fails the current relaxation

GRS[256,136], w=16,s=8, keep17 has raw expansion 32/17, R=7744, and output
length 2,106,368. Its fixed-window expectation is 243.4556 against exit
capacity 2048: this particular lower bound is inconclusive, not a safety
certificate. The new helper `theory/hd_fast_rectangular.py` extends the
selected-zero full-Cauchy formulas to seventeen emitted coordinates and
screens this one geometry. A fixed tall Cauchy map uses outputs+state rows
and width+state columns; normalizing its bottom state block preserves the
ranks of selected output rows together with all state rows.

The initial GF(2^16), 40% samples give A=16 margin -1380.0058 bits and
raw density .03 margin -80572.8663 bits. The A=16 optimizer converged;
the middle search was bounded at fourteen iterations. The extended local
composition proof has not been independently audited, and no dense inverse-
membership argument is claimed for the rectangular map. These failures do
not justify native implementation. A full-output endpoint may instead need
a direct rank-potential argument; simply reusing square-inner inversion
would be incorrect.

Five exact geometry tests check the deterministic dimension inequality by
small exhaustive enumeration, the excluded-row expectation, both deterministic
high-rate failures, the eighteen-row-shortened probabilistic failure, and
the inconclusive next gate. The fifteen saved non-rectangular witnesses also
passed fixed-mark replay. No existing file was changed.

```text
python -B research/workstreams/char2_field_spin/theory/hd_fast_geometry_tests.py
python -B research/workstreams/char2_field_spin/theory/hd_fast_rectangular.py --iterations 14
```

The concrete next research target is the 32/17 rectangular geometry only
if a sharper joint support/state proof can close its low-middle deficit.
Otherwise the next iteration should change the causal mixing topology,
not implement a high-rate expanding inner already rejected by the window
test. The old 4K control remains the only positive proof direction considered
here, without a speed claim.

## Global-prefix strengthening

A high-rate outer cannot be rescued by changing the route alone. Let the
routed outer image be any linear subspace of F^L with dimension K. The route
may be arbitrary; no MDS or independent-row property is needed. A causal
linear inner consumes t symbols and emits n symbols per packet, carries
s symbols of state, and starts at zero. Assume t divides L and retain all
N=nL/t outputs.

Restrict the routed outer word to a prefix of W input coordinates, where
t divides W. Setting the remaining L-W coordinates to zero leaves source
dimension at least K-(L-W). Impose zero state after this prefix. These
s additional linear equations leave a nonzero source whenever

    W > L-K+s.

All later inputs and states are then zero, so every later output is zero.
If the encoder is injective, this gives a nonzero output of weight at most
nW/t and hence relative distance at most W/L. Otherwise it already fails
the required injectivity. The bound is independent of the packet's dense
or sparse input-to-output map and all sampled coefficients.

The first aligned choice is

    W = t * (floor((L-K+s)/t)+1),

provided W<=L. For fixed-width and fixed-state families, the distance ceiling
therefore approaches the outer redundancy 1-K/L. Periodic output counts
give the same argument after aligning W to a complete period.

At K20, an outer of rate2/3 with L=1572864 and (t,n,s)=(24,32,4) has
W=524304 and an output-weight cap of699072, about33.3344% of N=2K.
For padded GRS[256,192], R=5462, L=1398272, and actual source dimension
K=1048576, the (16,24,4) inner has W=349712 and cap524568. Removing the
256 excess output symbols to restore N=2K cannot increase that absolute
weight, giving a cap of about25.0134% or an injectivity failure.

Thus neither allocation can attain40% distance, even with a fully random
global permutation and a dense MDS feedthrough. This is a structural
rejection, not a failure of the current certificate algorithm. At40%, an
outer rate above roughly60% is unavailable to this causal, small-state,
proportionally expanding architecture.
