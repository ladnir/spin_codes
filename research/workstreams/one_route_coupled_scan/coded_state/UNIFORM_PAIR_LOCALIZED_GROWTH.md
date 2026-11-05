# Local input occupancy limits growth to a full state

The unchanged three-shear construction now has a restricted 46% distance
certificate for every history with at most seven defective incoming
batches, with failure probability at most 2^-77. This includes histories
with full-state resets. The probability is over the fixed setup and
quantifies over all nonzero J-valued messages, J=GF(2^32).

Higher defect counts and extension-field closure remain unproved. This
is not a whole-code certificate. The prior episode proofs and the certified
uniform-feedback control are preserved.

The main new bound is local: every input batch contains only a bounded
number of coordinates from rows that could support a compatible nonzero
outer word. A small defect budget then limits how large the state can grow.

## Fixed parameters and the support set

Use the previous model: b=512 state coordinates, L=4096 input batches,
M=8192 outer rows, and 16 batches per region. Each outer row is an
RS[256,128] word over J. Every region independently permutes its M row
coordinates. Maps and routes are sampled once for all messages.

Fix a union I of complete input batches before sampling the route.
For row i, let f_i be the number of its coordinates routed into I.
A row can carry a nonzero word supported on I only if f_i>=129.
Call such a row eligible. Its shortened dimension is (f_i-128)_+.
Thus the message-space dimension on I is sum_i(f_i-128)_+.

The [regional-rank proof](UNIFORM_PAIR_REGIONAL_RANK.md) bounds this
dimension and the number of eligible rows simultaneously over interval
sets I. We retain its positive factorial-moment bounds and exact rational
arithmetic. The refinement below also fixes a single batch D contained in I.

## Last-full-state decay has a short window

Suppose an incoming full state is followed by a zero incoming state before
Q_L. Choose the last full state Q_f before the first subsequent zero Q_z.
All states between them are nonempty and proper. The terminal state Q_L
is unconstrained and does not define a full reset in this partition.

Let e be the available number of defective batches inside live intervals.
The earlier clean-run argument bounds each proper clean run by eight
incoming states. At most e defective batches divide these runs into
at most e+1 pieces. Consequently

    z-f-1 <= 8(e+1)+e,
    z-f <= 9(e+1) =: R.

Every coordinate initially active at Q_f must die in this R-step window.
For a fixed message, route, start f, and choice of one death time for each
coordinate, the uniform-image comparison costs at most q^-b. At a selected
time, require the coordinate to be active and then zero; a currently
inactive coordinate cannot satisfy that event. Exposing updates in order
therefore retains the q^-1 cost for each selected coordinate. Two selected
coordinates in one pair jointly cost q^-2.

There are at most L choices of f and R^b death-time assignments. Including
the single global accepted-map density factor, the reset probability for
one fixed message is at most

    2 L R^b q^-b.

For c live intervals and initial-activity indicator a, the c-a later
activations consume that many defects. Hence e<=h-c+a when h is the total
defect budget. On a route event giving input dimension d(I), the bad-word
union has exponent at most

    1 + 32 d(I) + 12 + ceil(log2(R^512)) - 32*512,

before interval placement. The old full-reset bound used L^512 instead
of L R^512. It discarded the defect budget in the final decay.

## Exact input alignment avoids unnecessary dimension

The c incoming-state intervals expand backward by their activation batches.
They use T+c input batches in total, where T counts nonzero incoming states.
The final input batch is additionally available if it is not already in
an expanded interval. Thus |I|<=512(T+c+1).

If a=1, the first interval begins at batch zero. The final input is always
included and ends at the final regional boundary. There are at most c+1
intervals. Their partial-region count is therefore at most

    2c+1-a.

When the final input is already included, pad by one batch to reach mass
(T+c+1)/16 in row-coordinate units. All residual cells used here have
T+c+1<L, so a padding batch exists. Extending an existing interval into a
gap preserves the endpoint anchoring and cannot increase the interval bound.
This avoids the old padding to a multiple of eight batches.

For no-full-reset histories with T exceeding the proper-state bound, some
incoming state is full. The final live interval must then reach Q_(L-1).
There is no separate terminal batch: the exact mass is (T+c)/16, and the
partial-region count is at most 2c-a-1.

These are upper envelopes for the actual support sets, not independence
assumptions about their regional occupancies.

## Eligible rows inside one batch

Fix D contained in region r. Its size is M/16. Write mu=|I|/M and let J_0
bound the number of partial regions of I. If region r contributes a/16
of its coordinates to I, then a is in {1,...,16}. An eligible row whose
coordinate lies in D already has one allowed coordinate. It needs at
least 128 allowed coordinates in the other regions.

Let beta_rest(mu-a/16,J_0-1[a<16],128) be the retained factorial-moment
upper bound for this latter event. Define

    beta_D = (1/16) * max_(a=1,...,16) beta_rest(...),

omitting impossible negative-mean or negative-partial-count arguments.
For s specified rows, their simultaneous batch membership has probability
(512)_s/(8192)_s <= (1/16)^s. Routing in the other regions is independent
of this region. The positive moment expansion then gives

    Pr[all s rows are eligible and have their coordinate in D] <= beta_D^s.

This is a joint-row bound, not independence of row events. Consequently
the number of these rows exceeds S with probability at most
binom(M,S+1) beta_D^(S+1). The verifier chooses S with exact directed
integer bounds, as in the retained row-cap routine.

Make this event simultaneous over all I with at most m=h+2 intervals and
all L batches D. There are fewer than 2m(L+1)^(2m) choices of I. Request

    26m + ceil(log2(2m)) + 12 + 80

bits per cap, leaving an exception probability at most 2^-80. The set I
and batch D may then depend on the message and sampled setup.

Every nonzero input coordinate belongs to an active RS row. Since that
row's word is supported on I, the row is eligible. Thus every actual
input batch has weight at most S outside this simultaneous route exception.

## Deterministic state-growth bound

Each live interval starts from an input batch, so its initial state has
weight at most S. At a nondefective update, incoming pairs are empty or
fully active, and no new coordinate is born from an empty pair. The next
support is therefore contained in the current support.

At any defective update, each active coordinate can influence at most
two coordinates through its pair map. The update Q_next=AQ+X gives

    wt(Q_next) <= 2 wt(Q)+S.

After at most e internal defective batches in the interval, every incoming
state therefore has weight at most

    (2^(e+1)-1) S.

Use e=h-c+a as an upper bound. If this quantity is below 512, no incoming
full state is possible on the good route. This excludes a full-reset
history directly. It also excludes a no-full-reset profile whenever
T>9(h+a)-c, since the retained proper-state bound then forces a full state.

The e=0 special case only needs the total eligible-row cap: reaching a
full state requires a full activation batch with 512 distinct active rows.
The per-batch cap is stronger and covers the later e=1 residuals.

## Finite results

The retained full-reset verifier at h=6 had 25 residual cells. The decay
window resolves 16. Exact alignment lowers the remaining c=6,T=1968
dimension from 609 to 358, resolving that cell. The eight c=7 residuals
have e=0 and fewer than512 eligible rows, excluding full activation.
All other cells retain their old bounds. The complete h<=6 full-reset
family has good-route exponent at most -229 and failure at most 2^-78.
Together with the older no-full-reset result, this first gives 2^-68
for every history with at most six defects.

At h=7, decay and exact alignment reduce 132 old full-reset residuals to
ten. The per-batch cap excludes all ten. The hardest surviving profile
before this final exclusion is a=1,c=7,T=1968: S=156 and e=1 imply state
weight at most468. The complete full-reset good-route exponent is -263.

The no-full-reset verifier at h=7 checks61376 count cells. The local
state-growth bound excludes all379 old residuals. Its remaining old
good-route exponent is at most -310. Both branches include their route
exceptions:

| Family at h<=7 | Good-route exponent | Complete upper bound |
| --- | ---: | ---: |
| Full reset | -263 | 2^-78 |
| No full reset | -310 | 2^-78 |
| Union | — | 2^-77 |

This seven-defect result supplies an alternative proof of the already
closed six-defect cases. The more detailed episode proofs remain preserved.

## First limitation at eight defects

A bounded h=8 no-full-reset replay checks69551 cells. Of706 old residuals,
the new exclusion resolves661 and leaves45. Their families are:

| a,c | Residual T values | Internal budget e | Largest local cap S |
| --- | --- | ---: | ---: |
| 0,5 | 1950--1951 | 3 | 36 |
| 1,5 | 1940--1960 | 4 | 39 |
| 1,6 | 1943--1959 | 3 | 59 |
| 1,7 | 1955--1959 | 2 | 87 |

The largest remaining old exponent is +14000. These are failed bounds,
not exhibited low-weight words. In the first residual, S=35 and e=3
give the envelope 15*35=525, just above512. At this point the state-growth
test no longer proves that a full state is impossible.

The next refinement should charge the internal defects needed to reach
full support against the same budget used for early closed episodes.
The current separate envelopes allow all e defects to maximize state
growth and also use them to enlarge early intervals. Retain the minimum
number of growth updates implied by S, couple it to the early-episode
budget, and re-run the hardest classes before changing the construction.
If that remains insufficient, inspect the route constraints on the few
input batches responsible for growth. No redesign requirement is yet proved.

## Tests and reproducibility

`diagnose_pair_full_reset_frontier.py --test` passes54241 last-full decay
checks,2047 anchored input-geometry checks, and400 exact joint localized-row
moment checks. `--integrate --defects 6` and `--integrate --defects 7`
recompute the full old cover and refine every residual. Use `--detail`
for the compact per-residual records.

`diagnose_pair_localized_frontier.py --defects 7` checks the complete
seven-defect no-full-reset cover; `--defects 8` reproduces the residual map.
Full-reset h7 takes about30s locally, no-full-reset h7 about24s, and the
h8 diagnostic about40s. These are verifier times, not encoder timings.
All runs were serial, Python-B, with bounded caches and explicit45s caps.
Tests have a20s cap. No result files or bytecode caches were generated.

The state size, RS encoder, one regional route, pair distribution, and reset
are unchanged. The certified control and frozen earlier verifiers are intact.
