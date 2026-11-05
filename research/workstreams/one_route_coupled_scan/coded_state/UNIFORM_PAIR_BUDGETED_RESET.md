# Charge full-state intervals and their growth to the same defect budget

This continuation keeps the RS outer, regional route, 512-coordinate state,
and rejected-uniform SL2 pair maps unchanged. It refines the full-reset
branch of the [eight-defect analysis](UNIFORM_PAIR_GROWTH_BUDGET.md).
All probability bounds concern setup sampled once for all J-valued
messages, where J=GF(2^32). They are not extension-field or whole-code claims.

The earlier full-reset count discarded two useful facts. Full support must
first be created, and the updates that create it already supply random
inverse-input coordinates. Both facts use the existing defect budget.

## Input caps constrain the state intervals

Fix the incoming-state interval endpoints before sampling the route.
As before, h bounds defective batches, a indicates Q_1 nonzero, c counts
nonzero incoming-state intervals, and T counts their total length.
Activations after Q_1 consume c-a defects, leaving e=h-c+a internal defects.

On the localized good-route event, every input batch has weight at most S.
Suppose 0<S<512. Every interval then starts proper: its initial state is
its activation input. A clean update cannot increase support. Consequently,
each entry into full support consumes an internal defect, and there are
at most e full-state runs. This replaces the earlier h+1-run envelope.

Let l=floor(log2 S). Before its first internal defect, an interval has
at most l clean incoming states. A longer clean run would require at least
2^(l+1) coordinates at activation, by backward closure under the cycling
matchings. After a defect, a clean run has at most eight proper incoming
states. It may contain an arbitrarily long full prefix, but its proper
suffix still has length at most eight. Including the defective incoming
states gives

    P <= c*l + 9e.

Thus at least F=max(0,T-c*l-9e) incoming states are full. Select the earliest
F full states; this preserves the bound of e runs. Their coordinates touch
at most floor((F+30e)/16) regions. If this is at most128, all their RS parity
columns are independent, and their parity-check rank is exactly512F.
Otherwise the existing simultaneous regional-rank bound applies with
e runs, not h+1 runs. The full-state subset is chosen after fixing an output;
the route event is simultaneous over all admissible interval sets.

The cap S=0 excludes a full state. For S>=512 none of these stronger
activation-based conclusions is used.

## Keep the coordinates that create the first full state

Consider the first full state in the history. Every preceding incoming
state is proper or zero. To reach state weight at least R>S from an
activation of weight at most S, some preceding internal update is defective.
Choose the last such update. Subsequent clean steps do not increase support.

Let f be the number of coordinates in its live pairs and w its incoming
state weight. The actual update Q_next=AQ+X gives

    R <= f+S,    f <= 2w.

The integer f is even. Hence this update supplies at least
2*ceil((R-S)/2) random inverse-input coordinates, and its incoming state
has weight at least ceil((R-S)/2). Repeat this argument backwards until
the required incoming weight is at most S. The selected updates are
distinct and all precede the first full state.

Starting at R_0=512, define

    R_(j+1)=ceil((R_j-S)/2),    b_j=2R_(j+1),

while R_j>S. Let g count these steps. This is also the minimum integer
with (2^(g+1)-1)S>=512. If g>e, the full-reset history is impossible.
Otherwise the selected proper batches supply at least sum_j b_j additional
random coordinates, disjoint from the F selected full batches.

Each precursor batch adds at most one region. When

    floor((F+30e)/16)+g <= 128,

the combined columns remain independent in each RS row. The inverse count
may therefore use rank at least512F+sum_j b_j. This is a deterministic
rank argument on the good-route event. It neither assumes new independent
state bits nor multiplies independent probabilities for growth events.

For example, the remaining h=8,a=1,c=7,T=1967 profile has S=133 and e=2.
Backward requirements are512 ->190 ->29, giving380+58=438 additional
random coordinates. Its full-batch rank is972800; the strengthened rank
is973238. The combined region envelope is124, below the128-column threshold.

## Growth also limits the reset window

Choose a full-reset interval. Reaching its first full state uses at least
g internal defects. At most e-g remain after its last full state before
the subsequent zero. The retained decay argument therefore gives window
R=9(e-g+1), rather than9(e+1). Its fixed-message bound remains
2 L R^512 q^-512. The message union uses the supported outer-space dimension
and the same interval endpoints. The verifier takes the better of this
reset bound and the inverse rank bound; it does not multiply the two.

## Route exceptions and finite checks

For fixed a,c there are fewer than2^(13(2c-a)) endpoint tuples: the first
start is fixed when a=1. The free final input batch adds no variable
endpoint. Input support has mass at most(T+c+1)/16 and at most2c+1-a partial
regions. Use the previous deterministic padding convention if needed.

Let E=13(2c-a)+ceil(log2(2(h+1))). Request E+12+1+80 bits in each localized
batch cap and E+7+1+80 bits in each dimension-layer cap. The factors cover
4096 batches and128 dimension layers respectively. The two event types
together have exception at most2^-80. The old route event and the newly
used full-run rank event each have their own2^-80 allowance. Thus, if the
summed good-route contribution is at most2^-80, the complete full-reset
branch has failure at most4*2^-80=2^-78.

`diagnose_pair_budgeted_reset.py --test` checks31239 support transitions
for the state-count bounds and526 first-full precursor cases on four
coordinates. These exhaustive small cases test the formulas; the arguments
above establish the general bounds.

The first full replay, without precursor rank or exact-volume refinement,
left five of341 old residual profiles. Four of those close under the
existing output-volume bound using an exact binomial coefficient alone. The remaining profile
is the S=133 example above: its exponent changes from+12851 to-1165 after
adding the438 precursor coordinates. The coarse reset bound there is+2484.

The verifier uses exact binomial volume only where the retained coarse
bounds fail. `--coarse --no-precursors` preserves the initial ablation.
`--no-precursors` isolates the volume refinement. No frozen helper is edited.

The complete replay covers all341 old full-reset residuals, with none
remaining. Of these,72 retain the old decay-window bound; four need the
new exact-binomial calculation. The maximum good-route exponent is-108,
at a=0,c=5,T=1963, handled by the old decay window. Adding the three route
exceptions gives full-reset failure at most2^-78. Together with the
preserved no-full-reset result, every history with at most eight defective
incoming batches has restricted failure at most2^-77 at46% distance.

The successful serial replay took46.15 seconds under a60-second cap.
An earlier attempt hit that cap because it recomputed exact coefficients
even for cases already closed by the old decay bound. That incomplete run
was not treated as a certificate. Skipping those unnecessary computations
preserves all bounds and gave the complete replay above. Tests and numerical
checks used Python-B and wrote no result files. The frozen defect helper's
SHA256 remains7238fec9e5553e7b4d1ce3edb3089cec06abcc34cb090c029d0f3d6b54bbaed3.

## Scope and next discriminating step

The improvement is a relation between routing, growth, and inverse rank,
not a change to the encoder. It does not establish that defect-rich
histories are covered. Before extending the count one budget at a time,
test whether the precursor rank and full-run count continue to control
substantially larger defect budgets. If not, retain the actual pair support
lost by the batch-level defect summary, and compare that loss with the
cost of a small inner modification. A failed envelope is not itself a
counterexample to the construction.
