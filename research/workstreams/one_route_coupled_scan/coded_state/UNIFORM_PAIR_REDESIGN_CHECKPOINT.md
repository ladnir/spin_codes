# Redesign checkpoint: proper-state runs, not full-state randomness

The fixed three-shear construction retains its restricted46% certificate
for histories with at most eight defective incoming batches, at failure
at most2^-77 over setup. The [eight-defect proof](UNIFORM_PAIR_BUDGETED_RESET.md)
is preserved. No whole-code or extension-field certificate is claimed.

This pass tests whether that argument scales beyond the immediate frontier.
It identifies a specific reason to discuss an incremental inner redesign:
the current local update permits long proper-state runs, and each restart
or internal defect can introduce another run. The count loses too many
random coordinates when these runs accumulate. Stronger local mixing can
shorten this delay. Nothing here exhibits a bad codeword in the current code.

The construction and implementation remain unchanged. The comparison maps
below are mathematical controls, not promoted kernels or measured designs.

## What the current proof loses

The512 state coordinates are indexed by nine bits. A pair update mixes
along one bit, and nine successive phases cover all nine bits. A clean
proper run can last eight incoming states. This length is attainable as
support geometry: the256 coordinates with top bit zero are closed under
the first eight bit matchings, but not under the ninth. This does not mean
the outer code can realize the corresponding input sequence.

For h=12,a=1,c=12,T=1967, the localized route cap is S=424 and the internal
defect budget is e=h-c+a=1. The current bound permits105 proper states,
leaving1862 guaranteed full states. Counting at least two random coordinates
in every proper state gives only953554 guaranteed random coordinates.
The allowed output weight is W=floor(.46*2^21)=964689, which is larger.
The fixed-message Hamming-ball bound on those random coordinates is then
vacuous, even before the message union.

This is not repaired merely by counting the minimum support within each
clean run. A clean proper run of length l<=8 supplies at least

    2^l + 2^(l-1) + ... + 2 = 2^(l+1)-2

random coordinates, by backward closure. Its loss relative to l full
batches is at most3586. An internal defective proper batch loses at most510
coordinates and can start another clean run. Therefore a stronger total
lower bound is

    D >= 512T - 3586c - 4096e.

For this profile that is959976, still below W. The calculation is a limit
of this aggregate summary, not an upper bound on actual randomness and
not an impossibility theorem for another proof. A further proof would need
to exclude the low-mass histories jointly with the outer constraints, or
retain considerably more state geometry. Sharper parity rank alone cannot
make a forward volume bound useful when its retained random mass is below W.

The other available bound is also insufficient here: supported outer
dimension3619 makes the full-reset exponent+101404. In contrast, a diagnostic
that supplies512 random coordinates at every nonzero incoming state gives
an inverse exponent below-992000 on the same interval parameters. The
outer and route are therefore not, by themselves, preventing this profile
from closing in a stronger-inner model.

## Compare stronger local maps without changing the outer or route

For a control of block size b in {4,8}, sample independent uniform b-by-b
matrices A over J at setup, rejecting unless A and A+delta I are
superregular. Keep the same recurrence

    Y=X+(A+delta I)Q,    Q_next=delta Q+Y.

Each accepted map is fixed for all messages. These are not the current
SL2 maps. Before rejection, Au is uniform on J^b for each fixed nonzero u.
A union bound over square minors rejects with probability at most

    epsilon_b = (2/q) sum_(k=1)^b k*binom(b,k)^2.

Conditioning raises each point probability by at most1/(1-epsilon_b).
Across N/b maps, the verifier bounds the logarithm of total inflation by
ceil((3/2)(N/b)epsilon_b/(1-epsilon_b)). This is one bit for b=4 and ten
bits for b=8. The controls thus have an explicit point bound; they are not
just assumptions that state coordinates become independent.

For b=4, cycle through index-bit pairs (0,1),(2,3),(4,5),(6,7),(8,0).
Any five successive phases cover all nine bits, so a clean proper run
has length at most four. For b=8, use triples (0,1,2),(3,4,5),(6,7,8),
giving at most two proper clean states. The generic proper-state envelope
becomes P<=d(h+a)-c, for d=5 and3 respectively. An occupied block supplies
b random inverse-input coordinates. The previous regional rank, input
dimension, growth, and decay calculations then apply to these controls.

The notion of a defective batch uses the control's own block partition.
Equal h values compare analogous count profiles, not identical sets of
messages across constructions. A fixed-code certificate cannot be
transferred by just relabeling its defect count.

## Selected numerical comparisons

The table gives good-route log2 upper bounds for full-reset count profiles,
including their endpoint-placement allowance. Positive entries mean the
bound fails; negative entries are not whole-code certificates. Each test
uses a=1,c=h,T=1967 and the same RS outer and regional route.

| Defect budget h | Current pair model | Uniform4 control | Uniform8 control |
| --- | ---: | ---: | ---: |
| 12 | +101404 | -264685 | -631118 |
| 16 | +244421 | +85505 | -435116 |

The current column takes the best retained pair, inverse, forward, and
reset bounds. In particular, a cap S>=512 must not cause the comparison
to discard the available reset bound and exaggerate the pair-model gap.

A second test checks T in {1936,1952,1967,1984,2000} at both budgets.
The pair model passes one of these ten selected profiles; uniform4 passes
six; uniform8 passes ten. This is deliberately a small diagnostic grid,
not an exhaustive cover or a new finite certificate.

There is an important nonmonotonic case. At h=12,T=1936 the pair model
excludes full support because S=169 and e=1 imply3S=507<512. The uniform4
model can grow faster and does not inherit that exclusion; its current
bound is+17097. Stronger mixing changes both propagation and which
histories become possible. Each candidate needs its own full proof.

## A cheaper next candidate and its known obstruction

The smallest natural implementation change is a second crossing pair
layer inside an update, without another global permutation. On four
coordinates, first mix pairs(0,1),(2,3), then(0,2),(1,3). Reuse the accepted
three-shear maps independently at each layer.

This product is not a uniform four-coordinate map. For input(1,1,0,0),
the first active pair can collapse onto one coordinate. The crossing
layer then produces only two nonzero coordinates. Exact small-field
enumeration gives probabilities2/7 over GF8 and2/15 over GF16 for that
event. By contrast, input(1,0,0,0) always reaches all four coordinates:
superregularity forbids cancellation in a one-active-coordinate pair.

Thus the next proof should separate ordinary four-coordinate spreading
from an internal cancellation event. A first-layer zero in a live pair
costs roughly one field constraint; it must be retained, not hidden by
claiming a uniform four-coordinate output law. The existing local point
bounds can count such events, but their aggregate cost through the inner
and outer remains unproved. This is the concrete risk to examine first.

Arithmetic also matters. Ignoring common input/output work, the present
three-shear pair uses1.5 field multiplications per state coordinate. Two
pair layers use3. Dense uniform4 and uniform8 maps use4 and8 respectively.
These are operation counts, not encoder timings. In particular, uniform8
is a proof control, not a recommendation to pay its dense arithmetic cost.

## Decision and handoff

Stop this fixed-construction exploration at the requested redesign
checkpoint. Preserve the eight-defect result and the uniform-feedback
certified control. The current evidence supports testing an additional
local crossing layer, with an explicit cancellation count, before spending
another sequence of turns expanding the same batch-level frontier.

First establish the actual four-coordinate transition bound and its cost.
Use uniform4 and uniform8 as diagnostic controls, not substitutions for
the candidate's distribution. If the exceptional cancellation count loses
the gain, test a less frequent stronger local refresh rather than silently
changing the outer or adding another global route. No such redesign is
implemented by this checkpoint.

## Reproduction and checks

`diagnose_pair_scaling.py` runs six separated (h,c) profiles in at most45s.
`--window` runs the ten selected length profiles, also capped at45s; the
final replay took31.14s. `--test` checks8704 backward-closure cases,
20212 toy transitions for the clean-run mass inequality, and the exact
accepted-pair cancellation counts above. All tests pass. The tests use
Python-B, run serially, and create no numerical files. They do not replace
the general arguments or provide missing profile coverage.
