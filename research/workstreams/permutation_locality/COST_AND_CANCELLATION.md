# Sequential cost diagnostic and grouped cancellation

2026-09-23. The 2x goal remains open. Production code and the paper are unchanged.

## Cost diagnostic

The column-major candidate's slow inner phase is not explained by the tested
aliasing annotation. Marking its separately allocated scratch pointer `restrict`
gives 12.929 ms versus 12.904 ms without it. These are medians of three
101-call process medians, alternating order. Phase means remain about 7.9 ms
for inner/routing and 5.1 ms for transpose/BCH. Input and output may alias;
neither receives the annotation.

To remove routing from the experiment, `sequential-control` sends each inner
output directly to the next packed BCH input position. It uses the existing
unrolled IMT and four-row BCH kernels. The reference check uses the matching
materialized permutation. This permutation is only a cost diagnostic: complete
outer rows are confined to a short interval, so it is not a distance candidate.

Matched measurements use K=2^20, 128-bit elements, CPU 15 on Peach, the same
compiler options and buffers, three warmups, and three 101-call processes.
All benchmarks acquire the three shared locks and execute serially. Setup,
allocation, initialization, and correctness checks are excluded.

| Full encoder | Process medians (ms) | Median (ms) |
|---|---|---:|
| Production distribution and kernel | 9.834, 9.677, 9.651 | 9.677 |
| Sequential cost diagnostic | 6.094, 6.069, 6.100 | 6.094 |

A separate instrumented process measures 6.055 ms overall, with phase means of
2.398 ms for the inner and 3.666 ms for BCH. Phase means include warmups;
they explain this diagnostic rather than predict other pipelines.

The 2x target is about 4.84 ms against this control. Even this sequential
implementation would need roughly another 21% reduction to reach it.
This is not a lower bound on all implementations. It shows why cheaper routing
alone, followed by the current separate BCH stage, has not reached the target.
Fusion, arithmetic changes, or a different layout may change these costs.

## Exact cancellation boundary for four-row groups

Use the c=1 distribution in ANALYSIS.md and the actual 19-by-128 feedback
matrix B. An epoch contains 32 aligned four-position windows. A nonempty
support in a window is one of its 15 nonempty lane subsets.

The new enumeration checks all 480 window/subset choices. Two distinct windows
never have equal nonzero feedback syndromes, including choices with different
Hamming weights. Together with the existing one-window check, this proves:

> Every nonzero epoch input confined to at most two aligned four-position
> windows has nonzero feedback B*x.

This is an exact finite computation on the selected matrix, not a sampling
claim. The XOR of the two windows' syndromes is zero exactly when they agree.
The computation checks every such pair.

Consequently, a nonzero message confined to at most two fixed four-row groups
cannot keep the forward inner state identically zero, for any route in this
c=1 family. At the first active epoch, the input occupies at most two windows.
Its feedback is nonzero, so the next state is nonzero. A nonzero BCH row
occupies at least 38 regions, ensuring that this first epoch is not the final
one. This conclusion applies to independent or shared row-coordinate shuffles.

Three windows already admit cancellation. Exact enumeration finds 27 unordered
triples of nonempty supports in distinct windows whose syndromes XOR to zero:

| Total input weight | Number of triples |
|---:|---:|
| 4 | 11 |
| 6 | 10 |
| 8 | 6 |

These are local patterns, not BCH messages realizing a zero-state trajectory.
They prevent extending the two-window statement to three windows. For aligned
eight-position windows, even two windows admit cancellation: there are 19
unordered pairs of nonempty supports with equal syndromes.

## A conditional bound for proof development

The four-row census also supplies an atom bound for an epoch with q active
groups, where 1 <= q <= 32. Condition on their assignment to this epoch and
on each group's active-lane count. Leave their within-epoch window positions
and lane permutations random. Fix one selected group with active count a>0.

Reveal the other q-1 active groups' windows and lane subsets. The selected
group has m=33-q remaining windows, and each of its binomial(4,a) lane subsets
is equally likely. The census proves that these m*binomial(4,a) choices have
distinct feedback syndromes. Hence, for any target fixed independently of
the selected group's remaining randomness,

    Pr[B*x_selected = target] <= 1 / ((33-q) * binomial(4,a)).

The target may incorporate the revealed groups' feedback. It may also include
an incoming state contribution, provided that contribution is independent of
the unrevealed within-epoch choices under the stated conditioning. This is
not a bound under arbitrary additional conditioning on the route or trajectory.

This local bound does not count messages, control output weights, or certify
all epochs jointly. Those steps still require the grouped occupancy process
and suitable bounds on joint BCH supports. No existing 40-bit certificate is
transferred to this distribution.

## Replay and next step

`run_alias.sh` compares the archived pre-annotation binary
`build/locality-before-restrict` with `build/locality`. `run_floor.sh` runs the
matched sequential diagnostic. Raw receipts remain outside version control
under `/tmp/spin-locality-JjY7yR/measurements/{alias,floor}` on Peach.

`locality 14 4 1 pair-census 17` reproduces the four-row pair/triple census;
replace 4 by 8 for the eight-row pair census. The seed and message length do
not affect this local matrix enumeration. `verify-column` passes the dense
transpose, adjoint, and suffix checks; all four existing library tests pass.
Full sanitizers for the experimental kernels remain outstanding.

After the final diagnostic-only change, executable SHA256 is
`fc27c5b2bd78e227e2991304bb6106516554f3fa635ec9e1b4b7158a45e78bdf`.
The production library remains unchanged with SHA256
`ecbebbfa16d15b4b5460c8811053720d1dd534f2aa81ea5b82757bebd153f0d7`.

Next, test whether fusing a cache-local route with BCH can beat the separate
stages in the sequential diagnostic before investing in another broad sweep.
Keep four-row groups as the preferred proof starting point. They have a useful
two-group cancellation exclusion; larger locality units already encounter
the obstructions recorded in SECOND_ITERATION.md.
