# One unrestricted closed episode with a shared region

For the fixed three-shear model, the selected five-episode geometry now
closes with arbitrary internal activity in one closed episode. The restricted
failure bound is at most 2^-99. This joins the defect-envelope argument to
shared-region routing; it does not establish whole-code distance.

The geometry has five proper closed crossing episodes of length three,
followed by a live tail. Exactly one region is shared between intervals;
131 regions are touched in total. Four closed episodes are clean. The
remaining closed episode may have any births or one-sided pairs consistent
with its length and endpoint reset. The tail's internal activity is unrestricted.

## Cover charge intervals, not isolated charges

For an episode of length l and initial weight w, the previous identity gives
C=w+e, where e counts empty-pair births plus one-sided pair occurrences.
Since C counts zero coordinates across l live-map updates, C<=bl and
e<=bl. Thus l=3 and b=512 give a safe finite cover e=0,...,1536. Including
the impossible top endpoint only enlarges the family used in the bound.

For a band a<=e<=d, use these bounds simultaneously:

- live-map probability at most q^(-w-a);
- input-envelope size at most min(b,w+d);
- initial-support count at most A_b(w,d);
- history count at most T(l,w,d).

All support histories in the band are included by the birth-event count
for B<=d. This is a bound on their total number, not one bound per charge.
It therefore needs no additional factor for the number of integer charges
inside the band.

There is also a useful cap independent of d. Once initial support and the
zero terminal state are fixed, there are only l-1 intermediate supports.
Their total number is at most 2^(b(l-1)). Replace the history-count bound by

    min(T(l,w,d), 2^(b(l-1))).

This avoids letting a loose birth/death description exceed the count of
all possible support histories when the charge is large. It does not count
field-valued state amplitudes; those remain integrated in the forward
pair-image probability.

Substitute these quantities in the defective activation-route subset count.
The actual RS zero count, virtual shared-region deficit, finite-population
correction, and sum over random fragment allocations are unchanged from
UNIFORM_PAIR_SHARED_REGIONS.md. In particular, the actual message dimension
is three per active row, not four after the virtual region split.

## Deterministic replay

The helper diagnose_pair_defect_bands.py checks five contiguous bands and
the separate clean case. For every band, it covers all source counts
s=8,...,48, interval endpoints, and all five choices of exceptional episode.
The four clean episodes enforce the lower source count eight; this lower
bound is not assumed for the defective episode itself.

| Defect charge e | Good-route upper exponent |
| --- | ---: |
| 0 | -623 |
| 1--2 | -503 |
| 3--8 | -447 |
| 9--32 | -352 |
| 33--128 | -480 |
| 129--1536 | -3292 |

The six contributions total at most 2^-349. They use the same eligible-row
route event, so its failure probability 2^-100 is added once. The combined
restricted-family bound is at most 2^-99. No output-weight restriction is
used. The profile has arbitrary internal charge in one episode, not merely
one or two defective coordinate events.

Run `python -B diagnose_pair_defect_bands.py --replay` to reproduce the
contiguous charge cover and union. The replay took about 5.3 seconds on
Peach. Each profile has a 35-second source cap; the replay ran under a
40-second process timeout and 256 MiB address-space limit. Tests check 72
charge-band component inequalities. All work was serial and source-only.

## Coverage audit: why this is not the six-defect certificate

The global partition counts defective batches h. It records whether Q_1
is active (a), the number c of live intervals, and their total length T.
There are c-a activations from zero, leaving at most h-c+a internally
defective batches. One such batch can contain many coordinate events e.
The all-charge cover above removes that particular small-e limitation,
but does not remove its geometric restrictions.

For the retained no-full-reset h=6 frontier:

| Global branch | Internal defective-batch budget | Geometric gap in the new family |
| --- | ---: | --- |
| a=0,c=6 | 0 | Five early clean episodes may stay within one region or share more than one region. |
| a=1,c=7 | 0 | Six early episodes; the first activation is aligned at batch zero and does not cross a boundary. |
| a=0,c=5 | 1 | Four early closed episodes, not five; lengths need not equal three. |
| a=1,c=6 | 1 | Five early episodes, but the first is aligned; lengths and shared regions vary. |

Other h=6 branches, the full-reset branch, and defect-rich histories also
remain in the partition. The earlier complete h<=5 certificate is preserved.
No updated global residual count is claimed from the new restricted replay.

The next high-value integration step is a single-region/aligned activation
factor, followed by an explicit test of the retained a=1,c=7,T=1938 residual.
That branch has no internal defects, so more charge refinement cannot by
itself settle it. A coverage pass must include all permitted early lengths,
noncrossing intervals, sharing patterns, and the correct endpoint union.

The new probability/envelope identity remains useful for the one-internal-
batch branches once their geometries are represented. Keep the outer,
route, state size, and maps fixed. No structural counterexample or necessary
extra mixing operation has been established. Whole-code and extension-field
closure remain open; the certified control is unchanged.
