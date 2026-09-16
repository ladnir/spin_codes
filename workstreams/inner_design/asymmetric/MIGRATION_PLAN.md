# Migrating from RM2Sub to IMT

The active [finite replacement effort](../finite_migration/README.md) now
tracks every paper certificate, diagnostic slice, and SPIN timing that must
move to IMT. It has added full half-rate K22 and K24 certificates with
48.3904916829 and 46.4572549297 margin bits at 10% distance. Both have
complete occupancy coverage and 512-bit replay. K16/K18 and the manuscript
and performance migration remain open; the historical checkpoint below
records the earlier gates and K20 result.

The replacement inner is **Independent-Map Transvection (IMT)**; see the
[terminology and source-name guide](../IMT.md). The complete construction is
SPIN with a BCH outer, randomized bit-transpose permutation, and IMT inner.
The outer instance and inner parameters remain separate configuration choices.

IMT retains the structured SPIN route and recursive architecture.
It replaces the tied expansion/feedback maps with a balanced expansion A and
an independent sparse feedback map B. Each step emits X+Aq and updates the
state to Mq+BX, where M is a fresh sampled transvection in the shared setup.

The first replacement target is the measured quarter-rate instance:
BCH [128,32,32], t=128, s=19, K=2^20, feedback map `greedy3_2`.
The two targets are 16.5% relative distance with 40 margin bits and 19% with
30 margin bits. The supported default remains unchanged during this work.

## Gates for adopting the candidate

1. Produce an outward certificate for both operating points, replay at higher
   precision, and bind the exact independent maps to the tested implementation.
2. Extend the certificate producer to the other required finite instances.
   Prioritize the paper's BCH-256, t=128,s=19 cells at K=2^16,2^18,2^20, then
   its larger certified sizes. Reuse rigorous outer-spectrum envelopes there;
   an exact BCH-256 spectrum is not a prerequisite.
3. Adapt and re-certify the asymptotic inner argument. This gate has passed
   in-session review: the manuscript now uses IMT at 11% and c=39/4.
4. Implement and measure a complete forward encoder, including the forward
   BCH circuit and the inverse direction of the routing schedule. Compare
   against the original forward inner under the same workload and setup policy.
5. Switch the default and revise the paper only after the required proof and
   performance cells pass. Keep the old implementation as an artifact baseline.

These gates do not require certifying every possible t,s combination. An
unmeasured or uncertified parameter choice must not inherit a claim merely
because another instance of the family passes.

## Finite certificate work

The first gate is complete. Both operating points passed the 256-bit outward
producer, 512-bit replay, exact union and occupancy checks, and implementation
binding. `CERTIFICATE_RESULT.md` records the result and reproduction commands.

The weight-five feedback variant now has a full BCH-256 finite certificate
at K=2^20 and 10% relative distance, with 50.0620882264 margin bits. All
occupancies passed 512-bit replay and exact union checks. Its isolated
optimized transpose is now bound to the certificate, passes twelve release
tests and five sanitizer tests, and measures 10.155 ms versus 11.169 ms for
the supported build. The routing-tuned supported control measures 10.658 ms.
See `bch256/weight5/implementation/README.md`; no default has changed.

The original weight-three BCH-256 extension remains partially certified.
Its replayed sparse ranges are Q=1..91 at K=2^16, Q=1..255 at K=2^18, and
Q=1..511 at K=2^20. None is a full distance certificate yet. The first K20
dense cover has an intermediate-density bottleneck that persists at a point;
an exhaustive cross-map overlap audit improves that point by only 1.68 bits.
The zero-state contribution limits this fixed proof expression. Two
weight-five feedback candidates passed the point; the seed-0 search then
closed K=2^20. K=2^16 and K=2^18 remain unresolved.
Balanced-feedback control points also fail the present shorter-length bound.
The next engineering step is opt-in integration of the two certified
rate-specific profiles. The shorter lengths need further proof work.
The historical overlap diagnosis
is in `bch256/overlap/README.md`.

`certify.py` reconstructs the spectra of A and B^T separately, obtains the
B-kernel spectrum by exact MacWilliams inversion, and rebuilds integer fiber
caps and low-input cancellation histograms. It reuses only the generic
fixed-weight transfer from `certify_no_constant.py`.

Q=1 uses the independent-map fixed-weight matrices for input weights zero and
one. It does not use the older symmetric Q=1 shortcut. The dense Fourier
transfer bounds overlaps between supports of Aq and B^T a by their feasible
endpoints, without identifying the two maps. Each saved witness selects one
complete transfer family; entries from different representations are not mixed.

The producer uses 256-bit Arb arithmetic and upward dyadic bounds. Its
`--verify` mode recomputes all terms at 512 bits and checks enclosure by the
saved bounds. `verify_certificate.py` separately checks exact union arithmetic,
occupancy coverage, the map header, source hashes, and existing correctness
receipts. `test_asymmetric_outward.py` checks the new outward transfer against
exact rational toy-state distributions, including nonzero BA.

## Asymptotic audit

Current follow-up: the [IMT asymptotic investigation](../imt_asymptotic/README.md)
now has a complete 11% proof draft and passing numerical gates for all
occupancies, retaining c=39/4. This includes 1,023 dense boxes, an exact sparse
polynomial certificate, fixed-Q inequalities, a full-state continuum lemma,
and a summable union. The written reductions passed an in-session analytic
review and are integrated into the asymptotic section and appendix.
An IMT-only paper remains the objective, but existing finite RM2Sub statements
retain their scope until the required finite replacements close.

The pre-migration theorem was not an abstract theorem for an arbitrary
state mixer: it instantiated RM2Sub-S19 with independent nonzero field
multipliers. The following historical dependencies have now been replaced:

| Current dependency | What changes | Required replacement |
|---|---|---|
| `structured_spin.tex`, structured inner: CA=0 and field multiplication | Independent B need not satisfy BA=0; one transvection has a lazy component | General A/B/M interface; a separate instantiation statement |
| `structured_appendix.tex`, three-state envelope and positive occupancy | Exact refresh and the old association matrix are used | Independent-map weighted-state transfer and a new outward dense-domain cover |
| Same appendix, fixed and growing sparse occupancy | Continuum limits and a four-state Bernstein certificate use exact-refresh termination | New sparse limits and interval inequalities that retain lazy-state information |

For fixed nonzero q, the candidate has the law

    Mq ~ (1/2) delta_q + (1/2) Uniform(nonzero states).

The old proof instead uses exact uniform refresh. For example, a nonzero
syndrome terminates a refreshed state with probability 1/(2^s-1). With the
candidate, a lazy branch can terminate whenever BX=q. Its probability cannot
be replaced by the old uniform-refresh factor after conditioning on output
weight. The finite certificate controls precisely this distinction.

The route reduction, outer-spectrum event, and injectivity argument can be
retained at the architectural level. The new state update still yields a
block-triangular inner with identity diagonal blocks. With fixed t and s,
both directions still require O(N) arithmetic work. These observations do
not establish the old 0.11 distance threshold or its 39/4 block-growth constant
for the new maps.

Those three occupation regimes now have replacement arguments in
`../imt_asymptotic/ASYMPTOTIC_PROOF.md`, extended to 11% by
`../imt_asymptotic/d11/PROOF_UPDATE.md`. The subsequent review is recorded in
`../imt_asymptotic/d11/PAPER_REVIEW.md`, especially the uniform sparse route
comparison and componentwise full-state continuum limit.
Refining the outer-spectrum majorant recovers
the old theorem's exact 11% without changing IMT. A finite K=2^20 certificate remains distinct
from these uniform-in-length inequalities.

## Forward implementation requirements

The exact transpose uses A^T on the raw transpose input, not its emitted
value; the implementation already tests that identity with BA nonzero.
The forward implementation must compute feedback from raw X as well.
It should retain fixed-size, generated XOR circuits and preallocated buffers.
Compare complete K-to-N encoders, not a forward inner microbenchmark against
a complete transposed encoder. Include setup separately and run timings serially.

No forward performance claim is made by the current transpose measurements.
There is currently no optimized forward quarter-rate BCH circuit exposed by
the checked-in `QuarterCircuit.h`; generating and validating that circuit is
part of the forward gate, not an assumed existing baseline.
