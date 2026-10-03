# A positive all-occupancy width-eight proposal

The active objective remains a complete 10% distance / 40-bit setup-failure
bound at K=65,536, with precomputed transposed encoding near or below 100 us.
The [preceding iteration](../iteration3/README.md) improved the proof but did
not close it. Its strongest tested actual24 q=16 proposal is 26.75 bits.
The measured 91.12 us kernel uses the original 16-bit maps, not actual24.

This iteration obtains a **68.893994681-bit all-occupancy floating bound**
for a wider outer and the actual 24-bit inner. This is the first positive
complete occupancy screen in this width-eight study. It is not yet an
outward-rounded certificate, and the changed encoder has not been timed.
The active goal therefore remains incomplete.

## Construction selected for the next implementation gate

At K=65,536 and N=131,072, use 256 outer groups. Each group maps 256 bits
to 512 bits through eight parallel GF16 RS[16,8] rows. Identify each aligned
eight-nibble column with a 32-bit symbol, and apply an independent uniform
nonzero-transitive 32-bit linear randomizer to each symbol. Split these into
64 bytes and shuffle them uniformly within the group. The regional route
then independently shuffles the 256 group slots in each of the 64 regions.
Two independent 16-bit symbol randomizers do not realize this distribution.

The physical inner still consumes eight bytes per step. Its 24-bit state
is (a,b,c) over the AES polynomial-basis GF256, with

```
A(a,b,c)_h = a + h*b + h^2*c,       h=0,...,7,
C(X) = (sum X_h, sum h*X_h, sum h^2*X_h),
Y = X + A(state),                  next = M*state + C(X).
```

There are 32 physical steps per region and 2,048 in total. State starts at
zero, continues across all regional boundaries, and has no final flush.
The literal local moments were computed with independent uniform GL3 updates.
For implementation, fresh uniform nonzero GF(2^24) scalar updates, including
their binary adjoints, preserve the fixed-message distribution used by the
proof. They do not produce the same ensemble of full encoder matrices.
The [cost note](COST_OPTIONS.md) gives the six-product scalar circuit.

## Evidence and numerical scope

The proof conditions on potential occupancies before clipping the conditional
first moment. It groups four chronological physical operators before taking
entrywise fractional powers. The group-subset union factor remains outside
those powers. Regional and global products preserve continuous state.
For q=1, the calculation uses exact expected outer support counts separately.

`wider_actual24_allq_v2.json` checks every q from 1 through 256; it does not
interpolate between the initially sampled points. It sums their probability
upper-bound proposals at cutoff 13,107, corresponding to the 10% target.

| Quantity | Floating margin |
|---|---:|
| Combined bound over all nonzero occupancies | 68.893994681 bits |
| Weakest individual occupancy, q=3 | 68.894675 bits |
| q=1 from exact expected outer support counts | 110.785884910 bits |

The initial, coarser tilt collection gave 55.538329393 combined bits.
Adding two cached tilts improved coverage between the sampled points.
[Independent review](AUDIT_WIDER.md) checks the construction and conditional
bound and replays the numerical witnesses using a different global recurrence.
These checks do not replace outward arithmetic for a final certificate.

The wider outer with the original 16-bit state does not close the sampled
small-support cases: q=8 gives -19.08 bits and q=32 gives -27.11 bits.
Expansion scaling improves these only to -18.09 and -6.63 bits. Negative
values mean vacuous proof bounds, not observed low-distance words.

## Directions tested and stopped

- [A-only byte scaling](SCALED24.md) gains at most 0.375 bits at actual24
  q=16 in the four-step gate. Fresh literal censuses preserve the feedback
  map and reproduce the unscaled baseline exactly. Stop this search.
- [Local filling slack](LOCAL_SLACK.md) can recover less than 1.5 bits at
  that q=16 fractional proof setting. This is insufficient for the missing
  13.25 per-class bits, before reserving the full union budget. Deleting all
  zero returns would remove genuine transitions and is not a valid fix.
- The old t32/s16 design was rechecked under fractional routing. It gives
  -126.06 bits at q=16 with four-step grouping and -109.11 with eight-step
  grouping. The q=119 four-step bound improves to +85.89, but the smaller
  occupancies still fail. Its lower arithmetic cost does not justify
  implementation before a new proof idea.

## Performance gate and remaining work

The positive proposal's wider outer and scalar24 inner add about 30.4% to
the baseline whole-call GFNI count, plus XOR work. Routing volume and RS
parity work do not grow. This does not predict a proportional slowdown or
establish that the code can meet 100 us.

Next, implement this candidate in an isolated experiment. Preserve the
measured byte-packed layout and shared-parity circuit; importing an older
wrapper would add avoidable repacking and parity work. Check literal
forward/transpose adjoints and scalar/SIMD equality before serialized A/B
timing against the frozen 91.12 us byte kernel and certified nibble control.
If the cost is competitive, finish outward evaluation and independent
certificate review. Otherwise retain the positive proof as a checkpoint
while revisiting the construction. Do not silently trade away the time target.

## Replay and preservation

`wider_fractional.py` covers the two 16-bit-map wider controls.
`cached_wider24.py` combines authenticated 24-bit local caches with the
wider outer. `wider_coverage.py` replays all occupancies from those witnesses.
`narrow_fractional.py` records the t32 controls. All scripts require fresh
output paths and preserve earlier sources and receipts.

End-of-iteration verification: 122 portable tests pass across the width-eight
study, including 25 new tests here. The third and fourth iterations contain
80 completed receipts with 1,403 matching source pins and no mismatches.
The certified K16 checkpoint still matches all 220 source pins and its archived
receipt. No encoder benchmark, production edit, commit, or promotion occurred.

Each changed map requires its own literal map identity and local moments.
Diagnostic optimistic kernels are not constructions or probability upper
bounds. Floating proposals are not certificates. Preserve all earlier source
pins and certified checkpoints; write new work here. No production change,
promotion, commit, or parallel encoder benchmarks are authorized by this work.
