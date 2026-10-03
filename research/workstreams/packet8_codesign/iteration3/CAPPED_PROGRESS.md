# Third capped-occupancy statistic

The first selected point is positive for the **actual 24-bit-state** width-eight
construction. At outer occupancy q=119 and weight tilt theta=0.5, the H3 route
condition gives a floating good-message margin of 168.231009 bits. The exact
rational bad-route term limits the combined bound to approximately 61.990946
bits. This is not a whole-code certificate: the good-message calculation is
not outward-rounded, and other occupancies have not been checked in this run.

## Unchanged construction and new condition

The calculation uses the authenticated local operators in
`../larger_state/trajectory_v1.json`. They describe the GF(256)^3 state, not the
16-bit-state implementation. The outer has 512 groups and 32 byte regions at
K=65536, N=131072. Each region has 64 physical steps with eight packet slots.
The state starts at zero, continues across regions, and has no final flush.

Fix q group identities. Let j be their potential-slot occupancy in a physical
step, before any zero byte labels are removed, and define

```
H3 = sum over physical steps of min(j, 3).
```

The route condition requires H3>=3183 for every 119-group subset. With rational
Chernoff variable x=77/324, the exact integer calculation in the frozen
`../route_conditioning/cap_counts.py` proves that the probability of any
violating subset is at most 2^-60. Its rational bound has estimated margin
61.990946 bits and certified integer margin 61 bits. The union includes
binomial(512,119) group subsets, not their message labels. At this same rational
variable, threshold 3184 does not pass the 60-bit check.

For the good-message term, the driver first averages the local active-packet
operators over Binomial(j,255/256). It then multiplies each resulting operator
by exp(nu*min(j,3)), performs the ordered fixed-q regional placement, and applies
the scalar factor exp(-nu*3183). The message factor remains
binomial(512,q)*beta^q, where beta=2^256/65535^8. There is no second binomial
mixture over the potential occupancy q. Splitting the bad-route event from
the good-message term does not require independence between the good event
and the message weight.

## Selected results

All entries below use q=119 and theta=0.5. The good-message column is floating.
The combined column adds the exact rational route bound, with its logarithm
evaluated in floating point for this proposal.

| Route tilt nu | Good-message margin (bits) | Combined margin (bits) |
|---|---:|---:|
| 0 | -3810.646333 | -3810.646333 |
| 1 | -2068.740125 | -2068.740125 |
| 2 | -739.601013 | -739.601013 |
| 3 | 28.564871 | 28.564871 |
| 4 | 168.231009 | 61.990946 |
| 5 | -209.248324 | -209.248324 |

For comparison, the previously saved H2-only gate at this same construction
and weight tilt gave -258.391144 bits at nu=3.5. This comparison concerns proof
bounds, not measured failure probabilities or measured encoder performance.

## Reproduction and scope

Run from the repository root, with a fresh output path:

```
python -m unittest discover -s research/workstreams/packet8_codesign/iteration3 -p test_capped_gate.py
python research/workstreams/packet8_codesign/iteration3/capped_gate.py --nus 0 1 2 3 4 5 --output NEW_RECEIPT.json
```

`cap3_saved24_v1.json` records the completed run, the rational route witness,
all local-source hashes, the saved operator receipt hash, and a successful
final source-pin verification. No new 24-bit-state census was performed.
The three tests check ordered support placement with joint cap markers,
reduction when a second marker is disabled, and invalid-marker rejection.

## Expanded gates

Independent review confirmed the potential-slot conditioning, the strict
bad-route cutoff, and the exact rational witness. A separate batched replay
of the q=119, nu=4 point differed from the original log moment by only
2.73e-12. This agreement checks the implementation, not outward rounding.

The expanded H3 screen reuses eight authenticated local caches with weight
tilts 0.1 through 0.8. `cached_screen.py` batches candidate evaluations in
scaled floating arithmetic, then replays every selected occupancy winner
with the frozen all-logarithmic evaluator. The search backend permits
subnormal underflow; the reported winner must agree with its log replay to
within 1e-7 in log moment. Every listed winner passed that check, with the
largest discrepancy below 4e-11, including the later extreme-tilt checks.

| q | Selected theta | Selected nu | Good-message margin (bits) |
|---:|---:|---:|---:|
| 16 | 0.1 | 1 | -1070.574900 |
| 32 | 0.1 | 1 | -1500.370733 |
| 48 | 0.2 | 2 | -1951.632053 |
| 64 | 0.3 | 3 | -1824.713364 |
| 80 | 0.4 | 3.5 | -1248.122683 |
| 96 | 0.5 | 3.5 | -601.690015 |
| 119 | 0.5 | 4 | 168.231009 |
| 128 | 0.5 | 3.5 | 339.936749 |
| 144 | 0.6 | 4 | 1108.343787 |
| 160 | 0.6 | 4 | 1487.649592 |
| 192 | 0.7 | 3.5 | 2114.805563 |
| 224 | 0.8 | 3.5 | 2466.681157 |
| 256 | 0.8 | 3.5 | 2371.218449 |
| 384 | 0.8 | 3.5 | -4930.588309 |
| 512 | 0.8 | 1 | -19765.608039 |

These are selected grid values, not optimized bounds for every q. In
particular, the smallest and largest occupancies select boundary weight
tilts. Smaller tilts for q=16,32 and larger tilts for q=384,512 are being
checked separately. At q=512 every physical step has eight potential
packets, so H3=6144 deterministically and the bad-route term is zero.

At q=64, a separate comparison used theta in {0.2,0.3,0.4} and each route
tilt in {0,1,2,3,4}. H2 gave the best selected bound:

| Route condition | Thresholds | Selected theta and route tilts | Good-message margin (bits) |
|---|---|---|---:|
| H1 | 1057 | 0.2; nu1=3 | -1325.523113 |
| H2 | 1615 | 0.3; nu2=3 | -951.525638 |
| H3 | 1849 | 0.3; nu3=3 | -1824.713364 |
| H1 and H3 | 1057,1849 | 0.3; nu1=2,nu3=2 | -1294.182404 |

The joint calculation sums the two exact rational bad-route bounds. It
does not assume the route conditions are independent. Every individual
condition has an exact 60-bit witness. None of these q=64 bounds closes.

The unchanged 16-bit-state control also remains negative at q=119. With
theta in {0.4,0.5} and nu in {3,4,5}, its best H3 good-message margin was
-602.299125 bits at theta=0.5,nu=3. The positive 24-bit-state points therefore
must not be reported as a certificate for the existing 16-bit-state kernel.

Receipts for these additions are `cap3_saved24_q64_v1.json`,
`cap3_saved24_q128_v1.json`, `cap3_saved24_q192_v1.json`,
`cap3_adapted_v1.json`, `cap3_coarse_v1.json`, `cap3_baseline16_v1.json`,
and `q64_cap2_v1.json` / `q64_caps1_v1.json` / `q64_caps3_v1.json` /
`q64_caps13_v1.json`. Every receipt verifies its source pins at completion.

## Extreme-tilt checks and stopping point

The second cache contains small weight tilts 0.01,0.02,0.04,0.06,0.08 and
large tilts 1,1.4,1.8,2.2. Each point below used route tilts {0,1,2,3,4}.

| q | Selected theta | Selected nu | Good-message margin (bits) | Combined margin (bits) |
|---:|---:|---:|---:|---:|
| 16 | 0.04 | 0 | -665.080292 | -665.080292 |
| 32 | 0.08 | 0 | -1506.200644 | -1506.200644 |
| 384 | 1.8 | 1 | 4670.919010 | 63.111300 |
| 512 | 2.2 | 0 | 4064.288499 | 4064.288499 |

The earlier theta=0.1 result remains slightly better for q=32. These checks
show that the negative high-q points were weight-tilt grid effects. The
low-q gap remains. At q=16 the new best choice uses a zero route marker, so
H3 itself does not improve that point. Receipts are
`cap3_small_tilts_v1.json` and `cap3_large_tilts_v1.json`; both finish with
valid source pins and independent all-log winner replays.

The root's separate `q64_h2_diagnostic_v1.json` studies the current best
q=64 comparison, at theta=0.3 and nu2=3. In its normalized positive comparison
path sum, about 1540.67 of 2048 steps are empty zero-to-zero transitions.
There are about 19.23 births and 18.23 returns; occupied zero-feedback
zero-to-zero transitions contribute only about 0.00258 steps. The comparison
has about 1561.26 active byte packets out of 2048 potential packets. These
are diagnostics of the upper-bound path sum, not the actual encoder's
failure distribution. They point toward the remaining treatment of zero
labels and long zero-state stretches, not another small adjustment of the
H3 tilt.

No further cap grids were run after these extreme checks. The recommended
next proof work is to understand or tighten the low-occupancy comparison,
using the q=64 H2 diagnostic as a concrete target, before attempting a full
occupancy certificate. A positive middle or high occupancy does not establish
the requested 10% distance with 40-bit margin for the complete construction.
No full occupancy sum or outward certificate has been produced.

## Separate fractional follow-up

The root subsequently requested a different bound, without a capped-route
condition. For each ordered tuple of four potential occupancies, first
multiply its four already-thinned local matrices. Take an entrywise alpha
power, then average the tuples conditional on their total occupancy. The
message prefactor is C(512,q)*beta^(q alpha); the group-subset union is not
raised to alpha.

The actual24 cached G4 gate gives the following selected floating margins:

| q | theta | alpha | Original birth basis | Potential-birth basis |
|---:|---:|---:|---:|---:|
| 4 | 0.01 | 0.7 | 65.794678 | 65.808729 |
| 16 | 0.06 | 0.35 | 15.661667 | 16.426915 |
| 64 | 0.3 | 0.35 | 123.647248 | 132.500268 |

The potential-birth basis replaces the mixture of physical-occupancy birth
families by one family per potential occupancy. The checked identity
Tnew_j Q=Q Told_j, with e0 Q=e0 and Q1=1, preserves the unpowered comparison
moment. Fractional path bounds can differ because their hidden-state
decomposition changes. This is not a new encoder or a claim that the filled
comparison matrices describe the exact physical state law.

The runs check every alpha=1 global control and replay selected winners in
the all-log backend. Source pins and local intertwining checks pass. Receipts
are `actual24_fractional_g4_v1.json` and
`actual24_fractional_g4_rebased_v1.json`; drivers are `cached_fractional24.py`
and `cached_fractional24_rebased.py`.

The subsequent bounded refinement gives these potential-birth-basis winners:

| q | theta | alpha | Margin (bits) |
|---:|---:|---:|---:|
| 8 | 0.03 | 0.5 | 30.699442 |
| 16 | 0.06 | 0.4 | 17.354183 |
| 32 | 0.12 | 0.3 | 34.045246 |

For q=16, the grid was theta in {0.04,0.05,0.06,0.07,0.08} and alpha in
{0.25,0.3,0.35,0.4,0.45}. For q=8 it was theta in {0.02,0.03,0.04};
for q=32 it was theta in {0.09,0.12,0.15}. Both latter points used alpha
in {0.3,0.4,0.5}. Source pins, alpha=1 controls, and selected all-log replays
passed for both the original and rebased calculations.

The added receipts are `actual24_g4_q16_refined_v1.json`,
`actual24_g4_q16_rebased_refined_v1.json`, `actual24_g4_q8q32_v1.json`, and
`actual24_g4_q8q32_rebased_v1.json`. These results singled out an actual24
G8 calculation at theta=0.06,alpha=0.4, targeting q=16. The selected G4
point was about 23 bits below the desired per-occupancy 40-bit margin.
These selected G4 results do not establish a complete certificate.

Before a larger fixed group was evaluated, the root requested an event-aligned
control. It groups two successive nonempty potential steps, including their
preceding zero runs, and retains a terminal zero- or one-event block. The
actual24 q=16 result at theta=0.06,alpha=0.4 was 14.715358493 bits, below
the rebased G4 result. Full regional alpha=1 agreement and source validation
passed. The receipt is `actual24_event_q16_v1.json`, generated by
`cached_event24.py` using the independently tested `event_aligned.py` helper.
No larger event-aligned search was run. This selected point does not motivate
expanding that branch without another concrete tightening.

## Final bounded G8 check

The root approved exactly one actual24 rebased G8 macro calculation at
theta=0.06,alpha=0.4. Its q=16 margin is **26.752381158 bits**, a gain of
9.398198 bits over rebased G4. It remains below the per-occupancy 40-bit
target. The same saved macro was evaluated for every q from 1 through 512;
it is positive on q=12..20, with no occupancy reaching 40 bits. These are
results for this one parameter point, not optimized bounds for the other
occupancies.

The frozen helper enumerated 43,046,721 ordered tuples. Its simultaneous
alpha=1 macro sum agreed with ordinary grouping to relative error
7.09e-14. Independent all-log comparisons covered q=0..512: the maximum
regional discrepancy was 2.50e-12, and the maximum global discrepancy was
6.50e-11. Strict NumPy underflow checks passed during the macro calculation.
The conservative positive-product logarithmic lower bound for the supplied
floating matrices was -148.25, well above the binary64 normal threshold.
These checks do not turn the floating inputs or powers into outward endpoints.

`cached_long24.py` wrote `actual24_rebased_g8_theta006_v1.json`, including
the powered macro, all occupancy results, the original actual24 map identity,
the rebasing checks, and successfully verified source hashes. Total elapsed
time was about 91.54 seconds. The helper's 64 MiB setting is a nominal chunk
heuristic, not a measured peak-RSS bound; overlapping temporaries and runtime
allocations are additional.

No second G8 point, new state census, encoder change, or benchmark was run.
The bounded cycle stops here. The remaining gap warrants a concrete local
bound or map improvement before expanding the expensive G8 grid. A complete
certificate also needs the occupancy union budget and outward evaluation;
the per-occupancy 40-bit threshold alone is not sufficient.

`NUMERICAL_CERTIFICATION.md` describes an outward-verification design for
the eventual fine-grouped bound. It separates exact local arithmetic,
audited positive floating arithmetic, and final outward interval evaluation.
It does not certify any of these NumPy proposal receipts.
