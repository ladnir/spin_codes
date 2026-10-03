# Width-eight packets: matched implementation and outward verification

This checkpoint tests whether byte-sized routing packets can retain a complete
10%-distance proof near the 100-microsecond K16 target. The faithful isolated
encoder measures **98.1625 microseconds**, including inner, route, and outer.
The complete outward certificate gives **68.893748662558 bits** of whole-code
margin at 10% distance. It covers every nonzero occupancy, not a sampled grid.
The earlier floating proposal is not used as a certificate input.

The [performance record](PERFORMANCE.md) retains every fresh-seed run median.
The [implementation](../../../../spin/experiments/packet8_wider24/README.md)
is research code, not a promoted backend. Existing certified four-bit
baselines and production code remain unchanged by this investigation.

## Construction and claimed event

Fix K=65536 and N=131072. Divide the message into 256 groups of 256 bits.
Each group contains eight parallel GF(16) RS[16,8] rows. At each of their
16 output positions, collect the eight four-bit symbols into one 32-bit symbol.
Independently randomize each symbol by the binary adjoint of multiplication
by a uniformly sampled nonzero GF(2^32) scalar.

Split each randomized symbol into four consecutive bytes. Independently
shuffle each group's 64 bytes, assigning one to every region. Within each
of the 64 regions, independently shuffle the 256 group slots. Thus routing
moves eight binary coordinates together, not eight 128-bit payload lanes.

The inner processes eight bytes per step, using a 24-bit state. Over
GF(256) in the AES polynomial basis, write the state as s=(a,b,c) and define

\[
(As)_h=a+hb+h^2c,\qquad
Cx=\left(\sum_hx_h,\sum_hhx_h,\sum_hh^2x_h\right),\qquad h=0,\ldots,7.
\]

A step emits y=x+As and updates s'=M s+Cx. At setup, independently sample
each M as the binary adjoint of multiplication by a nonzero scalar in
GF(256)[z]/(z^3+z+1). State starts at zero, continues across all regions,
and has no final flush. All setup choices are fixed for subsequent messages.

Let E denote the resulting binary encoder. The bad setup event is

\[
\exists m\in\mathbb F_2^{65536}\setminus\{0\}:\quad
\operatorname{wt}(E(m))\le13107.
\]

The certified upper bound on this probability is the exact dyadic number
`0x11390c2f99a70907b2cb6a53aa73f5deffcdcdba114b * 2^-241`, which is below 2^-40.
Excluding this event gives minimum distance at least 13108,
which exceeds 10% of N. This is a whole-code claim, not a per-message tail
bound or the measured distance of one benchmark seed. The numerical argument
uses the explicit binary64 arithmetic contract in [LOCAL_ARITHMETIC.md](LOCAL_ARITHMETIC.md).

The scalar maps send every fixed nonzero input uniformly to a nonzero
output. This suffices for the proof; a uniformly sampled full binary matrix
is unnecessary. The [implementation audit](IMPLEMENTATION_AUDIT.md) checks
the adjoints, field bases, literal maps, chronology, and byte layout.

## Why the proof now closes the difficult cases

The original width-eight experiment retained a smaller outer and a 16-bit
state. Its ordinary first-moment bound counted rare routing concentrations
once for each message. That bound was far too large. A larger return
denominator alone did not repair it.

The current construction uses 64 regions and the actual 24-bit maps above.
For q>=2 active groups, the proof first conditions on their potential
occupancies in every inner step. It then clips the conditional first-moment
bound before averaging over routing. A fractional power implements this
clipping: min(1,x)<=x^alpha for 0<alpha<=1.

Four chronological local matrices are multiplied before taking their
entrywise fractional power. Averaging the fine occupancy tuples happens
afterward. The resulting regional matrices retain the inner state through
all 64 regions. Neither the route probabilities nor the choice of q groups
is raised to alpha. The [mathematical audit](../iteration4/AUDIT_WIDER.md)
details this conditioning argument and the wider outer's pointwise envelope.

The q=1 contribution uses exact expected outer shell counts instead of that
uniform envelope. Its outward bound is **110.785884909590 bits**. The weakest
previously identified q=3 witness retains **68.894429041914 bits** outward;
the floating counterpart was 68.894675 bits. Every other integer q is
covered, and their exact union gives the margin stated above.

The receipt `whole_wider24_v1.json` has SHA256
`505152bf4593cb506880ba195fbce955819dbe9e4aa2ccd5ccc4c6b60b55595f`.
Its 68 authenticated source/file pins include the local and global proof
producers, all numerical witness inputs, and the implementation's local
include dependencies. The full global replay took about 180 seconds;
the 12 additional local/macro witnesses took about 224 seconds with one
shared 2^24-state census.

An independent global replay uses 64 sequential row-vector updates instead
of binary matrix powers. It recomputed all 255 selected q>=2 bounds and the
q=1 term in about 34 seconds. Its exact union also passes 2^-40, with the
same displayed margin. The largest individual margin difference is below
2.9e-14 bits. This replay shares the authenticated local/macro endpoints and
positive arithmetic primitives; it is not a second local census implementation.
All 71 pins of `whole_wider24_independent_v1.json` authenticate; its SHA256 is
`0b295ec3e267929424f7079407dad57f70ef1bb6c50462b0a313a48d90b45cf8`.
All 36 proof tests and six native encoder configurations pass. The frozen
four-bit K16 checkpoint's 220 source pins and archived receipt still match.

## Numerical verification

The [certificate plan](CERTIFICATE_PLAN.md) separates the proof identities
from floating arithmetic assumptions. The implementation follows these steps:

1. `local_outward24.py` computes exact profile polynomials and bounds every
   Walsh-transform and profile-summation error. A fresh census covers all 2^24 states.
2. `macro_outward.py` thins active bytes, defines an exactly normalized rational
   birth basis, and verifies the change-of-basis identity before rounding upward.
3. `outward_positive.py` uses positive rounding bounds and integer-certified
   concave tangents. It does not trust a library fractional power as an endpoint.
4. `scaled_positive.py` carries an integer scale through positive matrix
   operations, including explicit allowances for gradual underflow.
5. `global_outward.py` evaluates exact hypergeometric placement weights and
   the continuous-state product across 64 regions.
6. `q1_outward.py` handles the single-group support polynomial and exact shell counts.
7. `certify_wider24.py` covers q=1,...,256 and compares the exact sum of
   dyadic endpoints with 2^-40. Displayed bit margins do not decide acceptance.

All 13 local witnesses have matching positive supports. Their largest
relative interval width is below 4.07e-7. The runtime check rejects FTZ/DAZ;
the arithmetic contract requires binary64 round-to-nearest operations.
Exact small-field and small-matrix tests check the bound directions,
structural zeros, underflow, rebasing, chronology, and final union.

## Reproduction

Run from the repository root. Python requires NumPy and SciPy through the
retained mathematical source dependencies. Native performance measurements
are separate from these proof computations.

```powershell
python -B -m unittest discover -s research/workstreams/packet8_codesign/iteration5 -p 'test_*.py' -v
```

The first local parameter is the exact dyadic approximation selected near
exp(-.03). Its exact value, not an exponential approximation, enters the proof:

```powershell
python -B research/workstreams/packet8_codesign/iteration5/local_outward24.py --z 8740996286544847/9007199254740992 --output research/workstreams/packet8_codesign/iteration5/local_theta003_v1.json
python -B research/workstreams/packet8_codesign/iteration5/macro_outward.py --local research/workstreams/packet8_codesign/iteration5/local_theta003_v1.json --alpha 1/2 --output research/workstreams/packet8_codesign/iteration5/macro_theta003_alpha1_2_v1.json
python -B research/workstreams/packet8_codesign/iteration5/local_batch.py --directory research/workstreams/packet8_codesign/iteration5/cache_v1
```

Those output paths must be fresh, except that the batch reuses authenticated
existing receipts. Do not overwrite a source-pinned receipt after changing
its producer. Use fresh versioned paths for subsequent experiments.

```powershell
$packetMacros = @('research/workstreams/packet8_codesign/iteration5/macro_theta003_alpha1_2_v1.json') + @(Get-ChildItem -LiteralPath 'research/workstreams/packet8_codesign/iteration5/cache_v1' -Filter 'macro*.json' | Sort-Object Name | ForEach-Object { $_.FullName })
python -B research/workstreams/packet8_codesign/iteration5/certify_wider24.py --macros $packetMacros --output research/workstreams/packet8_codesign/iteration5/whole_wider24_v1.json
python -B research/workstreams/packet8_codesign/iteration5/verify_whole.py --receipt research/workstreams/packet8_codesign/iteration5/whole_wider24_v1.json --output research/workstreams/packet8_codesign/iteration5/whole_wider24_independent_v1.json
```

Generated numerical receipts and timing logs remain ignored. The scripts,
tests, construction notes, and measured source identities are the reproducible
record. No benchmark should run concurrently with another benchmark.

## Engineering interpretation and next step

The new construction is about 7.4% slower than the uncertified smaller-outer
byte kernel, but about level with the certified four-bit control. This is a
new proof/performance checkpoint, not evidence of a large speedup.

The diagnostic split attributes most added latency to the wider outer:
about 52.0 microseconds versus 46.4 for the earlier byte kernel. Inner and
routing rise only from about 44.9 to 45.9 microseconds. The optimized inner
has no observed register spills. Preserve this matched implementation and
certificate as a research checkpoint. Subsequent speed work should target
the outer's 32-bit randomization and parity evaluation.
