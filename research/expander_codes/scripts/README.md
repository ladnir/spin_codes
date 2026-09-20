# EA finite certificates

`prime_field_ea_diagnostic.py` evaluates the proved asymptotic labeled
finite-field Expand--Accumulate conditions.  It reports the $p$-ary GV distance,
the sparse accumulator rate, and the limiting Bernoulli row-weight constant.
It is an asymptotic diagnostic, not a finite-length certificate.

`prime_field_regular_ea_diagnostic.py` evaluates representative exact-shell,
Poissonized, and dense Fourier terms for the labeled region-stratified regular
ensemble.  Its output is a floating-point parameter probe; it does not cover
all message supports or use interval arithmetic.

`ea_certificate.py` generates and verifies rigorous finite first-moment bounds
for the binary Expand--Accumulate ensemble.

Generation uses SciPy only to select Chernoff markers. Verification treats
those markers as fixed decimal values and recomputes the bound with
outward-rounded Arb balls from `python-flint`. The floating-point optimizer is
therefore outside the trusted calculation.

To verify the checked-in 20-bit table serially, run:

```powershell
python ea_certificate.py table `
  certificates\ea-k20-r5-d005-w50-s20.json `
  certificates\ea-k25-r5-d005-w56-s20.json `
  certificates\ea-k30-r5-d005-w62-s20.json
```

The certificates cover

```text
k = 2^20, 2^25, and 2^30
n = 5k
distance cutoff = 0.05n
expected Bernoulli-expander row weights = 50, 56, and 62
failure target = 2^-20
```

The verifier evaluates message weights 1 through 128 with the exact two-state
generating matrix. It covers all remaining weights with 211 monotone spectral
blocks. It rejects missing, overlapping, or out-of-order ranges.

To regenerate the marker choices, run:

```powershell
python ea_certificate.py generate `
  --row-weight 50 `
  --output certificates\ea-k20-r5-d005-w50-s20.json
```

Regeneration can change the stored decimal markers slightly across SciPy
versions. Any generated file remains sound if the Arb verifier accepts it.

`exact_row_ec_bounds.py` is the corresponding floating-point diagnostic for
wrapped Expand--Convolute. It composes the positive exact-row shell recurrence
with the exact two-marker convolution transfer matrix. For example,

```powershell
python exact_row_ec_bounds.py --k-log2 20 --row-weight 15 --memory 21
```

`exact_row_ec_certificate.py` verifies all message weights with
outward-rounded Arb balls. Small weights use the positive shell recurrence,
intermediate weights use uniform conditioned-Bernoulli transfer blocks, and
dense weights use the invertibility/L2 bound. The three certified parameter
sets are:

```powershell
python exact_row_ec_certificate.py `
  --k-log2 20 --row-weight 15 --memory 21 --exact-limit 59
python exact_row_ec_certificate.py `
  --k-log2 25 --row-weight 17 --memory 21 --exact-limit 240
python exact_row_ec_certificate.py `
  --k-log2 30 --row-weight 20 --memory 21 --exact-limit 529
```

At rate `1/5` and relative distance `0.05`, their verified failure bounds are
`5.1603074722e-7`, `4.2598255031e-7`, and `4.8307909353e-8`. These correspond
to `20.8860396`, `21.1627023`, and `24.3031653` certified bits.

`wrapped_ec_asymptotic.py` evaluates the closed sparse exponent from the
asymptotic wrapped-EC theorem. For the paper's primary parameters, run:

```powershell
python wrapped_ec_asymptotic.py --memory 21 --relative-distance 0.05 --rate 0.2
```

This gives sparse exponent `0.899585774977`, density threshold
`C > 1.111622735504`, and a positive Gilbert--Varshamov margin of
`0.356002501102` nats per output coordinate.

`prime_field_biregular_ec_diagnostic.py` explores the nonwrapping
prime-field EC constraint-trace bound with a two-sided regular expander. The
convolution samples a fresh feedback vector from `F_p^m` at each output
position. The trace counts every zero event that imposes an equation in fresh
feedback coefficients or fresh edge labels; it does not apply the projective
cap to an already averaged Markov probability. For a floating-point scan,
run:

```powershell
python prime_field_biregular_ec_diagnostic.py --degree 28 --memory 3 --exact-limit 8
```

`binary_biregular_diagnostic.py` explores two-sided regular binary EA and
wrapped EC.  Small supports use the exact parity-shell law.  Larger supports
use a positive coefficient saddle for the even and odd group polynomials.
For the first rate-half candidates, run the commands serially:

```powershell
python binary_biregular_diagnostic.py --code ec --k 1048567 --left-degree 26 --right-degree 13 --cutoff 230729 --memory 9
python binary_biregular_diagnostic.py --code ea --k 1048575 --left-degree 62 --right-degree 31 --cutoff 230729
```

Binary right degree must be odd.  If it is even, the all-one message maps to
zero and the code is not injective.

The selected rate-half binary EC profile has a checked-in marker certificate.
First validate its schema and support partition with the standard-library
checker.  Then evaluate the probability bound with outward-rounded Arb
arithmetic:

```powershell
python check_binary_biregular_ec_certificate.py ../results/binary_biregular_ec_rate_half_d10_m15_gv.json
python binary_biregular_ec_certificate.py verify ../results/binary_biregular_ec_rate_half_d10_m15_gv.json
```

The fixed target has left/right degrees `10/5`, memory `15`,
`k=1,048,575`, and cutoff `230,729`.  Exact transfers cover supports `1..16`.
Positive coefficient blocks cover the two outer ranges.  A degree-five
local-limit bound covers the central range.  At 192-bit precision, the total
is `1.6474999656e-10`, or `32.4990025` failure bits.

The same verifier handles all odd right degrees.  Two additional frozen
profiles use `k=1,048,635`, `n=2,097,270`, and cutoff `230,742`:

```powershell
python check_binary_biregular_ec_certificate.py `
  ../results/binary_biregular_ec_rate_half_d14_m6_gv.json
python binary_biregular_ec_certificate.py verify `
  ../results/binary_biregular_ec_rate_half_d14_m6_gv.json
python check_binary_biregular_ec_certificate.py `
  ../results/binary_biregular_ec_rate_half_d18_m4_gv.json
python binary_biregular_ec_certificate.py verify `
  ../results/binary_biregular_ec_rate_half_d18_m4_gv.json
```

Degrees `14/7` with memory `6` give `26.0875` failure bits.  Degrees `18/9`
with memory `4` give `21.6746` bits.  For a general odd right degree, the
central bound factors each parity-conditioned group polynomial into
Bernoulli generating polynomials and bounds their total variance.

The first command imports neither SciPy nor Arb.  It checks the parameter
identities, decimal marker domains, and a disjoint cover of supports `1..k`.
The second command consumes only the stored markers; it does not optimize
during verification.  To regenerate the marker file separately, run:

```powershell
python binary_biregular_ec_certificate.py generate `
  --output ../results/binary_biregular_ec_rate_half_d10_m15_gv.json
```

`binary_biregular_outer_scan.py` investigates the lower-degree `6/3`
profile.  It evaluates one exact support block with a shared output marker
and reuses the expensive uniform-slice transfer family:

```powershell
python binary_biregular_outer_scan.py `
  --support-start 80 --support-limit 95 --output-marker 0.99775366
```

The default memory is `80`.  The calculation is a binary64 diagnostic, not an
interval certificate.  Run separate blocks serially.

`binary_biregular_activation_diagnostic.py` marks transitions that leave the
inactive convolution state.  It explains why the ordinary coefficient saddle
is loose at small support, but it does not close the degree-`6/3`, memory-`60`
candidate.

`binary_biregular_shell_conditioned_diagnostic.py` conditions first on the
exact right-degree-three parity weight.  It can then split a fixed-weight
convolution slice by activation count, inactive duration, or the time of the
first inactive wait.  These positive splits diagnose the multimodal
coefficient distribution at memory `80`, but none closes the analytic
support-256 bound.  The numerical findings are summarized in
`../notes/BINARY_DEGREE_SIX_OUTER.md`.

`binary_biregular_hybrid_diagnostic.py` carries one logarithmic scale per
conditional uniform slice.  This prevents absolute underflow when the output
marker makes every entry of an otherwise regular transfer matrix smaller than
binary64 can represent.  For example:

```powershell
python binary_biregular_hybrid_diagnostic.py `
  --message-weight 320 --output-marker 0.9917376
```

This tool is also diagnostic.  A proof must reproduce the scaled recurrence
with outward-rounded interval arithmetic.

The same tool evaluates a scaled exact support block with one shared marker:

```powershell
python binary_biregular_hybrid_diagnostic.py `
  --support-start 256 --support-limit 383 --output-marker 0.992
```

`binary_biregular_ec_d6_certificate.py` and
`../results/binary_biregular_ec_rate_half_d6_m79_gv.json` certify degree
`6/3`, memory `79`.  Run the structural checker and numerical verifier
serially:

```powershell
python check_binary_biregular_ec_d6_certificate.py `
  ../results/binary_biregular_ec_rate_half_d6_m79_gv.json
python binary_biregular_ec_d6_certificate.py verify `
  ../results/binary_biregular_ec_rate_half_d6_m79_gv.json
```

The verifier uses entry-scaled exponent layers for exact supports `1..383`,
103 positive-coefficient blocks, 28 complemented central blocks, and a
separate all-one calculation.  At 192-bit Arb precision, the total bound is
`3.5584438609e-7`, or `21.42225019` failure bits.  The final exact block is
computationally expensive; run no other benchmark concurrently.

`binary_biregular_tail_scan.py` partitions a coefficient-bound tail into
multiplicative blocks and reports both the worst block and their sum.  Run it
serially; each block retries several optimizer starting points.  These marker
selection calculations remain separate from certificate verification.

`prime_field_biregular_ec_certificate.py` closes the rate-`1/2` finite result
at the floored `p`-ary GV cutoff for `p=2^127-1`, `n=2097144`, left degree
`28`, right degree `14`, and convolution memory `3`:

```powershell
python prime_field_biregular_ec_certificate.py verify ../results/prime_field_biregular_ec_p127_rate_half_d28_m3_gv.json
```

At 256-bit Arb precision, the certified ensemble failure probability is at
most `4.752556241e-9`, or `27.64864915` bits. Eight exact bands cover supports
`1..192`; 32 convex saddle blocks cover `193..k-1`; full support is evaluated
directly. Floating point selects the stored decimal markers but is not part
of verification.

The singleton-refined transfer closes two lower-degree points. Degree `26`
with memory `4` has `71.47243412` certified bits. Degree `24` with memory `6`
has `39.50030421` certified bits. Verify the latter certificate with:

```powershell
python prime_field_biregular_ec_d24_m6_certificate.py verify ../results/prime_field_biregular_ec_p127_rate_half_d24_m6_singleton_gv.json
```

For the degree-`24` certificate, 17 exact bands cover supports `1..700`, 30
convex saddle blocks cover `701..k-1`, and the verifier evaluates full support
directly. The interval `225..256` supplies the largest contribution.

Degree `22` with memory `12` has `27.29518303` certified bits.  Its verifier
separates traces with fewer equations than the message support.  A
singleton-free-region bound handles those traces through support `48`; exact
field and transfer bounds cover the complementary traces.  Verify it with:

```powershell
python prime_field_biregular_ec_d22_m12_certificate.py verify ../results/prime_field_biregular_ec_p127_rate_half_d22_m12_region_gv.json
```

The exact transfer covers supports `49..1400`, and 26 convex saddle blocks
cover `1401..k-1`.  The interval `481..520` supplies the largest contribution.
The verifier keeps an advisory content-addressed cache under
`../.verification_cache/d22_m12`.  The key includes the certificate, Arb
precision, verifier code, and python-flint version.  A normal rerun can use a
completed result immediately.  To recompute every interval without reading
or writing the cache, run:

```powershell
python prime_field_biregular_ec_d22_m12_certificate.py verify ../results/prime_field_biregular_ec_p127_rate_half_d22_m12_region_gv.json --no-cache
```

Use `--refresh-cache` instead to recompute the proof and populate resumable
per-band entries.  Cache entries accelerate local work; they are not part of
the trusted certificate.

The obsolete degree-`18`, memory-`1` and degree-`24`, memory-`1` certificates
used a probability transfer followed by a projective cap. That order is not
valid when different projective messages impose different feedback equations;
the files were removed rather than retained as apparent certificates.

`bernoulli_ec_certificate.py` verifies the same Bernoulli, memory-5,
zero-initial-state ensembles reported in Table 4 of the original EC paper:

```powershell
python bernoulli_ec_certificate.py --coefficient 3 --memory 5
python bernoulli_ec_certificate.py --coefficient 2.5 --memory 5
python bernoulli_ec_certificate.py --coefficient 2.3 --memory 5
```

For `k=2^20`, rate `1/5`, and relative distance `0.05`, the certified bounds
are `2.187816572e-10`, `9.485867880e-8`, and `1.124903039e-6`. The verifier
checks weights `1..32` individually, covers the intermediate range with
uniform transfer blocks, and uses invertibility from `k/4` onward. The dense
range is split into ten activation-probability blocks. Each block uses an
outward-rounded Hamming-ball bound, so the verifier retains the finite-length
prefactor that the entropy bound discards.

The optional `--cutoff` argument selects an exact integer output-weight
threshold. It takes precedence over `--relative-cutoff`; this avoids binary64
rounding when certifying a finite distance frontier.

The optional `--dense-start` argument moves the verified switch from transfer
blocks to invertibility blocks. The rate-`1/2` certificates use this switch
because the transfer relaxation becomes loose well before `k/4`.

For the finite GV-nearness point, run:

```powershell
python bernoulli_ec_certificate.py --k-log2 20 --rate-denominator 5 --coefficient 5 --memory 9 --cutoff 1274032
```

Here `n=5*2^20`, so the certified relative distance is
`1274032/n = 0.2430023193359375`. The all-weight bound is
`3.819697244e-7`, or `21.320038` bits. At rate `1/5`, the binary GV distance
is `0.2430038538089538`; the certificate reaches `99.99937%` of it. The
smaller cutoff `1272015` gives `29.140327` bits at `99.8411%` of GV.

For the rate-`1/2` Pareto profiles, run the commands serially:

```powershell
python bernoulli_ec_certificate.py --k-log2 20 --rate-denominator 2 --coefficient 3.15 --memory 5 --cutoff 230741 --dense-start 26363
python bernoulli_ec_certificate.py --k-log2 20 --rate-denominator 2 --coefficient 3 --memory 7 --cutoff 230741 --dense-start 26363
python bernoulli_ec_certificate.py --k-log2 20 --rate-denominator 2 --coefficient 2.9 --memory 9 --cutoff 230741 --dense-start 27462
python bernoulli_ec_certificate.py --k-log2 20 --rate-denominator 2 --coefficient 2.85 --memory 13 --cutoff 230741 --dense-start 28046
```

In increasing memory order, the certified `(expected row weight, bits)` pairs
are `(45.8517, 20.320353)`, `(43.6683, 20.853424)`,
`(42.2127, 20.778210)`, and `(41.4849, 20.675516)`. The relative distance is
`0.11002588272094727`, which reaches `99.99820%` of the rate-`1/2` binary GV
distance `0.11002786443829`.

## Aggressive paper-EC collision finder

`ec_collision_finder.py` attacks the archived one-pass systematic paper
ensemble at rate `1/2`. The old expander samples its five coordinates with
replacement. A repeated coordinate cancels over `F_2`, so approximately ten
of the `2^20` nominal weight-five rows have effective weight three. The finder
streams the expander rows and evaluates only these collision candidates under
one shared wrapping convolution. Its bit-sliced evaluator handles up to 63
candidates together.

Run a full paper-parameter trial with:

```powershell
python ec_collision_finder.py `
  --k-log2 20 --nominal-weight 5 --memory 25 --seed 1 `
  --output ..\results\ec_collision_k20_w5_m25_seed1.json
```

Here wrapping means that the oldest of the 25 feedback taps is fixed to one.
The block boundary is not cyclic. The NumPy generator samples the same
idealized ensemble but does not reproduce libOTe's AES-based seed stream.

The checked-in 16-seed sweep used the same command with `--trials 16`. It found
an average of `8.5` collision candidates per code. The mean best relative
weight was `0.1057841`; the smallest was `0.03628016`, from seed `7`. Thus this
linear-time collision attack alone reproduces the scale of the paper's
reported `0.1` empirical pseudo-distance.

`ec_region_attack.py` attacks the cleaner paper-like ensemble. It samples every
expander row as a uniform five-subset without replacement. The attack first
evaluates rows whose first edge occurs latest. It then buckets rows by region
multiset and solves the actual convolution-state equations at gaps between
selected regions.

Run the paper-scale experiment with:

```powershell
python ec_region_attack.py `
  --k-log2 20 --row-weight 5 --memory 25 --regions 16 --seed 1 `
  --max-buckets 16 --kernel-samples 32 --late-row-candidates 256 `
  --output ..\results\ec_region_k20_w5_m25_r16_seed1.json
```

For seed `1`, the best row has five distinct edges in the final sixteenth of
the parity block. Its systematic codeword has weight `21344` out of `2097152`,
or relative weight `0.0101776123`. The receipt includes the message row,
expander support, exact parity weight, and a SHA-256 digest of the packed
codeword.

The mechanism is independent of convolution randomness. An upper-triangular
convolution preserves every terminal suffix. For constant row weight `w`, a
linear scan finds a suffix row of length
`O(k^(1-1/w) log(k)^(1/w))` with high probability. The corresponding relative
codeword weight is `O((log(k)/k)^(1/w))`.

## Exact-row diagnostics

`exact_row_bounds.py` evaluates the exact-row ensemble through a positive
Hamming-shell Markov chain and a hypergeometric accumulator tail. It is a
floating-point diagnostic tool; rigorous exact-row certificates will use a
separate verifier.

`exact_row_certificate.py` is that verifier. It combines exact positive shell
recurrences, a conditioned-Bernoulli bound, and a dense Fourier/L2 bound. Its
floating-point optimizer only selects fixed decimal Chernoff markers; Arb
checks every resulting inequality.

To reproduce the exact-row table serially, run:

```powershell
python exact_row_certificate.py --k-log2 20 --row-weight 31 --exact-limit 7
python exact_row_certificate.py --k-log2 25 --row-weight 35 --exact-limit 16
python exact_row_certificate.py --k-log2 30 --row-weight 39 --exact-limit 48
```

## Streaming field heuristic audits

`streaming_ec_heuristic_audit.py` checks rank defects in compact regular
expanders. `streaming_ec_pair_collision_audit.py` counts the shared neighbors
that one projective row pair can cancel. Run the deployed-size topology audit
with:

```powershell
python streaming_ec_pair_collision_audit.py --seeds 1000
```

`streaming_ec_distance_ablation.py` compares regular topology, signed stripes,
shared labels, and periodic taps on reduced 26/13 instances. The exact tiny-code
control enumerates every nonzero ternary message:

```powershell
python streaming_ec_distance_ablation.py --seed 1
python streaming_ec_exact_ablation.py --seeds 32
```

`streaming_ec_mitm_search.py` performs a larger probabilistic search by joining
two half-message tables on selected zero coordinates.  The following reduced
rate-half comparison uses a collision-free prime-region construction:

```powershell
python streaming_ec_mitm_search.py `
  --left-degree 6 --right-degree 3 --region-size 7 `
  --memory 4 --prime 3 --zero-coordinates 16 --trials 32 `
  --variants proved_like collision_free_topology `
    collision_free_heuristic collision_free_edge_signs
```

The output includes the exact matrix rank and the lightest word found.  The
rank is a correctness check.  The reported weight is only an experimental
upper bound on distance.

The current libOTe heuristic regenerates taps at every position. The periodic
variant remains in these scripts as a negative control.

The same scripts compare full-field, random-sign, and all-one column labels.
The label ablation and the corresponding Goldilocks timings are recorded in
`../notes/STREAMING_LABEL_ABLATION.md`.  The C++ benchmark uses the streaming
heuristic with random coordinate signs by default:

```powershell
frontend_libOTe -bench -regularEc -goldilocks -t 10
```

Add `-reference` to benchmark the reference implementation. It materializes the
balanced assignments, independently generated nonzero edge labels, and
full-field convolution taps.

`regular_ec_streaming_spectral_audit.py` replays the exact signed-stripe
topology produced by libOTe for the deployed region size 80,659. Run the three
audits serially with:

```powershell
python regular_ec_streaming_spectral_audit.py
```

The script checks the order-1, order-79, and order-1021 components exactly
over the Goldilocks field. It also checks all 80,659 complex frequencies
numerically and all frequencies exactly over a compatible auxiliary prime.
The script does not check the six degree-13,260 Goldilocks factors exactly.
