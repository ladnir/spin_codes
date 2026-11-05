# Paper and research

Start with the [SPIN manuscript](paper/README.md). The production library is
separate, in [`../spin/`](../spin/README.md), and does not depend on this tree.

The [one-route coupled-scan study](workstreams/one_route_coupled_scan/README.md)
preserves the fixed field-quartet construction and its global inverse/RS-rank
proof work. The target is rate 1/2, 45% distance, and a 40-bit setup-failure
margin at K20. Several complete support classes are bounded, but the low-rank
remainder and whole-code certificate remain open. The goal is paused; the
[checkpoint](workstreams/one_route_coupled_scan/GOAL.md) records the latest
restricted diagnostics and the unreported full-population control.

The new [characteristic-two field family](workstreams/char2_field_spin/README.md)
explores genuine GF(2^m)-linear SPIN codes separately from the binary profiles.
Its [certified GF256 high-rate winner](workstreams/char2_field_spin/GF8_HIGH_RATE_WINNER.md)
is preserved. The earlier [full-alphabet distance ladder](workstreams/char2_field_spin/HD_FIELD128_LADDER.md)
lifts the leading 35%-distance paired-route design to GF(2^128), then seeks
higher distance at K20 and rate 1/2. Its literal 35% lift has a complete
108.117-bit ideal-setup certificate, but measures about 22.7 ms. A separate
gamma30 construction now has a [complete 40% certificate](workstreams/char2_field_spin/HD_FIELD128_GAMMA30_CLOSED.md)
with 49.171 bits of ideal-setup failure margin. Its 274 event bounds cover
all 8,192 active-group weights. This construction is not yet timed;
45% and higher-distance rungs remain open.
The GF256 3.80-ms checkpoint remains the faster certified reference.

The earlier [high-distance blockwise restart](workstreams/char2_field_spin/HD_BLOCKWISE_RESTART.md)
targeted rate 1/2, 40% distance, and sub-5ms forward encoding at K20.
A [sixteen-state refinement](workstreams/char2_field_spin/HD_SUMRANK_FEEDBACK_CERTIFICATE.md)
now has a complete 40%-distance certificate with more than 59.87 bits of
setup-failure margin. The [encoding refinement](workstreams/char2_field_spin/HD_ENCODING_OPTIMIZATION.md)
now reaches about 6.1 ms using a [shared-diagonal RS inner](workstreams/char2_field_spin/HD_SHARED_DIAGONAL.md),
with a complete 40%-distance certificate above 49.33 bits. A newer
[alternating-feedback refinement](workstreams/char2_field_spin/HD_BLOCK_FEEDBACK.md)
measures about 5.5 ms with streaming output and has complete directed
coverage above 51.04 bits, independently replayed at 384-bit precision.
The [rank-two alternative](workstreams/char2_field_spin/HD_BLOCK_EMISSION.md)
also passes independent replay, with more than 50.99 bits and about
5.4--5.5 ms full encoding. A period-four feedback schedule has a
complete certificate above 56.65 bits at similar speed, independently
replayed at 384-bit precision. The independent-map
[6.3-ms checkpoint](workstreams/char2_field_spin/HD_SYSTEMATIC_MDS.md) is preserved.
The [optimization checkpoint](workstreams/char2_field_spin/HD_OPTIMIZED_CHECKPOINT.md)
records the validated approximately 5.5-ms frontier. Sub-5ms encoding
at 40% remains unmet; no library promotion has been made. At 35%, the
[grouped small-RS outer](workstreams/char2_field_spin/HD_GROUPED_OUTER.md)
now has a GF256 candidate with a complete 72.3975-bit certificate and
about 4.0 ms K20 forward encoding. It combines invertible Toeplitz column
mixers, compact byte coefficients, and an exact-map striped inner schedule.
Repeated shared-host runs compare against the preceding 4.4--4.5 ms version
and the approximately 4.76 ms RS256 control; production is unchanged.
The newer [paired-route checkpoint](workstreams/char2_field_spin/HD_SHIFTED_CLOSED_CERTIFICATE.md)
now has a complete 35%-distance certificate above 55.09 bits. Matched
serial timings give 3.80 ms versus 4.05 ms for the retained reference;
native output layout alone gives 3.87 ms. The proof covers every nonzero
active-group weight, with matching 256/384-bit replays and an independent
complete probability sum. The next structural target is cheaper column
mixing; production and paper results remain unchanged.
The [randomized-outer baseline](workstreams/char2_field_spin/HD_RANDOMIZED_OUTER.md)
restores the binary proof's measure comparison at a measured0.48-ms cost:
7.02ms versus6.53--6.55ms for its fixed-outer control. The
[exact sparse ladder](workstreams/char2_field_spin/HD_RANDOMIZED_OUTER_SPARSE_LADDER.md)
now covers every word with one through five active rows at a 56.51-bit
union margin. Additional quiet and dense class bounds are recorded;
the complete certificate for that four-state native baseline remains open.
The earlier [native research control](workstreams/char2_field_spin/HD_BLOCKWISE_NATIVE.md)
takes about6.35ms; the [quiet-subcode theorem](workstreams/char2_field_spin/HD_BLOCKWISE_QUIET_SUBCODE.md)
records which cheaper state choices cannot attain the target distance.
The current [shared local mixing circuit](workstreams/char2_field_spin/HD_BLOCKWISE_SHARED_FACTOR.md)
targets the structural bottleneck with exact local rank proofs. Its
[native comparison](workstreams/char2_field_spin/HD_BLOCKWISE_SHARED_NATIVE.md)
shows that59 packed GFNI calls do not beat82:7.65ms versus6.54ms; a
whole-vector rewrite takes7.00ms. It has no whole-code certificate.
The measured6.35-ms implementation remains a frozen reference.
The preceding
[accumulator detour](workstreams/char2_field_spin/ACCUMULATOR_DETOUR_CHECKPOINT.md)
is parked and archived, not promoted.

The [width-eight packet study](workstreams/packet8_codesign/README.md) is
parked: it has a complete K16 proof, but only modest small-size gains and
a K20 regression. Its sources and results are retained separately from
production. The certified four-bit design remains the main path; see the
[K16 co-design checkpoint](workstreams/k16_codesign_100us/README.md).

| Directory | Contents |
|---|---|
| [paper/](paper/README.md) | Current manuscript, figures, and supporting paper scripts. |
| [artifact/](artifact/README.md) | Paper artifact packaging and author-side reproduction tools. |
| [workstreams/](workstreams/) | Proof development, experimental implementations, and selected results. |
| `scripts/`, `constructions/`, `explorations/`, `bch_spectrum_work/` | Supporting research and historical experiments. |
| `BA_paper/`, `enumerator_paper/` | Related manuscripts. |
| [Expander codes](https://github.com/ladnir/expander_codes) | Separate repository for the expander-code manuscript, source, and certificates. |
| `collaboration/` and top-level Markdown notes | Historical plans, audits, and handoffs. |

## Running historical commands

The research tree was moved together, preserving its internal relative paths
and the bytes of frozen sources and manifests. Run commands in its historical
notes from this directory, not the repository root:

```sh
cd research
python -B -m unittest discover -s artifact -p 'test_*.py'
```

Old absolute checkout paths and generated build directories may need to be
reconfigured. Historical GitHub links using paths at the former repository root
can be read at commit `557e40119f58323ad4e10238532d131be0a17f02`; current paths
have the `research/` prefix. Frozen manifests retain their original relative
names and hashes.

See the [publication policy](GITHUB_PUBLISH_POLICY.md) before adding research data.
Do not commit experiment dumps or build products.
