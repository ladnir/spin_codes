# Paper and research

Start with the [SPIN manuscript](paper/README.md). The production library is
separate, in [`../spin/`](../spin/README.md), and does not depend on this tree.

The new [characteristic-two field family](workstreams/char2_field_spin/README.md)
explores genuine GF(2^m)-linear SPIN codes separately from the binary profiles.
Its [certified GF256 high-rate winner](workstreams/char2_field_spin/GF8_HIGH_RATE_WINNER.md)
is preserved. The active [high-distance blockwise restart](workstreams/char2_field_spin/HD_BLOCKWISE_RESTART.md)
targets rate 1/2, 40% distance, and sub-5ms forward encoding at K20.
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
remains unmet; no library promotion has been made.
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
