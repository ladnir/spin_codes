# Cheaper randomization and packed outer layouts

The selected K16 implementation takes **93.895 us**, compared with
**97.8175 us** for the frozen wider24 implementation in the same fresh-seed
campaign. This is a 4.01% time reduction. It computes exactly the same
binary code, so it retains the complete **68.893748662558-bit margin at
10% distance**. Setup and allocation are excluded; payload elements are
128 bits. This is an isolated research implementation, not a production
promotion.

The selected change flattens the nine-product GF(2^32) randomizer circuit
and shares its output sums. AVX-512 ternary XOR combines three inputs in
one instruction. Each payload half still uses nine GFNI byte products,
but now needs eleven logical instructions rather than fifteen. It keeps
the small nine-byte coefficient table, the RS circuit, route, state maps,
and setup distribution unchanged.

## Results

| Candidate | Full transpose, us | Complete margin at 10% | Decision |
|---|---:|---:|---|
| Frozen wider24 control | 97.8175 | 68.893749 bits | Preserve |
| Flattened field multiplier | **93.8950** | **68.893749 bits** | Select |
| Flattened multiplier, one identity symbol per group | 93.8205 | 68.893749 bits, by aggregate-measure equality | Retain; indistinguishable timing |
| Seven-diagonal MDS sandwich, three fixed GFNI maps | 96.8310 | 68.350804 bits | Provable alternative, not faster |

These are medians of eight run medians per candidate. Every run is retained
in [PERFORMANCE.md](PERFORMANCE.md). The 0.08% difference between the two
flattened variants is not sufficient to select a changed construction.

The packed flow already converts only at the external boundaries:

```
128-bit elements -> byte planes -> inner^T + route -> randomizer + RS^T -> output
                                  packed throughout
```

There was no intermediate unpack/repack left to remove. We tested reduced
outer workspace, parity-first output fusion, and six route-prefetch distances.
None beat the selected kernel. Half-payload layouts added masked-store work;
prefetch added instructions without improving the full call. The parity-first
layout also failed to improve the flattened randomizer. Smaller live workspace
does not by itself imply less traffic or a faster encoder.

## What changed in the proof exploration

The [structured-randomizer note](STRUCTURED_RANDOMIZERS.md) proves a useful
alternative to a full field scalar. It uses seven nonzero byte scalars around
a fixed four-byte MDS map. Every fixed nonzero symbol image is pointwise
dominated by

```
gamma = (2^32-1)/255^4 = 16843009/16581375
```

times the uniform-nonzero distribution. This also holds for the binary-adjoint
family used by the forward code. The resulting correction to the existing
certificate loses only about 0.543 bits. All 256 occupancy endpoints pass
after outward correction; the minimum-distance conclusion remains 13,108
at N=131,072.

The fixed MDS map is not free: the measured version uses seven random GFNI
affines plus three fixed ones, and a 224 KiB hot coefficient table instead
of 36 KiB. Replacing those fixed maps with shifts and masks is slower still.
Neither fewer random parameters nor fewer field multiplications alone predicts
the runtime.

A separate argument permits one preselected symbol randomizer per outer
group to be the identity. Global row symmetry of the eight parallel RS rows
makes the aggregate expected word measure exactly unchanged for the original
uniform-transitive family. This is **not** a fixed-message uniformity claim.
The equality holds before independent routing and fractional clipping, so
the existing whole-code certificate carries over. An independent reviewer
checked this conditioning point. The argument does not justify omitting
multiple symbols.

For the structured family, one identity symbol reduces the group domination
factor from gamma^16 to gamma^15. Its separate transported certificate gives
68.384761536057 bits. That combined structured/identity implementation was
not benchmarked in this campaign.

## Verification and reproduction

- Six native configurations pass: K=2048, 6144, and 65536, two seeds each.
  Each checks all 18 modes, literal scalar equality, 128 independent adjoint
  lanes, four alignments, in-place operation, guards, and route padding.
- Ten portable finite-field/algebra tests pass. They check the field towers,
  flattened multiplier, binary adjoints, all 69 MDS minors, the normalized
  diagonal family, and counterexamples to two overly light shortcuts.
- [verify_transport.py](verify_transport.py) independently checks every
  correction using exact rational powers, without the producer's floating
  or tangent arithmetic. Both 256-endpoint unions and their 71/73 source
  pins authenticate and pass 2^-40. The original 68-pin receipt also
  authenticates. This is an independent transport check, not a new local census.

From the repository root:

```sh
python -B -m unittest discover -s research/workstreams/packet8_codesign/iteration6 -p 'test_*.py' -v
python -B research/workstreams/packet8_codesign/iteration6/structured_randomizers.py --output research/workstreams/packet8_codesign/iteration6/mds_sandwich_transport_v1.json
python -B research/workstreams/packet8_codesign/iteration6/structured_omit_transport.py --output research/workstreams/packet8_codesign/iteration6/mds_sandwich_omit_one_transport_v1.json
python -B research/workstreams/packet8_codesign/iteration6/verify_transport.py research/workstreams/packet8_codesign/iteration6/mds_sandwich_transport_v1.json research/workstreams/packet8_codesign/iteration6/mds_sandwich_omit_one_transport_v1.json
```

Producer outputs must be fresh. Reuse the authenticated iteration5 inputs
or reproduce them using its [instructions](../iteration5/README.md).
Native build and benchmark commands are in the
[experiment README](../../../../spin/experiments/packet8_wider24_opt/README.md).
Generated receipts, binaries, and timing logs remain ignored. No frozen proof
source or production kernel was changed by this iteration.

## Next step

Retain the flattened exact-map kernel as the byte-packet implementation
checkpoint. For further speed work, optimize randomization and RS parity as
one circuit: the structured family provides a rigorously bounded alternative
if its fixed mixing can be made useful to the outer circuit. Do not select
another changed distribution from operation counts alone. Keep the present
kernel and the old frozen kernel in every matched full-call benchmark.
