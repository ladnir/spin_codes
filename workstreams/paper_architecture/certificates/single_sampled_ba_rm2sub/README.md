# One-sampled Golay--BA-3/RM2Sub-S19 certificate bundle

This directory is the paper-owned snapshot of the scalable Structured SPIN
certificate. It supports the following claim:

- rate exactly `1/2` at native admissible lengths;
- relative minimum distance `0.11` with probability `1-o(1)`;
- block schedule `B=(39/4) log_2(N)+O(1)`; and
- `O(N)` ordinary and transposed encoding work.

The certificate does not concern the frozen `N=2^21` ParityFanout member.

## Provenance

The files in the directory root were copied byte-for-byte from
`C:/Users/peter/.codex/worktrees/3061/permute_conv/workstreams/finite_asymptotic_theory/`.
The files in `linear_time_audit/` were copied byte-for-byte from
`C:/Users/peter/.codex/worktrees/3ef9/permute_conv/workstreams/linear_time_audit/`.

`SINGLE_SAMPLED_BA_RM2SUB_CERTIFICATE_MANIFEST.json` is the upstream binding
manifest. Its SHA-256 hash in this snapshot is
`cd5afa44c245bc46b6c55376d5f86908c5c0d069d17b1b5669c64a28ed466786`.
The manifest binds 23 theory files, five companion linear-time files, and
three frozen dependencies. `verify_imported_manifest.py` checks the relocated
files and the repository-relative frozen dependencies.

## Verification

Run the manifest check from the repository root:

```text
python workstreams/paper_architecture/certificates/single_sampled_ba_rm2sub/verify_imported_manifest.py
```

The publication-facing receipts are reproducible with the imported verifier
scripts. Run them sequentially. The principal checks are:

```text
python workstreams/paper_architecture/certificates/single_sampled_ba_rm2sub/certify_golay_ba_concave_majorant.py --output tmp/golay_ba3_concave_majorant.json
python workstreams/paper_architecture/certificates/single_sampled_ba_rm2sub/certify_golay_ba_rm2sub_joint_interval.py --majorant workstreams/paper_architecture/certificates/single_sampled_ba_rm2sub/golay_ba3_concave_majorant.json --output tmp/golay_ba3_rm2sub_joint_interval_d11.json
python workstreams/paper_architecture/certificates/single_sampled_ba_rm2sub/certify_rm2sub_dense_small.py --output tmp/rm2sub_dense_small_exact_d11.json
python workstreams/paper_architecture/certificates/single_sampled_ba_rm2sub/certify_rm2sub_uniform_fixed_occupation.py --output tmp/rm2sub_uniform_fixed_occupation_d11.json
python workstreams/paper_architecture/certificates/single_sampled_ba_rm2sub/certify_golay_ba_rm2sub_weight_coupled_fixed.py --majorant workstreams/paper_architecture/certificates/single_sampled_ba_rm2sub/golay_ba3_concave_majorant.json --output tmp/golay_ba3_rm2sub_weight_coupled_fixed_d11.json
```

The sparse verifier consumes those receipts. The imported JSON files record
the exact accepted results and remain the hash-bound publication artifacts.
