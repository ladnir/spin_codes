# Parked small-state research checkpoint

Parked on 2026-10-21 at the user's request. The active task is now to push
the certified GF(2^128) RS/stratified-route/Toeplitz-window construction
toward 50% distance. Completing another 35% certificate is not an objective.
All small-state sources and proof notes remain in place; nothing is deleted
or promoted. This is a local research checkpoint, not a Git commit.

## Preserved construction and evidence

The [specification](SMALL_STATE_PROOF_GATE.md) retains four RS[64,32] rows,
GF256 Toeplitz column mixers, the paired regional route, and 64-input/16-state
maps. Its two refresh diagonals use independent nonzero GF(2^32) coefficients.
Setup is fixed for all messages. There is no competitive timing for this
exact proposal and no whole-code distance certificate.

At 35% distance, directed results cover active-group counts 1..128,
2048..2073, and 4096..4108. The first range has combined failure below
2^-44; each other range is below 2^-100. The remaining ranges are open.
The g=1024 sixteen-step score near -450.52 is FLOAT-only, not certified.
See the [coverage record](SMALL_STATE_EXACT_ROUTE_GATE.md) and
[progress ledger](GOAL.md).

Selected source SHA-256 checkpoints:

- `theory/verify_small_state_grouped_interval.py`:
  `79644db7e509a017340918e847b0685f51fc8910c72edb84bcbbf7fd385711ff`.
- `theory/probe_small_state_window_route.py`:
  `0a7f941e8c0e94c937bc7f5982622eab1ac85a93da9dec00baaf313da0e9ac8f`.

The exact interval verifier uses directed arithmetic and complete run counts.
The FLOAT window helper still needs certified upward grid assignments and
positive arbitrary-exponent arithmetic before any promotion.

## Why the present proof is not the high-distance path

The [fourth-root bridge](SMALL_STATE_RANK_BRIDGE.md) reduces the full-field
existence count to a GF256 auxiliary count. For all L=8192 source groups,
every allowed input batch is full. Put t(z)=(1+255z)/256 and N=2^21.
The exact two-state operator contains the first-birth/all-live path t(z)^N.
With beta=256^160/255^32 and kappa=256^256/beta, the current joint-density
bound is therefore at least

    [(2/255) beta^L kappa^(L/256) z^(-e) t(z)^N]^alpha.

Here 2 is the verifier's bound on Gamma. Optimizing this surviving term
uses z=e/[255(N-e)]. Its exponent changes sign at approximately 37.8596%
relative distance. This is a barrier of this bound, not an upper bound on
the actual code distance. More grouping or numerical precision cannot
repair this earlier loss of the larger-field exponent.

The cost review also needs the right comparator. The direct GF32 refresh
uses 32 GFNI vector multiplications and 24 additional word rotations per
64-record step, versus 8 multiplications and no rotations in the retained
GF256 refresh. Its prepared coefficient layout uses 16 MiB at K20, versus
1 MiB for the compact GF256 keys. These are operation and layout counts,
not timings. The preserved 3.798 ms GF256 result remains a separate code.

Resume this work only if it supplies an explicit route to the near-50%
target, rather than merely extending its 35% source-class coverage.
