# Paper version 2: plan (draft; nothing run yet)

## Scope decisions (user, 2026-11-11)

- Paper keeps four constructions: Accumulator SPIN and Random SPIN (warmups),
  Scalable SPIN (theory only, no performance), Turbo SPIN (all measurements).
- Remove BCH-IMT everywhere. An explicit "Version 2" note at the start says
  version 1 also described a BCH-IMT instantiation.
- Turbo is Lean-certified at K = 2^14, 2^16, ..., 2^26. The implementation
  is used anywhere in that range and may run beyond it. Paper wording must
  keep the theorem at the certified lengths; other lengths are use, not proof.
- Applications (Silent OT regular + stationary, end to end; standalone PCS;
  Flock) move to Turbo. The transposed-comparison table keeps BAA/EC/RAA and
  re-measures the SPIN row as Turbo inside that harness.
- Rate 1/4: screen easy first-moment bounds for a Turbo variant numerically
  (no Lean), paper-level argument only.

## Code facts that shape the work

- Hypercat main `82015655` selects `PacketRsT64S15K16` (paired-s15, four
  RS[16,8] rows, GF(2^16) multipliers, K = 2^16 only). That is not the
  paper's Turbo (`PacketRs8T64`: eight RS[8,4] rows, GF(2^32) columns, s22).
- `PacketRs8T64` forward encoding is the scalar fallback (Code.h:122); only
  its 128-bit transpose is optimized. PCS and Flock need forward encoding.
- libOTe pins spin_codes `f2010e03` and calls BCH-IMT only. Its stationary/regular
  harness was deleted in libOTe `5d976988`; restore from `c575f160`.
- spin_codes main `f5f232d8` and this branch (`42dd4bac` s22 change) have
  diverged; the applications need one pinned revision containing both.

## Work order

1. Paper cut: remove BCH-IMT sections, figures, certificates, rate-1/4 BCH,
   and application text tied to BCH-IMT; add the version-2 note.
2. Code: merge spin main into the Turbo branch; add an optimized PacketRs8T64
   forward kernel (pending decision); add Turbo to libOTe (regular and
   stationary) and to Hypercat PCS/Flock.
3. Rate-1/4 Turbo screen (Python, local or Peach as a build-class job).
4. All builds on Peach, then one serial benchmark campaign (core 15, flocks,
   no concurrent builds): Turbo forward/transpose, transposed comparison,
   Silent OT regular + stationary, PCS, Ligerito, Bolt, Flock.
5. Update tables and text from receipts.
