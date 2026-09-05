# Riffle LDPCSplitState PacketMul g=4 t=256 s=64

**Status: paused for low expected value as of 2026-08-29.** PacketMul gives a
clean conditional law for each nonzero packet, but it does not remove the
shared-profile correlation across 256 regions. Its implementation cost is
also unmeasured and may consume the packet route's approximately 1 ms gross
advantage. The proof and performance receipts are preserved for possible
reuse.

This candidate adds independent nonzero `GF(16)` packet multipliers to
`Riffle LDPCSplitState g=4 t=256 s=64`.  It is a distinct construction.

PacketMul does not improve the current worst occupation-128 profile.  The
all-quad input law, its 7,105.481-bit margin, and its termination behavior are
unchanged.  PacketMul instead gives an exact survivor-count factorization:
conditioned on being nonzero, every multiplied packet has the same uniform
nonzero-nibble law.

That factorization enabled the first complete scan of all 16,335 shared-group
profiles at occupation 128.  Every profile is positive; the all-quad profile
is worst.  The construction remains exploratory.  It has no implementation,
no fixed multiplier schedule certificate, and no complete all-occupation
distance proof.

The main receipt is
`receipts/occupation128_full_profiles_fixed_point.json`.

`ACTIVE_PACKET_COLLAPSE.md` records an attempted proof simplification that
removes `n1,n2,n3,n4` from the inner interface.  The reduction is exact for
one region.  A rigorous profile-free exponential envelope was also tested,
but it loses 19,140 bits relative to the exact occupation-128 profile sum and
does not certify the target.  Its failure is a late-start bounding loss, not
a counterexample to the construction.
