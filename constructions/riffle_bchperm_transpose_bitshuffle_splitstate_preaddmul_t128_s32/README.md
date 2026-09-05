# Riffle BCHPerm-TransposeBitShuffle-SplitState-PreAddMul t=128 s=32

**Status: fixed-inner regular proof closed.** This variant moves the existing
state multiplication before the syndrome addition. The change adds no
logical operation. It removes the structured termination output that blocks
the parent SplitState proof.

At 9% relative distance, every regular outer occupation closes under the
recorded constituent profiles. The aggregate regular margin is 80.8139 bits.
The proof uses a fixed-code shell average; it does not assume random affine
cosets.

The inner maps are now fixed. Exact enumeration proves that the selected
map `A` has distance 28 and satisfies the coordinate conditions. Exact dual
enumeration followed by the MacWilliams transform gives the complete kernel
spectrum of the selected map `B`. Its nonactivation probability is below the
profile used by the proof for every nonzero input weight.

The remaining substantive tasks are to replace the modeled outer spectrum
with an explicit outer code and to cover mixed cases with at least two
all-one outer words and more than 64 total active blocks. The numerical
certificate also needs outward-rounded arithmetic.

See `CONSTRUCTION.md` and `PROOF_STATUS.md` for the exact claim.
