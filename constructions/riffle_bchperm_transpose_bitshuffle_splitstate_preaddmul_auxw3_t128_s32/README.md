# Riffle BCHPerm-TransposeBitShuffle-SplitState-PreAddMul-AuxW3 t=128 s=32

**Status: fixed-inner regular proof closed.** This variant keeps the state,
the fixed map `B`, and the two accumulator layers of its parent. It replaces
the 32 weight-seven auxiliary columns of `A` with fixed weight-three columns.

Complete enumeration of the selected image of `A` gives minimum distance 24.
The map satisfies `BA=0`, all 128 coordinate forms are nonzero and distinct,
and no three coordinate forms are dependent. It therefore satisfies the
distance-20 profile used by the parent's 9% regular proof. The inherited
regular margin is 80.8139 bits.

Under the existing logical accounting, the change removes 128 sparse-map XOR
contributions per 128-bit epoch. The inner proxy decreases from 9.3657 to
8.3657 XORs per output bit. This is a logical estimate; the variant has not
yet received a wall-clock implementation.

The modeled outer spectrum, large mixed all-one tail, and outward-rounded
calculation remain open exactly as in the parent.

See `PROOF_STATUS.md` and the linked receipts for the exact scope.
