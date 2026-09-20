# Riffle FullLaneFieldConv g=4 t=s=256

This folder records a distinct inner construction obtained by strengthening
the four-long-accumulator checkpoint.  The 256-bit state occupies all four
lanes, and a fresh nonzero multiplication in `GF(2^256)` mixes the complete
epoch output into the next state.

The strongest current diagnostic reaches relative distance 0.09 with 55.95
bits of regular-word margin under the modeled `[256,128,38]` outer spectrum.
It covers every maximally packed rank profile `4q+r`.  The dense profile has
about 119,987 bits of margin.  These numbers are not yet a certificate: the
main missing lemma must show that rank packing bounds every group-rank
histogram of the same total rank.

The full lane state matters.  Two and three state dimensions per four-bit
column improve the live-output exponent but still fail at high outer
multiplicity.  Four dimensions remove lane-mask dependence and make a live
epoch a uniform 256-bit affine coset.

See [CONSTRUCTION.md](CONSTRUCTION.md), [PROOF_STATUS.md](PROOF_STATUS.md),
and `receipts/profile_independent_k4_delta09.json`.
