# Riffle TransposePacketShuffle-RandomStepConv

This folder records the transpose-first construction that preserves fixed
`g`-bit packets. It is the active member of the transpose-first comparison.

The evaluator uses a packed active-group upper model. Packing domination is
proved in `proof/PACKING_DOMINATION.md`. The complete \(g=4,B=1024,\sigma=15\)
occupation sum is certified in `proof/DISTANCE_CERTIFICATE.md`.

`TRADEOFF_G4.md` records three certified 40-bit trade points down to
\(B=256,\sigma=18\). `TRADEOFF_G8.md` records the distinct bulk-packing
obstruction at width eight and six certified boundary points.
