# Random-outer parameter landscape

This exploration measures the relationship among outer length `B`, packet
width `g`, and inner memory `sigma`. Each outer block is an independent
random rate-half linear code. The base model has no global parity blocks.

For fixed parameters, let `X(B,g,sigma)` count the nonzero messages whose
encoded output has relative weight below 9%. The reported margin is

```text
lambda(B,g,sigma) = -log2 E[X(B,g,sigma)].
```

where the expectation covers the random outer codes, the global packet
permutation, and the independent random inner transitions. The target contour
is `lambda=40`.

The current receipts cover the first target contours for `g=1,2,4` and a
full-support obstruction at `g=8`. The calculations are diagnostic. Boundary
cells require exact outer coefficients, a complete support cover, and outward
rounding before they become proof claims.

See `REPORT_G4_INITIAL.md` and `REPORT_PACKET_WIDTH.md` for the results.
