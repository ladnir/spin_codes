# Riffle FieldPair512-P4-ScalarStepConv g=4 sigma=16

This entropy-optimization candidate uses the same independently sampled
FieldPair512 data maps and global four-bit packet permutation as the P2
candidate. It uses four field-parity symbols and a one-lap scalar random
convolution with 16 retained state bits.

The evaluated central support interval has first-moment exponent
`-44.905837`. This is evidence for the target `2^-40`, not yet a certificate:
the remaining support tails and floating-point rounding still need bounds.

Attempts to reduce the state to 15 bits by using five or six parity symbols
fail in the current model. Their sampled exponents are respectively
`+87.981540` and `+79.371762`. In both cases, a one-active-data-block family
is already sufficient to create the obstruction.

See `REPORT.md` and `manifest.json` for the entropy ledger and claim scope.
`OUTER_INNER_CURVE.md` records the first constituent-size versus state curve.
