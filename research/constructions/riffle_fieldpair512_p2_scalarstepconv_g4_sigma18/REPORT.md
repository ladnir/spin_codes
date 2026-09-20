# Initial entropy-optimization report

The dense random matrices were not entropy-efficient. Their proofs used only
the marginal image of each fixed nonzero vector.

FieldPair512 replaces every random \(512\times256\) outer matrix by two field
multipliers. ScalarStepConv replaces every random inner matrix by one field
multiplier. These replacements preserve the exact marginal laws and retain
independence between blocks or steps.

Reducing the state from 20 to 19 bits leaves a sampled pointwise maximum of
`-69.1951`. The crude support union gives `-50.195`, which clears the target.

Reducing the state to 18 bits requires exact outer coefficients. The dominant
exact point occurs near support 759. The three evaluated intervals sum to
`-43.5986`. This is the minimum state size supported by the current analysis.

At 17 state bits, the evaluated central interval is positive. A different
inner argument or additional outer distance would be required to use that
state size.

The sigma-18 conclusion remains diagnostic. Its 3.60-bit margin must absorb
the omitted support tails, numerical rounding, and the negligible correction
to uniform global parity. Goal 02 should turn these items into a certificate
before optimizing the permutation family.
