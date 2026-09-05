# Goal 01: packetized dense-outer warmup

## Question

Does replacing the uniform bit permutation and scalar accumulator by a uniform
four-bit packet permutation and four parallel accumulators preserve a provable
positive relative distance?

## Target

Prove both of the following statements.

1. For every fixed outer word of binary weight \(w\), bound its low-output
   probability by \(z^w\). The bound must be uniform over all packet values.
2. Compose that bound with the dense outer generating function to obtain
   linear minimum distance with high probability.

The proof may lose constants. It must not assume a favorable packet-value
profile or enumerate prefix-XOR orbits.

## Pass condition

Goal 01 passes if it gives an explicit constant \(\delta_4>0\) such that

\[
\Pr[d_{\min}\le \delta L]=o(1)
\]

for every fixed \(\delta<\delta_4\) and a stated logarithmic outer-memory
constant.
