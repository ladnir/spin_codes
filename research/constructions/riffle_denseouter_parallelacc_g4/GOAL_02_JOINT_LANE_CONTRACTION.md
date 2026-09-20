# Goal 02: joint-lane contraction

## Question

Can the four accumulator lanes be analyzed jointly without classifying packet
values or adding randomization?

## Target

At binary length \(2{,}097{,}408\) and distance threshold \(40\):

1. derive a low-output moment bound that uses all four state-bit weights;
2. make the bound uniform over every packet multiset and active-value order;
3. compute a constant \(\rho<\sqrt2-1\) such that a fixed outer word of bit
   weight \(w\) fails with probability at most \(\rho^w\); and
4. insert \(\rho\) into the repaired dense-outer generating function.

The result must retain the same construction. Numerical optimization may
choose proof parameters separately for each fixed profile because those
parameters do not change the encoder.

## Pass condition

Goal 02 passes if it reduces the outer memory required for either 20 or 40
failure bits under the finite first-moment bound.
