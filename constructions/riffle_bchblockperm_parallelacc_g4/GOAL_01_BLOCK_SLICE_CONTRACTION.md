# Goal 01: one-block slice contraction

## Question

Does a full bit permutation inside each BCH block remove the worst-case
four-bit grouping loss before the global packet permutation?

## Target

For one fixed BCH word of binary weight \(h\):

1. derive its exact packet-support law after the 128-bit permutation;
2. average the joint four-lane gap moment over the complete weight-\(h\) slice;
3. compare the resulting bound with the old one-lane projection and a full
   global bit permutation; and
4. identify the obstruction at the intended relative-distance scale.

## Pass condition

Goal 01 passes if the block-averaged joint-lane calculation gives a valid
improvement for the minimum BCH weight \(h=22\). A full-code distance theorem
is outside this goal.

