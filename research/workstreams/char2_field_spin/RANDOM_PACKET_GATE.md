# Random Toeplitz packets: exact local law, incomplete finite coverage

Fresh random rectangular Toeplitz maps give a simple exact packet transfer over GF16.
They remove the former local-envelope slack without requiring an input diagonal, a random state refresh, or a full random invertible packet matrix.
The tested rate-2/3 profiles still have unclosed middle classes.
A rate-5/8 backup has promising sampled bounds, but no whole-code certificate or competitive runtime claim.

## Construction and probability space

Let q=16 and use a systematic outer [b,k] code independently in each of L source rows.
Its parity map is an independent uniform (b-k)-by-k Toeplitz matrix.
Route the b outer coordinates through b independent uniform permutations of L positions, in one regional routing stage.
Assume t divides L.

The inner state is R in GF16^s, initially zero.
For every width-t packet, sample fresh independent uniform Toeplitz matrices A of shape t-by-s and C of shape s-by-t.
The primary version emits and updates

    Y = X + A R,       R' = R + C X.

It carries the state across all packet and region boundaries and discards the final state.
The setup samples all outer matrices, permutations, and packet matrices independently.
There is no input diagonal and no GL state refresh.

The map from all inputs to all outputs is invertible because it is causal with identity feedthrough.
Given earlier inputs, the current state is known and X=Y-A R recovers the current input.
The augmented map on (X,R) need not be invertible, and no Schur-complement condition is required.
The update uses C X, not C Y; replacing it by C Y changes the transfer below.

An optional version replaces X in the emission by B X.
It samples B independently from invertible matrices whose action sends every fixed nonzero X to a uniform nonzero vector.
Full GL sampling is unnecessary, but the primary B=I version avoids this operation entirely.

## Exact packet transfer

First fix the input word, outer matrices, route, and all earlier packet randomness.
The current X and R are then fixed; the current A and C remain independent and fresh.
For every fixed nonzero vector u, a uniform rectangular Toeplitz matrix T makes T u uniform over its output space.
To see this, fix the first nonzero coordinate j of u and condition on all but the output-indexed Toeplitz diagonals.
The remaining linear map is triangular with nonzero diagonal u_j.
This argument applies to every rectangular shape, including inputs whose first nonzero coordinate lies beyond the output dimension.

Write

    Z(z) = ((1+15z)/16)^t,       p = 16^-s.

Rows and columns of each following matrix index zero and nonzero state.
If X=0, the state stays zero or nonzero as it was; a nonzero state produces a uniform output through A.
Thus the exact output-weight moment matrix is

    T_0(z) = [[1, 0], [0, Z(z)]].

For X of positive weight w, C X is uniform and independent of the emission.
The exact matrix for B=I is

    T_w(z) = [[p z^w, (1-p) z^w],
              [p Z(z), (1-p) Z(z)]].

These laws hold for every fixed incoming state, not merely for a uniform incoming state.
Empty packets never kill a nonzero state.
For the transitive-B version, replace z^w in the first row by

    ((1+15z)^t-1)/(16^t-1).

The scripts also return exact Fraction matrices at rational z.
The later optimization and coefficient calculations remain floating-point screens.

All generator coefficients lie in GF16.
A proved GF16 distance would extend to GF(2^m) for m in {8,32,64,128} without a union over extension-field messages.
For any nonzero extension-field word, choose a GF16-linear functional that preserves one nonzero coordinate.
Applying it coordinatewise gives a nonzero base-field codeword with contained support.
No such complete distance proof is claimed here.

## Source counts and the middle-class gap

For a fixed nonzero source row u, the rectangular outer parity is uniform over GF16^(b-k).
Consequently, its averaged complete value measure is 16^-(b-k) times counting measure on u!=0 and arbitrary parity values.
No input diagonal is needed for this identity.
The sparse calculation preserves the condition that every selected row has a nonzero systematic part.

For a region, let R_h(z) be its two-state moment matrix given h prescribed nonzero inputs at uniformly permuted positions.
It is obtained from the coefficient of degree h in

    (sum_(w=0)^t binom(t,w) u^w T_w(z))^(L/t),

divided by binom(L,h).
The source dynamic program tracks how many selected rows have appeared in the systematic coordinates.
Its positive transition contains

    H_seen(new) = sum_old binom(seen,old) 15^(old+new) R_(old+new).

The new implementation uses H_(seen+1)(new)=H_seen(new)+H_seen(new+1).
This positive Pascal recurrence reduces construction of the source transition from cubic to quadratic work in the active-row count.
It agrees with the existing source program in independent tests.

For a selected active-row count a, the complete class bound also includes binom(L,a), the outer factor 16^(-(b-k)a), and one projective division by 15.
At z=1, the source moment equals (16^k-1)^a, which the tests check independently.

Two cheaper middle screens are retained: a global coefficient tilt and separate coefficient tilts for each regional boundary-state entry.
Both drop the nonempty-systematic condition on individual rows.
A marker enforcing only total systematic weight at least a did not improve the tested witnesses; the optimizer chose marker one.
The difference from the exact-source calculation is substantial and must not be hidden by reporting only dense classes.

## Bounded K20 results

All following values are optimized floating margins in bits, not directed certificates.
The packet width is sixteen and B=I.
Negative values mean that this upper bound does not prove the target; they do not establish construction failure.

| Outer | Rate | State | a=1 | a=4 | a=16 |
|---|---:|---:|---:|---:|---:|
| [48,32] | 2/3 | 2 | 14.282 | -41.422 | not checked |
| [48,32] | 2/3 | 4 | 30.227 | 91.844 | not checked |
| [64,44] | 11/16 | 4 | 38.664 | 122.322 | not checked |
| [64,48] | 3/4 | 4 | 22.576 | 58.824 | not checked |
| [64,40] | 5/8 | 3 | 49.650 | 127.808 | 35.617 |
| [64,40] | 5/8 | 4 | 54.698 | 185.772 | 209.913 |
| [64,40] | 5/8 | 5 | 58.582 | 207.509 | 339.181 |
| [96,64] | 2/3 | 3 | 65.854 | 157.978 | -94.685 |
| [96,64] | 2/3 | 4 | 71.909 | 245.790 | 167.714 |
| [96,64] | 2/3 | 5 | 75.876 | 278.909 | 362.493 |

Larger source classes distinguish the last two useful-rate possibilities:

| Outer, state | a=32 | a=64 | a=128 | a=256 |
|---|---:|---:|---:|---:|
| [64,40], s4 | 131.376 | 37.321 | not checked | not checked |
| [64,40], s5 | 283.356 | 216.203 | 208.435 | 478.169 |
| [96,64], s4 | not checked | -582.412 | not checked | not checked |
| [96,64], s5 | not checked | -313.141 | -989.492 | -1867.465 |

The [64,40], s5 backup has K=1048960, L=26224, and N=1678336 after padding.
For a=128, the exact-source witness z=0.9703716207340908 gives about 208.435 bits, whereas the regional relaxation gives about -417.792 bits.
For a=256, z=0.9413912472726532 gives about 478.169 bits.
The regional screen becomes positive at the tested a=320, with about 77.973 bits; it gives 379.386 at a=384 and 1155.235 at a=512.
These are isolated class evaluations, not interval coverage.

The explicit remaining gaps are therefore:

- For the preferred rate-2/3 [96,64] profile with s<=5, the exact-source upper bound still fails at the tested middle classes, including a=64,128,256 for s5.
- For [64,40], s5, the positive sampled sparse and middle bounds suggest a possible overlap near a=256..320, but most individual classes and all covering intervals remain unchecked.
- No directed rational replay or sum over every active-row class has been completed for any new profile.

A scalable nonempty-source envelope is the next proof step.
It should be compared with the exact-source values above, especially at a=128, before another kernel is considered.
The rate-5/8 backup merits further work only if a separate full-cost analysis supports a genuine advantage.

## Reproduction and scope

```sh
python -B -m unittest discover -s research/workstreams/char2_field_spin/theory -p random_packet_tests.py -v
python -B research/workstreams/char2_field_spin/theory/random_packet_screen.py --outer 64 40 --width 16 --state 5 --identity --active 1 4 16 32 64 128 256
python -B research/workstreams/char2_field_spin/theory/random_packet_screen.py --outer 64 40 --width 16 --state 5 --identity --regional --active 320 384 512
python -B research/workstreams/char2_field_spin/theory/random_packet_screen.py --outer 96 64 --width 16 --state 5 --identity --active 64 128 256
```

Nine tests cover rectangular Toeplitz projections, binary exhaustive packet laws, two-dimensional state, output normalization, source counting, the Pascal transform, and a recorded floating witness.
No production or native code, remote execution, benchmark, or commit was added.
