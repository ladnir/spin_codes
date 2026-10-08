# Lower outer dimensions at the K20 target

GRS[16,6] with the existing four-state GF16 inner is a viable fallback in
numerical screens. GRS[16,7] with five states also has useful margin.
GRS[16,7] with four states remains short of the 40-bit target after summing
the sparse classes. These are numerical screens, not complete certificates.

Each candidate retains one routing stage: sixteen independent uniform
regional permutations. The inner retains width sixteen, independent nonzero
GF16 input scalars, uniform GL(s,16) refreshes, monomial A, and C=A^T.
Changing the outer dimension changes the code family.

## Exact sizes and the requested subcode

Let d denote the outer dimension and let L be the number of outer rows.
Each region must contain a multiple of sixteen symbols. For a requested
input dimension K_req, the screen sets

    L = 16 ceil(K_req/(16d)),    K_full = dL,    N = 16L.

For K_req=1048576, the resulting sizes are:

| Outer | L | K_full | Fixed zero inputs | N | Minimum distance target |
|---|---:|---:|---:|---:|---:|
| [16,8] | 131072 | 1048576 | 0 | 2097152 | 209716 |
| [16,7] | 149808 | 1048656 | 80 | 2396928 | 239693 |
| [16,6] | 174768 | 1048608 | 32 | 2796288 | 279629 |

The distance target is ceil(N/10). The numerical bounds concern the full
K_full-dimensional code. A future full-code certificate would also cover
the K_req-dimensional subcode obtained by fixing the indicated inputs to zero.
Restricting the message space cannot decrease minimum distance. The output
length stays N, so this restriction is distinct from puncturing outputs.

## Source bound and sparse results

For a prescribed nonempty outer support of weight w, its multiplicity is

    n_w = sum_(j=0)^w (-1)^j binom(w,j) 16^max(0,w-j-(16-d)).

It vanishes below weight 17-d. For rho>0, the screen uses

    eta_d(rho) = max_(17-d<=w<=16) n_w/(16^(w-(16-d))*rho^w).

This gives n_w <= eta_d(rho) rho^w 16^(w-(16-d)). The coupled regional
bound therefore divides by 16^((16-d)a) for a active rows. The previous
half-rate exponent 8a would be incorrect for either new dimension.
The source enumeration through four active rows remains exact. Each class
receives one projective division by 15 after choosing its active-row set.

The following values report minus log2 of a sum of numerical class bounds:

| Outer / states | Active-row classes included | Source tilt | Sum margin |
|---|---|---|---:|
| [16,7] / 4 | 1--128 | Optimize rho through class 16; then rho=0.96 | 37.8339 bits |
| [16,7] / 5 | 1--47 | rho=0.871 | 48.9835 bits |
| [16,6] / 4 | 1--47 | rho=0.871 | 58.7796 bits |

The first row includes every sparse class through 128. Its weakest class
is a=27, with margin 41.3884 bits; summing nearby classes loses the target.
The calculation does not establish a lower bound on code-failure probability.
It also does not exclude a sharper source bound for the same construction.

For [16,6]/state4, the coefficient-Chernoff screen has 94.3147 bits at a=48
and 170.7893 bits at a=64. A separate dense-tail screen begins at a=75735
with 40.4530 bits. These points suggest that a complete certificate is
feasible, but no interval cover or directed whole-code replay was produced.

## Cost tradeoff and next step

The [16,7] output has about 14.29% more symbols than the half-rate K20 output.
The [16,6] output has about 33.34% more. Five states also increase the inner
arithmetic and state-matrix storage. Four states can reuse the existing inner
kernel. These observations do not determine which full encoder is faster.
No benchmark was run for this screen.

Keep these as fallback candidates while assessing the separate half-rate
outer redesign. If a fallback is needed, complete the [16,6]/state4 directed
certificate before implementing or claiming its distance guarantee.

## Reproduction and checks

Run from `research/workstreams/char2_field_spin/theory`:

```sh
python -B outer_rate_screen.py --summary --dimension 7 --state 4 --active 128 --source-rho 24/25 --optimize-through 16
python -B outer_rate_screen.py --summary --dimension 7 --state 5 --active 47
python -B outer_rate_screen.py --summary --dimension 6 --state 4 --active 47
python -B outer_rate_screen.py --middle --dimension 6 --state 4 --active 48 64
python -B -m unittest test_outer_rate_screen -v
```

Five tests check natural geometry, the complete MDS enumerator, the exact
source envelope, directed source enumeration against a rational oracle,
and agreement with the preserved half-rate bounds. All five pass.
`--replay` checks individual selected classes with directed arithmetic;
it always labels them incomplete. Existing certificates and production
code remain separate from this new screen.
