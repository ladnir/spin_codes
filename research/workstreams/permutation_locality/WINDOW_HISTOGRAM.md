# Exact output moments from expansion-window histograms

The multi-window bound previously retained only the total expansion weight.
It bounded an average over distinct windows by a product bound and a chord
bound. Both relaxations lose information about the actual expansion map.
`window_histogram.py` computes that average directly. It changes proof
coefficients, not the encoder or its setup distribution.

## The finite calculation

Fix an entering nonzero state a in F_2^19. Split its expansion Ea into 32
four-bit windows. Let h_r count windows of weight r, for r=0,...,4, and let
v=sum_r r h_r. Exhaustive enumeration of the production expansion map gives
20 histograms across all 2^19 states, including the zero state. The code
checks their multiplicities against the existing expansion spectrum.

Fix a shape n=(n_1,n_2,n_3,n_4), with j=sum_b n_b<=32. Choose j distinct
windows uniformly, assign the weights with these multiplicities uniformly,
and choose independent uniform lane masks of the assigned weights. Let X
be the resulting input. This is the local input distribution under the
ideal grouped route, conditional on its active weights.

For 0<z<=1, define

    f_{r,b}(z) = sum_l binomial(r,l) binomial(4-r,b-l)
                      z^(b-2l) / binomial(4,b).

The sum ranges over max(0,r+b-4)<=l<=min(r,b). In a weight-r window,
z^r f_{r,b}(z) is the expected value of z raised to the output weight
after adding a uniform weight-b lane mask. Windows are disjoint, so

    E[z^wt(Ea+X)] = z^v [x_1^n_1 ... x_4^n_4]
                   product_{r=0}^4 (1+sum_{b=1}^4 f_{r,b}(z)x_b)^h_r
                   / (binomial(32,j) j!/product_b n_b!).

Here the bracket extracts the coefficient of the indicated monomial.
The denominator counts the possible assignments of weights to distinct
windows. Each coefficient term has exactly j lane factors. This formula
uses sampling without replacement; it introduces no independence between
the selected window locations.

For z=exp(-lambda), each lane factor is rounded upward to a multiple of
2^-48 using Arb. Polynomial coefficients are then computed with unbounded
integers. A degree-j coefficient is divided by 2^(48j) and the assignment
count above. Arb evaluates the remaining z^v factor outward. Positivity
of every polynomial coefficient justifies these upper bounds.

## Use in the state bound

Taking a maximum over nonzero-state histograms bounds the moment for an
arbitrary entering state. Averaging with the exact histogram multiplicities
within an expansion-weight class bounds the corresponding uniform source.
For a fresh source, the script separately enumerates the feedback of every
one-window input and averages with its histogram distribution. It takes
the maximum over the four possible preceding window weights.

The two-update kernel is lazy with probability alpha=1/4 and refreshes
uniformly over nonzero states with probability beta=3/4. Write m=2^19-1.
A moment bound M contributes at most alpha M to lazy outgoing mass and
beta M N_w/m to the uniform envelope for an expansion class of size N_w.
The arbitrary-state refresh contribution to zero is at most beta M/m.
These estimates use that output is emitted before the two updates.

The existing pointwise-density and lazy-cancellation bounds remain in
place. In particular, an average output moment alone does not bound the
event that a lazy state cancels the feedback. The optional refinement is
applied after converting the envelope to the requested update count and
before the full-feedback refinement reuses the improved refresh bounds.

All-one-column and input-weight penalties multiply each shape's coefficient
before taking maxima. Entrywise minima with the old coefficients preserve
the same state-envelope inequalities. The extra mature-tail coordinates
remain nonnegative upper bounds; they are neither new probability mass
nor new assumptions about the source distribution.

## Checks and reproduction

The coefficient formula passes 235 exact rational checks against distinct
window assignments and explicit four-bit masks on small examples. Fifteen
independent production-map checks enumerate actual masks in one or two
windows, then compare high-precision moments with the outward result.
The histogram census also checks the production maps and all four fresh
feedback distributions. These tests support the implementation; the
coefficient identity above gives the general justification.

    python -B research/workstreams/permutation_locality/window_histogram.py --maximum 8

The local calculation and direct-mask checks also pass at 384-bit precision:

    python -B research/workstreams/permutation_locality/window_histogram.py --maximum 8 --precision 384

The optional `--window-histogram 8` flag computes all 495 shapes through
degree eight. Higher local occupancies retain their earlier bounds.
At tilt .032 and all-one penalty .75, it tightens 252 local coefficients
of the existing two-update model. This count is not a distance certificate.

    python -B research/workstreams/permutation_locality/occupancy_cdf_cover.py --groups 64 --mixing-rounds 2 --shortening-moments --positive-shortening --dual-shortening --dual-sandwich --window-histogram 8 --fresh --window-average --multi-average --prefix-flags --fresh-collision --zero-moment --full-feedback 6 --mature-tail 64 --pair-tail --joint-witness --penalties .75 --tilts .032 --probe-supports 128 144 160

This command is a selected-point binary64 diagnostic. It does not provide
outward full-support coverage or a full-code failure bound. The optional
flag leaves the default certificate path, production code, and paper unchanged.

The resulting log2 upper bounds are:

| Equal support in all 64 groups | Previous bound | With histogram moments |
|---|---:|---:|
| 128 | +2839.547802 | +2824.859997 |
| 144 | +2408.779214 | +2391.926914 |
| 160 | +2259.108135 | +2240.114367 |

The improvement is only 14.69--18.99 bits. Thus the loss in this output
moment bound is not the principal source of the remaining gap at these
points. The new calculation leaves lazy feedback cancellation unchanged.
The next useful target is a joint bound on cancellation and output weight,
rather than extending the histogram calculation to higher local occupancy
solely for these points.

## Fourteen-group outward replay

With the optional refinement, the complete fourteen-group cover still has
190 leaves after 25 splits. Its 192-bit outward sum is below
1.402421e-46, giving 152.3207736818 bits of margin for output weight at most
209715. The earlier bound was 152.2843977144 bits. The small improvement
does not extend the occupancy range.

    python -B research/workstreams/permutation_locality/occupancy_cdf_cover.py --groups 14 --mixing-rounds 2 --max-splits 100 --window-histogram 8 --fresh --window-average --multi-average --prefix-flags --fresh-collision --zero-moment --mature-tail 64 --pair-tail --retain-parents --joint-witness --penalties .5 .75 1 --tilts .0032 .005 .008 .01

This covers every support, rank, and group location at occupancy fourteen,
not every occupancy. No full two-update certificate follows.
