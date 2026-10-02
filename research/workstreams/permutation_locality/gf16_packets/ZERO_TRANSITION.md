# Count the Weighted Zero Transition Exactly

The regional diagnostics exposed a loss before nonzero-state mixing
matters. The comparison operator bounded a zero-to-zero transition by
combining an unweighted feedback probability with an output-weight
bound. Repeating that loose entry can dominate a long product. Adding
transvections cannot improve this entry: every transvection fixes zero.

The GF character calculation already used for birth classes also gives
the weighted zero transition. This refinement changes the proof, not
the construction.

## Identity

Fix the feedback map C with s output bits and W four-bit input packets.
Conditional on j active packets, their positions are a uniform j-subset
of the W positions. Each active value is independently uniform among
the 15 nonzero four-bit values. Call the resulting input X. For
0<z<=1, define

    q_j(z) = E[z^wt(X) 1{CX=0} | j active packets].

An entering zero state outputs X and advances to CX. Therefore q_j(z)
is exactly its weighted zero-to-zero transition, for any update count.

For a feedback character u in F_2^s, let r_w(u) be the number of ones
in the four-bit vector u^T C_w. Character orthogonality gives

    q_j(z) = [v^j] sum_u product_w (1 + v a_(r_w(u))(z))
             / (2^s binom(W,j)),

    a_r(z) = ((1+z)^(4-r) (1-z)^r - 1) / 15.

Here [v^j] means the coefficient of v^j. The product chooses j packet
positions; a_r averages their weighted character over nonzero labels.
The sum over characters imposes CX=0. The indicator of the zero
feedback value has Walsh transform identically one, so character-profile
multiplicities provide its census row. No expansion-state maximization
or independence approximation is needed.

## Implementation and Scope

`occupancy_birth_classes.class_masses(..., include_zero=True)` prepends
this census row to the existing nonzero birth classes. Arb evaluates
the polynomial coefficients and returns outward upper endpoints.
`refine_zero` takes the smaller of this bound and the old zero entry,
leaving every other entry unchanged. Both bounds concern the same
scalar transition, so this minimum preserves state domination.

The selected regional verifier enables it with `--exact-zero`, recorded
as `regional_exact_zero` in the witness. It regenerates the count at
the current output tilt and precision. Old witnesses retain their old
meaning. Small exhaustive tests check the zero entry, the sum over all
feedback classes, and multistep noncommuting transition products. The
tests include full-rank, rank-deficient, and zero feedback maps.

The initial zero-path diagnostics used the old upper entries, not the
exact q_j(z). Their positive scores therefore did not establish a floor
for the relaxed counting measure, much less an obstruction for actual
SPIN. They identified an entry worth sharpening. Complete-domain
coverage and a summed sparse/dense bound remain necessary for a new
distance claim.

At three updates and 9.9% distance, fresh 384-bit evaluation verifies
log2 upper bound -177.7997549331 on the mean interval
1/9 +/- 1/1000000, including every variance part and q>=49. This
corresponds to reference activity 0.40. Independent 256-bit evaluation
agrees to the displayed digits. The activity-0.45 diagnostic
remains positive (+2255.3303 at the best tested output tilt). The scalar
search now includes the exact zero count in newly proposed regional
witnesses; saved witnesses without the option are unchanged. Combining
the exact zero entry with the nonlinear shared-mass bound improves the
activity-0.45 score to +1897.8712, still positive. No encoding parameters
were changed for these checks.
