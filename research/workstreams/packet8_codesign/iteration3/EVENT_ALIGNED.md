# Fractional blocks aligned with potential events

Fixed-width blocks can contain only one nonempty physical step when few
outer groups are selected. Pairing consecutive nonempty steps is another
way to sum hidden state paths before taking a fractional power. The
construction, setup distribution, and physical step length remain unchanged.

The recurrence below is valid. The first numerical controls did not improve
the existing fixed-four-step bound. These evaluations use floating arithmetic;
they are not outward probability certificates.

## Objects and conditioning

Fix a subset of q outer groups and a weight variable z in (0,1]. A region
has E=64 physical steps, each containing W=8 potential byte slots. Let
J=(j_1,...,j_E) record the number of selected slots in each step. A nonempty
step means j_i>0, even if every selected byte later receives label zero.

Let T_j be the nonnegative comparison operator after potential-label thinning.
The implementation uses the potential-birth basis from
`potential_birth_gate.py`. Its common stochastic change of basis preserves
every unpowered conditional comparison moment. This comparison still fills
some impossible state destinations; it is not the exact physical kernel.

Conditional on J, the positions inside each physical step are uniform.
The operators T_j already average those positions and the uniform byte labels.
For 0<alpha<=1, entrywise powers are applied only after multiplying the
operators inside a block. In particular, route multiplicities are not powered.

## Unique pair decomposition

Scan a region from its first physical step. End each complete block at the
second nonempty step since the previous cut. Such a block has the form

```
T0^g0 Tj1 T0^g1 Tj2,       g0,g1>=0,  1<=j1,j2<=W.
```

Its physical length is L=g0+g1+2 and its potential count is k=j1+j2.
Each complete block ends at a nonempty step. Zeros after that step belong
to the next block, or to the final tail.

After all complete pairs, the tail contains zero or one nonempty step.
Its forms are T0^L or T0^g0 Tj T0^g1, where g0+g1=L-1. The empty tail
has length zero and operator I. This rule gives exactly one decomposition
of every J, including all-zero sequences and sequences ending at a pair.

Write P^alpha for the entrywise power of a nonnegative matrix P, with
0^alpha=0. Define the unnormalized pair family

```
A[L,k] = sum_{g0+g1=L-2} sum_{j1+j2=k, j1,j2>=1}
             C(W,j1) C(W,j2)
             (T0^g0 Tj1 T0^g1 Tj2)^alpha.
```

The terminal family is

```
B[L,0] = (T0^L)^alpha,
B[L,k] = C(W,k) sum_{g0+g1=L-1}
                    (T0^g0 Tk T0^g1)^alpha,     1<=k<=W, L>=1.
```

Set B[0,0]=I. All undefined or infeasible entries are zero matrices.
Neither family contains a placement denominator.

## Ordered recurrence and probability bound

Let D[e,q] sum concatenations of complete pairs consuming e steps and q
potential slots. Initialize D[0,0]=I and every other entry to zero. Append
the last complete pair on the right:

```
D[e,q] = sum_{L=2}^e sum_{k=2}^{2W} D[e-L,q-k] A[L,k].
```

In particular, D[e,0]=0 for e>0. Adding zero-only prefixes here would count
the same sequence more than once. The full regional operator is

```
R_q = [sum_{L=0}^E sum_{k=0}^W D[E-L,q-k] B[L,k]] / C(WE,q).
```

The single denominator is the number of q-subsets in one region. The pair
and tail coefficients recover exactly the probability of each fine-count
sequence: product_i C(W,j_i)/C(WE,q).

For fixed fine counts across all regions, expand the original comparison
moment only at the selected block boundaries. Nonnegative path weights and
the inequality (sum_i x_i)^alpha<=sum_i x_i^alpha give the fractional
majorant represented by the block products above. Averaging fine counts
afterward gives R_q. Pairing starts anew in each region, but the state does
not reset: retain the complete matrix and evaluate e0 R_q^32 1.

Conditional clipping therefore gives the occupancy-q bound

```
C(512,q) * (beta^q * z^(-cutoff))^alpha * e0 R_q^32 1.
```

The outer group-subset union C(512,q) is distinct from the regional placement
denominator. It remains outside alpha. The cutoff and beta factor remain
inside alpha. This is the same conditional-probability argument used by
the fixed-block fractional gate; only the deterministic partition of J changes.

## Alpha-one identity

The identity at alpha=1 also has a noncommuting generating-function proof.
Use t to count steps and x to count potential slots. Define formal matrix
series

```
G = (I-t T0)^(-1),
V = t sum_{j=1}^W C(W,j) x^j Tj.
```

Complete pairs have series (GV)^2. The tail has series (I+GV)G. Thus

```
(I-(GV)^2)^(-1) (I+GV) G
    = (I-GV)^(-1) G
    = (I-t T0-V)^(-1).
```

The first equality uses powers of the same matrix series GV; it does not
commute G with V. The final series is ordinary chronological placement.
Taking coefficient t^E x^q and dividing by C(WE,q) proves the regional
alpha-one regression for every entry, not only the initial-state moment.

## Cost and numerical representation

At E=64 and W=8, the local block construction has 129,024 pair terms and
16,640 one-event tail terms. It also retains 65 all-zero tails. This is
much smaller than the 9^8 fine-count tuples used by a fixed eight-step block.
The subsequent pair DP is additional work; it is not included in these counts.

Long zero runs need special care. In the present comparison basis,
T0[0,0]=1 and each nonzero row has only a U-column entry v_r. Put u=T0[U,U].
For g>=1, T0^g[r,U]=v_r u^(g-1). These entries can underflow before raising
them to alpha, even when the powered values remain representable. Scaling
the entire matrix by its maximum does not resolve this: the 00 entry stays one.

`event_aligned.py` therefore represents all local products, block sums, and
DP coefficients as logarithms. Negative infinity denotes a structural zero.
Its matrix product first removes separate row and column scales. It uses
ordinary multiplication only if a conservative minimum-log check proves
that every positive scalar product is above the binary64 normal threshold.
Otherwise it uses the full log-semiring product. The logarithmic scales
are never exponentiated. This prevents positive matrix entries from being
silently replaced by zero; finite-precision rounding remains, so the output
is still only a proposal.

The generic interface is

```
event_regional(potential, degree, epochs=64, alpha=.4,
               check_alpha1=True, force_log_products=False, progress=None)
    -> (regional_logs[degree+1,dimension,dimension], metadata).
```

It accepts the actual 16-bit or saved 24-bit potential operators without
changing the geometry. The default alpha-one control shares the block
census and compares every regional entry against ordinary all-log placement.

## Bounded controls and next step

Six tiny tests cover independent noncommuting sequence parsing, empty and
one-event tails, scalar normalization, fractional path domination, long-run
underflow, and guarded-product agreement with the full log product.

The baseline 16-bit controls used the rebased potential family:

| q | Weight tilt | alpha | Event pairs | Rebased fixed four-step blocks |
|---|---:|---:|---:|---:|
| 16 | .06 | .4 | -117.280740418 bits | -111.731661406 bits |
| 64 | .25 | .4 | -266.470349392 bits | -203.932667140 bits |

The receipts are `event_pair_baseline16_q16_v1.json` and
`event_pair_baseline16_q64_v1.json`. Both retain source hashes, comparison
operators, the selected regional matrix, and the alpha-one regression.
The q16 evaluation took about 9.3 seconds; q64 took about 33.2 seconds.
These are local proof-evaluation times, not encoder benchmarks.

An independent caller reused the same helper with cached actual24 operators.
Its q16 margin was +14.715358493 bits at tilt .06 and alpha .4, versus
+17.354183 bits for fixed four-step blocks. That result is recorded in
`actual24_event_q16_v1.json`; it is not a whole-code certificate.

These partitions are not refinements of the fixed-four-step partition.
Consequently, fewer cuts in some sparse intervals do not imply a stronger
bound everywhere. The requested controls show no benefit at their tested
parameters. Keep the stronger fixed-block witnesses. Any further event
scheme should first explain which additional hidden paths it would combine;
a broad event-pair parameter sweep is not justified by these results alone.
