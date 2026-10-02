# Mixed packet types: retain uncertainty about their original counts

The [homogeneous comparison](BOUNDARY_DIAGNOSIS.md) extends to arbitrary
mixtures of packet types. The new inequality retains randomness in the
posterior type counts after observing packet categories. It accepts
multiple stochastic laws and constant categories that overlap their
supports. The encoder and setup distribution do not change.

The first application is to completely enumerated boundary cells at
5% distance. This is not yet a full dense-domain cover. The complete
four-bit certificate remains [0.5% with 43.744 bits](FIRST_CLOSURE.md).

## Comparison within one region

Fix n slots and categorical laws f_i on packet weights 0,...,4. Let n_i
slots have law f_i, with sum_i n_i=n. Independently sample their
categories, then uniformly shuffle the slots. Write P for this
distribution. We work after the positive category tilt from the dense
proof. Its normalization is restored below.

Set p_i=n_i/n, omitting types with n_i=0. Define a joint reference
experiment: at each slot independently sample a type I with distribution
p, then sample a category Y with distribution f_I. Its category marginal
Q is iid with common probabilities x_j=sum_i p_i f_i(j).

Let E be the event that the reference type counts equal the fixed vector
(n_i). This event has positive probability

    e = n! / product_i n_i! * product_i p_i^n_i.

Conditioned on E, the reference type labels form a uniform permutation
of the prescribed multiset. Thus the category distribution conditional
on E is exactly P. Bayes' rule gives, for every y with Q(y)>0,

    P(y)/Q(y) = Pr[E | Y=y] / e.

The earlier type-conditioning bound could replace the numerator by one.
We instead bound it using the uncertainty remaining in the type labels.

## Equal observations give a multinomial posterior

For each category j, let m_j be the number of prescribed slots whose
law is deterministically category j. Every y in P's support contains
at least m_j occurrences of j. If m_j>0, then x_j>0. Conditioned on the
whole observed category sequence y, the type labels remain independent;
at a slot with category j their probabilities are

    r_i^(j) = p_i f_i(j) / x_j.

Fix the first m_j positions with observed category j. Their type-count
vector S has distribution Mult(m_j,r^(j)). Let T count the types at the
remaining positions. Conditional on y, S and T are independent. Hence

    Pr[E | Y=y] = sum_t Pr[T=t | Y=y] Pr[S=(n_i)-t | Y=y]
                 <= M_j,

where

    M_j = max { Pr[Mult(m_j,r^(j))=s] : sum_i s_i=m_j, 0<=s_i<=n_i }.

The coordinate caps hold because T has nonnegative coordinates. This
proves the pointwise bound

    P(y) <= B Q(y),       B = min(1, min_{j:m_j>0} M_j) / e.

If there is no deterministic category, use B=1/e. Sequences outside P's
support need no posterior bound. The result applies to every nonnegative
function of the sequence; it makes no assumption about the inner state.

Equal laws may be merged before constructing this reference experiment.
The shuffle cannot distinguish their labels, so the category distribution
P is unchanged. This avoids paying for artificial distinctions between
identically distributed types. Original outer counting coefficients and
location multiplicities are not merged or reduced.

Unlike the earlier homogeneous formula, B need not be the exact maximum
of P(y)/Q(y). It is a general upper bound computed with exact rationals.
The verifier also retains the previous category-count bound and takes
the smaller of the two valid losses.

## Computing M_j exactly

For fixed m and probabilities r, the multinomial mass is

    m! * product_i r_i^s_i / s_i!.

Increasing s_i from k-1 to k multiplies the product by r_i/k. These
marginal gains decrease with k. Select the largest m gains, with at most
n_i gains from coordinate i. Ties can be resolved arbitrarily.

The implementation starts with s_i=min(n_i,floor(m r_i)). Every selected
gain is at least 1/m, and every available unselected gain is at most
1/m. Greedy completion is therefore exact. Coordinates with r_i=0 receive
no count in a positive-mass outcome. For m=0 the mass is one. If the
positive-probability coordinates cannot hold m counts under the caps,
the maximum mass is zero.

Floating arithmetic proposes witnesses only. Replay recomputes the mode
and its mass with rationals, then evaluates the inner moment outward.
No numerical optimizer or approximate mode is trusted by the verifier.

## Combining this comparison with the SPIN bound

For each original packet law pi_i and positive category tilt t, use
Z_i=sum_j pi_i(j)t_j and f_i(j)=pi_i(j)t_j/Z_i. Apply the preceding
inequality to f_i. Its iid reference has mean x; undo the tilt with
unnormalized reference weights w_j=x_j/t_j.

The existing proof then uses the same outer coefficients, full group
multinomial multiplicity, product_i Z_i^(256 n_i), and reference mass
(sum_j w_j)^524288. The iid inner kernel receives w/sum_j w_j.
Independence of the 256 regional comparison measures gives loss B^256
for the whole input sequence. Inner setup remains independent and
unchanged. Uniform lane shuffles preserve the pointwise comparison on
four-bit packets conditional on their weight categories.

## Numerical checks at 5%

All 14 selected compositions pass with the posterior comparison. Some
new mixed examples are listed below. Labels have the meanings from
[the boundary diagnosis](BOUNDARY_DIAGNOSIS.md#what-was-tested); omitted
groups are `0000`. The entries are log2 upper bounds on contributions
to the expected bad-message count at cutoff 104857, not whole-code
failure margins.

| Composition | Log2 upper |
|---|---:|
| 64 `0003`, one `0002` | -351.57 |
| 32 `0003`, 33 `0002` | -859.57 |
| 30 `0003`, 30 `0002`, five `0034` | -1143.49 |
| 32 `0003`, 33 `0034` | -3304.39 |
| 32 `0003`, 32 `0034`, one `1111` | -3443.53 |
| 64 `2222`, one `1111` | -1823.58 |

For comparison, 65 `0003` gives -415.66 under this general inequality,
versus -420.88 under the exact homogeneous formula. Thus the extension
loses little on that example while accepting the mixtures that the
homogeneous formula rejects.

## Complete cells, not interpolation between points

`atlas.py --posterior-shuffle` tries the new bound on completely
enumerated cells with at most 32 feasible compositions. It can also
retain the exact homogeneous witness when that is better. During replay,
the verifier regenerates the complete composition list and checks that
every listed term occurs exactly once. Omissions and duplicates fail.
An enumeration work limit returns no conclusion, never a partial list.

`boundary_cells.py` exercises this path on nested four-dimensional boxes
near the least-active component, with q>=65. For the base category tilt
(1,1/32,1/64,1/128,1/256), define a=65/(49*2048) and epsilon=2^-32.
Each tested box has

    x_1 in [a-epsilon, a+excess/2048],
    x_2, x_3, x_4 in [0,epsilon].

The small positive widths exclude other component types by exact
enumeration, not by assuming their absence. These boxes have nonzero
width in all four coordinates. They overlap and must not be counted as
disjoint pieces of a cover.

| Excess | Complete composition count | Aggregate log2 upper |
|---|---:|---:|
| 1/64 | 2 | -351.57 |
| 1/32 | 7 | -351.57 |
| 1/16 | 24 | -347.05 |

The largest box uses 20 posterior-comparison witnesses and four exact
homogeneous witnesses. The aggregate includes every composition's full
outer coefficient and location multiplicity. It is a bound for this
entire box, not an extrapolation from its corners.
All 14 selected-composition bounds and all three complete-cell bounds
were independently replayed at 256-bit precision. All 26 regression
tests pass.

The existing 268 pending regions of the partial 1% atlas all exceed the
current bounded-enumeration capability at q>=65. The new inequality has
not certified those broad regions. The next step is either a cover that
isolates manageable composition lists or a posterior bound uniform over
all compositions in a broader cell. Mean packet weights alone do not
specify the posterior type-count loss.

## Reproduction

Run from the repository root with the existing Python dependencies:

```sh
python -B -m unittest discover -s research/workstreams/permutation_locality/independent_rows/dense_closure -p 'test_*.py'
python -B research/workstreams/permutation_locality/independent_rows/dense_closure/boundary_probe.py --distances 1/20 --posterior-shuffle --output tmp/four-bit-boundary-posterior-shuffle.json
python -B research/workstreams/permutation_locality/independent_rows/dense_closure/boundary_cells.py --output tmp/four-bit-boundary-cells-posterior.json
python -B research/workstreams/permutation_locality/independent_rows/dense_closure/boundary_probe.py --replay tmp/four-bit-boundary-posterior-shuffle.json --precision 256
python -B research/workstreams/permutation_locality/independent_rows/dense_closure/boundary_cells.py --replay tmp/four-bit-boundary-cells-posterior.json --precision 256
```

The tests exhaust small capped multinomial modes and small heterogeneous
category-count distributions, including overlapping constant categories.
They also verify complete-cell replay and rejection of missing terms.
Tests and local-cell bounds do not replace a complete dense-domain cover.

The 192-bit selected-composition and local-cell records have SHA-256
identifiers, respectively:

    edb987a674bfeb0de6419b24ac5b9c3c4bf4a465286d6fa14bb47eba168854b4
    45fe68b7104f4912d9b84057afe4aac2299798561b9ad9e6231e02236ccd462f

Raw records remain under `tmp/`. The existing complete 0.5% witness has
not been modified.
