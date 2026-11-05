# A deterministic rank bound from regional counts

The full/touched-region bound loses information when a region is almost full.
Retaining every regional cardinality gives a stronger bound at essentially
no counting cost. It controls outer shortening, not the probability of a route.

## Statement and proof

Let the outer be a direct sum of M MDS [B,k] codes, with B=2k.
Each routed row places exactly one coordinate in every region. Fix any route
and any coordinate set T. Write a_j for the number of its coordinates in
region j, and T_i for its occupancy in outer row i. Its shortened dimension is

    D(T)=sum_(i=1)^M (T_i-k)_+.

For every set J of k regions, row i selects at most k coordinates inside J.
Its excess (T_i-k)_+ is therefore at most its selected count outside J.
Summing over rows and choosing the k regions with largest counts proves

    D(T) <= sum of the k smallest a_j,
    rank(P_T)=|T|-D(T) >= sum of the k largest a_j.          (1)

Here P is a parity-check matrix of the routed outer and P_T contains its
columns indexed by T. The ranks and dimensions are over the code's field.
Equation (1) holds for every packet route, without an independence assumption.

The bound is optimal among binary row/region incidence matrices with those
column counts. Select the first a_j rows in column j. After sorting columns
by decreasing a_j, each row's excess equals its selected count after column k.
This nested matrix attains (1). It need not be attainable by a particular
packet route: optimality is asserted only for the incidence-matrix relaxation.

Completing the k most populated regions and applying the earlier secant bound
also gives (1). No other completion can improve it, because every completion
must upper-bound the same attaining incidence matrix. Thus a completion search
is unnecessary. With M=4 and counts (3,3,1,1), the old cap is 4 and (1) gives 2.

## A separable counting interface

Let Top_k(a) denote the sum of the k largest counts. The integer identity

    Top_k(a)=min_(theta=0,...,M) [k theta+sum_j(a_j-theta)_+]

gives, for tau>=1,

    tau^(-Top_k(a))
      <= sum_(theta=0)^M tau^(-k theta)
                            product_j tau^(-(a_j-theta)_+).   (2)

Every term separates across regions. The sum in (2) is an upper bound,
not an equality: the corresponding maximum is equal to its left side.
No implementation or numerical contraction is asserted for this interface.

## Scope and checks

The [bad-route construction](BAD_ROUTE_OBSTRUCTION.md) proves that some allowed
routes have distance at most 1.171875%, whatever the inner coefficients.
Consequently (1) cannot replace route averaging in a 35% whole-code theorem.
Use it for concentrated or nearly full supports alongside the
[guarded random-route moment](GUARDED_REGIONAL_MOMENT.md).

Root's serial exact test of `theory/regional_shortening.py` passed 135424 binary
matrices, 7467 count classes, 14196 virtual completions, 37800 monotonicity checks,
and 32 large-integer fixtures. It also checks nested attainment and the
complement identity. These are local combinatorial tests, not a distance
certificate. Reproduce with `python -B theory/test_regional_shortening.py`
from this workstream, under a 30-second/256-MiB cap.
