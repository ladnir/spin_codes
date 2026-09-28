# Retaining the central component in the dense cover

The two-coordinate cover can lose the relation between a large outer
count and its input distribution. Its three separate categorical maxima
may come from compositions different from the one maximizing the count.
This refinement retains one additional coordinate to restrict those choices.
It uses the same comparison measures and the same SPIN ensemble as
[MIXTURE_BOUND.md](MIXTURE_BOUND.md).

## The extra coordinate

Each pair component i combines two row components from the positive BCH
measure bound. Let c_i be half the number of those rows that use its
Bernoulli(1/2) component. Thus c_i is 0, 1/2, or 1. This is a label in
the comparison measure, not a classification of actual BCH row weights.

Retain the earlier tilted category features a_i1,a_i2 and define

    f_i = (a_i1,a_i2,c_i),
    m = sum_i (n_i/G) f_i,    G=4096.

The inactive component has f_i=0. Every composition with at least q_min
active labels projects into [0,1]^3. Covering the entire cube, while
discarding only certified empty cells, therefore covers every such
composition. The extra coordinate neither restricts the ensemble nor
drops any comparison component.

The component counts are integers. Consequently the central coordinate is
an integer multiple of 1/8192. Before bounding a cell, the verifier rounds
its central-coordinate endpoints inward to this grid. If no grid point
remains, the cell contains no composition. Also, let epsilon_j be the
smallest positive category feature f_ij. A positive mean in coordinate j
is at least epsilon_j/4096, because at least one component contributes.
If the cell's upper endpoint is smaller, its only possible mean is zero.
These contractions preserve all integer compositions. They do not change
the split tree or its coverage of the original cube.

## The same bound in three dimensions

Let C be a rectangular cell for m. Write b_i for component i's active
label and rho=q_min/G. For any values v_i, vector eta, and mu>=0, put

    h = max_i (v_i - eta dot f_i + mu*b_i).

Every composition in C with at least q_min active labels satisfies

    sum_i (n_i/G) v_i <= h + max_{m in C}(eta dot m) - mu*rho.

The componentwise inequalities defining h prove this after averaging;
the occupancy inequality supplies the final term. All terms can be
checked as exact rational numbers. Numerical linear programs only propose
eta and mu. In particular, taking every v_i=0 and obtaining a negative
right-hand side proves that C contains no feasible composition. A solver's
infeasibility status alone is not used to prune a cell.

For an alternative positive input tilt z', let Z'_i=pi_i dot z' and
D'_i=d_i*(Z'_i)^256, where d_i is the pair-component coefficient.
Apply the preceding affine bound three times, to v_ij=pi_ij/Z'_i.
This yields upper categorical weights w_j for the whole cell. Let S be
their sum and let K(lambda) be the two-state epoch envelope under iid
packet probabilities w/S.

For a separate exponential witness eta and mu>=0, define

    P = sum_i D'_i exp(eta dot f_i + mu*b_i),
    ell = min_{m in C}(eta dot m).

The multinomial theorem gives the complete-cell first-moment bound

    exp(lambda*209715) R_C^256 S^(4096*256)
      * e_zero K(lambda)^16384 (1,1)^T
      * P^4096 exp(-4096*ell - mu*q_min).

Here R_C is the three-category density bound, optionally tightened by
the category-count caps from the earlier method. The affine and
exponential witnesses need not be the same. Both bound every composition
in the same cell, which is the correlation the extra coordinate retains.

The old two-coordinate bound is also valid on any lifted cell. The driver
uses it when it suffices. Summing repeated upper bounds on different cells
can lose precision but cannot omit bad messages. The terminal sum excludes
no feasible central-component fraction.

## Verification and current scope

`lifted_mixture.py` first subdivides the packet-feature coordinates. When
their widths reach 1/16, it also subdivides the central fraction. The split
rule is deterministic and is used unchanged during replay. The verifier
checks the entire split tree, recomputes exact separation inequalities,
and evaluates accepted bounds with 192-bit Arb arithmetic. It does not
trust stored floating scores or stored upper bounds.

The exponential witness objective is smooth within each of its eight
sign regions. The proposal routine optimizes the needed regions separately,
using the correct endpoint of the cell in each coordinate. This avoids
choosing a witness for an infeasible cell center and then paying a large
radius penalty. Any proposed witness remains subject to outward replay.

A pilot covered the previously failing rectangle

    m1 in [12287/65536,3/16],
    m2 in [13909/65536,6955/32768],
    central fraction in [0,1],    q>=1024.

Splitting the central fraction into 16 intervals left 11 feasible cells;
each was outward-verified, with the smallest margin exceeding 8480 bits.
The other five cells were excluded by exact separating inequalities.
This is a complete cover of one selected rectangle, not the full dense
range. The full-domain run reached its 10,000-cell budget with 4,993
accepted or empty leaves and 15 pending subtrees. Its accepted-cell sum
has more than 70.3541 bits of margin, but does not cover the pending
subtrees. Most of the root cube is still pending; the number of checked
cells must not be treated as a percentage of proof completion.
The saved partial tree can be resumed:

```sh
python -B research/workstreams/permutation_locality/two_bit/lifted_mixture.py --minimum-groups 1024 --max-cells 10000 --max-depth 48 --max-unresolved 1 --output tmp/two-bit-lifted-dense-1024.json
python -B research/workstreams/permutation_locality/two_bit/lifted_mixture.py --replay tmp/two-bit-lifted-dense-1024.json
python -B -m unittest discover -s research/workstreams/permutation_locality/two_bit -p 'test_*.py'
```

Replay refuses incomplete covers. Full-code closure also needs every
occupancy below q_min and an outward aggregate below 2^-40.
If a search exhausts its cell budget, `--resume` continues its saved
unresolved partition. It reconstructs every cell from the split tree and
rechecks all accepted bounds and empty-cell claims first. Stored cell
coordinates, proposal scores, and aggregate bounds are not trusted.
The new cell budget counts only further search work:

```sh
python -B research/workstreams/permutation_locality/two_bit/lifted_mixture.py --resume tmp/two-bit-lifted-dense-1024.json --max-cells 20000 --max-depth 48 --output tmp/two-bit-lifted-dense-1024-continued.json
```

The test suite currently has 49 tests. New exhaustive checks compare
three-coordinate cell bounds with exact small composition sums, retain
feasible boundary points during separation, and reject incomplete
three-dimensional split trees. Resume tests also reject missing or
overlapping cells and false empty-cell claims. Integer-contraction tests
retain exact composition means, including threshold and grid boundaries.
Production and four-bit files are unchanged.

The subsequent [affine moment bound](AFFINE_MOMENT.md) keeps the inner
moment correlated with the outer count on each cell. It uses convexity
of the exact input polynomial and outward bounds at cell vertices, not
an assumed convexity property of the two-state envelope.
