#!/usr/bin/env python3
"""Degree-type MIP relaxation for the fixed three-band connected tail.

For a connected colored partition pattern on ``s`` active blocks, block ``v``
has tile-cluster degrees ``(d0,d1,d2)``.  Pair capacity one implies
``d0+d1+d2 <= s+2``.  The count of degree-d vertices in a band is a multiple
of d, since those vertices are partitioned into d-clusters.  A d-cluster and
an e-cluster in two distinct bands intersect in at most one vertex.

The MIP maximizes the sum of the solver-free BCH cell-cap factor bounds over
integer degree-type counts subject only to these necessary conditions.  It
then applies the connected embedding bound and multiplies by Bell(s)^3,
counting even incompatible partition triples.  Thus the combinatorial MIP is
a relaxation in the safe direction.  By default its floating objectives and
solver make it a diagnostic.  ``--exact-moment-costs`` instead obtains exact
rational diagonal moments, rounds every rooted moment upward on a dyadic
grid, and rounds each log objective upward to an integer.  The exported LP can
then be solved by exact SCIP and independently checked with VIPR.

The optional degree grid rounds every cluster degree upward before evaluating
the objective.  This is safe because the rooted diagonal factor is an
``L_d`` norm: writing the OR kernel as an expectation over a Bernoulli mask
gives ``K_d(w)^(1/d) = (E p_M(w)^d)^(1/d)``, which is nondecreasing in ``d``.
The combinatorial constraints continue to use the original unrounded degree.
"""

from __future__ import annotations

import argparse
import itertools
import math
from decimal import Decimal, ROUND_CEILING, localcontext
from fractions import Fraction
from pathlib import Path

import numpy as np

try:
    import highspy
except ImportError as error:
    raise SystemExit(
        "tail degree MIP requires highspy in a temporary tool environment"
    ) from error

from probe_bch_three_band_cell_cap_motifs import (
    BAND_SIZES,
    cell_caps,
    factor_table,
    greedy_factor_bound,
)
from probe_bch_three_band_transport_lp import (
    exact_rows,
    remove_trivial_codewords,
    triples,
)
from outward_log2 import log2_fraction
from three_band_exact_moments import (
    factor_table_upper,
    greedy_factor_bound_upper,
)


BLOCKS = 16384
LANES = 64
GRAPH_CODIMENSION = 24


def bell(number: int) -> int:
    values = [1]
    for size in range(1, number + 1):
        values.append(
            sum(math.comb(size - 1, prefix) * values[prefix] for prefix in range(size))
        )
    return values[number]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--blocks", type=int, default=7)
    parser.add_argument("--tilt", type=float, default=0.075)
    parser.add_argument("--target-pole", type=float, default=0.181)
    parser.add_argument("--support-cutoff", type=int, default=106)
    parser.add_argument("--solver-output", action="store_true")
    parser.add_argument("--relax-integrality", action="store_true")
    parser.add_argument(
        "--pairwise-majorant",
        action="store_true",
        help="replace every three-degree factor by its pairwise half-sum upper envelope",
    )
    parser.add_argument(
        "--cauchy-majorant",
        action="store_true",
        help="use the rigorous best pair-versus-single Cauchy upper envelope",
    )
    parser.add_argument(
        "--joint-selector-lp",
        action="store_true",
        help="add pairwise joint selector hulls for a stronger LP relaxation",
    )
    parser.add_argument(
        "--integer-cost-scale",
        type=int,
        default=0,
        help="round every log2 objective coefficient upward on this integer grid",
    )
    parser.add_argument(
        "--cost-guard-bits",
        type=float,
        default=1e-7,
        help="conservative per-coefficient guard before integer rounding",
    )
    parser.add_argument(
        "--exact-moment-costs",
        action="store_true",
        help=(
            "derive integer objective costs from exact rational moments and "
            "outward dyadic roots (requires --integer-cost-scale)"
        ),
    )
    parser.add_argument(
        "--root-bits",
        type=int,
        default=160,
        help="dyadic precision for --exact-moment-costs",
    )
    parser.add_argument(
        "--degree-bins",
        default=None,
        help=(
            "comma-separated degrees for monotone upward objective rounding, "
            "for example 1,2,3,4,6,8,12,16,24,32,48,64"
        ),
    )
    parser.add_argument(
        "--write-model-prefix",
        type=Path,
        default=None,
        help="write one LP model per required band pair using this prefix",
    )
    args = parser.parse_args()
    if not (
        3 <= args.blocks <= 64
        and 0 < args.tilt <= args.target_pole < 1
        and args.support_cutoff >= 0
    ):
        raise SystemExit("tail degree MIP: invalid arguments")
    if args.exact_moment_costs and not args.integer_cost_scale:
        raise SystemExit("--exact-moment-costs requires --integer-cost-scale")
    if args.degree_bins is not None and (
        args.pairwise_majorant or args.cauchy_majorant
    ):
        raise SystemExit("degree bins are only implemented for direct factors")

    if args.degree_bins is None:
        degree_bins = tuple(range(1, args.blocks + 1))
    else:
        degree_bins = tuple(
            sorted({int(value) for value in args.degree_bins.split(",")})
        )
        if not degree_bins or degree_bins[0] != 1 or degree_bins[-1] < args.blocks:
            raise SystemExit("degree bins must start at 1 and reach the block count")
        if any(not 1 <= value <= 64 for value in degree_bins):
            raise SystemExit("degree bins must lie in 1..64")

    def rounded_degree(degree: int) -> int:
        return next(value for value in degree_bins if value >= degree)

    states = triples()
    exact = remove_trivial_codewords(states, exact_rows(states))
    caps, diagonal_mass, by_total = cell_caps(states, exact)
    factor_cache = {
        (band, degree): factor_table(BAND_SIZES[band], degree, args.tilt)
        for band in range(3)
        for degree in degree_bins
    }
    value_cache: dict[tuple[int, int, int], float] = {}
    array_value_cache: dict[tuple[tuple[int, int], ...], float] = {}

    def factor_log2(kind: tuple[int, int, int]) -> float:
        rounded_kind = tuple(rounded_degree(degree) for degree in kind)
        if rounded_kind not in value_cache:
            value_cache[rounded_kind] = math.log2(
                greedy_factor_bound(
                    states,
                    caps,
                    diagonal_mass,
                    by_total,
                    tuple(
                        factor_cache[band, rounded_kind[band]] for band in range(3)
                    ),
                )
            )
        return value_cache[rounded_kind]

    base_factor = factor_log2((1, 1, 1))

    def array_factor_log2(specification: tuple[tuple[int, int], ...]) -> float:
        """Specification entries are (cluster degree, integer power)."""
        if specification not in array_value_cache:
            arrays = tuple(
                factor_cache[band, degree] ** power
                for band, (degree, power) in enumerate(specification)
            )
            array_value_cache[specification] = math.log2(
                greedy_factor_bound(
                    states, caps, diagonal_mass, by_total, arrays
                )
            )
        return array_value_cache[specification]

    def objective_log2(kind: tuple[int, int, int]) -> float:
        if not args.pairwise_majorant:
            if not args.cauchy_majorant:
                return factor_log2(kind)
            candidates = []
            for singleton_band in range(3):
                pair_specification = tuple(
                    (kind[band], 0 if band == singleton_band else 2)
                    for band in range(3)
                )
                singleton_specification = tuple(
                    (kind[band], 2 if band == singleton_band else 0)
                    for band in range(3)
                )
                candidates.append(
                    0.5
                    * (
                        array_factor_log2(pair_specification)
                        + array_factor_log2(singleton_specification)
                    )
                )
            return min(candidates)
        first, second, third = kind
        return 0.5 * (
            factor_log2((first, second, 1))
            + factor_log2((first, 1, third))
            + factor_log2((1, second, third))
            - base_factor
        )

    size = args.blocks
    kinds = [
        kind
        for kind in itertools.product(range(1, size + 1), repeat=3)
        if sum(kind) <= size + 2
    ]
    kind_columns = len(kinds)
    cluster_keys = [
        (band, degree)
        for band in range(3)
        for degree in range(1, size + 1)
    ]
    cluster_column = {
        key: kind_columns + index for index, key in enumerate(cluster_keys)
    }
    selector_keys = [
        (band, degree, count)
        for band in range(3)
        for degree in range(2, size + 1)
        for count in range(size // degree + 1)
    ]
    selector_column = {
        key: kind_columns + len(cluster_keys) + index
        for index, key in enumerate(selector_keys)
    }
    joint_selector_keys = []
    if args.joint_selector_lp:
        joint_selector_keys = [
            (
                first_band,
                first_degree,
                second_band,
                second_degree,
                first_count,
                second_count,
            )
            for first_band, second_band in itertools.combinations(range(3), 2)
            for first_degree in range(2, size + 1)
            for second_degree in range(2, size + 1)
            for first_count in range(size // first_degree + 1)
            for second_count in range(size // second_degree + 1)
        ]
    joint_selector_column = {
        key: kind_columns + len(cluster_keys) + len(selector_keys) + index
        for index, key in enumerate(joint_selector_keys)
    }
    columns = (
        kind_columns
        + len(cluster_keys)
        + len(selector_keys)
        + len(joint_selector_keys)
    )
    raw_costs = [objective_log2(kind) for kind in kinds]
    if args.exact_moment_costs:
        pole_exact = Fraction(str(args.tilt))
        exact_factor_cache = {
            (band, degree): factor_table_upper(
                BAND_SIZES[band],
                degree,
                pole_exact,
                root_bits=args.root_bits,
            )
            for band in range(3)
            for degree in degree_bins
        }
        exact_value_cache: dict[tuple[int, int, int], Fraction] = {}

        def exact_factor(kind: tuple[int, int, int]) -> Fraction:
            rounded_kind = tuple(rounded_degree(degree) for degree in kind)
            if rounded_kind not in exact_value_cache:
                exact_value_cache[rounded_kind] = greedy_factor_bound_upper(
                    states,
                    caps,
                    diagonal_mass,
                    by_total,
                    tuple(
                        exact_factor_cache[band, rounded_kind[band]]
                        for band in range(3)
                    ),
                )
            return exact_value_cache[rounded_kind]

        objective_costs = []
        for kind in kinds:
            log_upper = log2_fraction(exact_factor(kind)).hi
            with localcontext() as context:
                context.prec = 100
                context.rounding = ROUND_CEILING
                objective_costs.append(
                    int(
                        (log_upper * Decimal(args.integer_cost_scale)).to_integral_value(
                            rounding=ROUND_CEILING
                        )
                    )
                )
    elif args.integer_cost_scale:
        if args.integer_cost_scale <= 0:
            raise SystemExit("tail degree MIP: integer cost scale must be positive")
        objective_costs = [
            math.ceil((cost + args.cost_guard_bits) * args.integer_cost_scale)
            for cost in raw_costs
        ]
    else:
        objective_costs = raw_costs
    costs = np.array(
        objective_costs
        + [0.0]
        * (len(cluster_keys) + len(selector_keys) + len(joint_selector_keys))
    )

    best = -math.inf
    best_required = None
    best_types = None
    best_solver_info = None
    for required_bands in itertools.combinations(range(3), 2):
        model = highspy.Highs()
        model.setOptionValue("output_flag", args.solver_output)
        column_upper = np.full(columns, size, dtype=float)
        column_upper[kind_columns + len(cluster_keys) :] = 1
        model.addCols(
            columns,
            costs,
            np.zeros(columns),
            column_upper,
            0,
            np.zeros(columns + 1, dtype=np.int32),
            np.zeros(0, dtype=np.int32),
            np.zeros(0),
        )
        if not args.relax_integrality:
            model.changeColsIntegrality(
                columns,
                np.arange(columns, dtype=np.int32),
                np.array([highspy.HighsVarType.kInteger] * columns),
            )

        lower: list[float] = []
        upper: list[float] = []
        starts = [0]
        indices: list[int] = []
        values: list[float] = []

        def add_row(items, row_lower, row_upper):
            indices.extend(index for index, _value in items)
            values.extend(value for _index, value in items)
            starts.append(len(indices))
            lower.append(row_lower)
            upper.append(row_upper)

        add_row([(index, 1) for index in range(kind_columns)], size, size)
        for band, degree in cluster_keys:
            add_row(
                [
                    (index, 1)
                    for index, kind in enumerate(kinds)
                    if kind[band] == degree
                ]
                + [(cluster_column[band, degree], -degree)],
                0,
                0,
            )
            if degree >= 2:
                add_row(
                    [
                        (selector_column[band, degree, count], 1)
                        for count in range(size // degree + 1)
                    ],
                    1,
                    1,
                )
                add_row(
                    [(cluster_column[band, degree], 1)]
                    + [
                        (selector_column[band, degree, count], -count)
                        for count in range(size // degree + 1)
                    ],
                    0,
                    0,
                )

        collision_degree = [sum(degree - 1 for degree in kind) for kind in kinds]
        add_row(
            [(index, collision_degree[index]) for index in range(kind_columns)],
            2 * (size - 1),
            size * (size - 1),
        )
        for band in required_bands:
            add_row(
                [(index, kind[band] - 1) for index, kind in enumerate(kinds)],
                2,
                highspy.kHighsInf,
            )

        for first_band, second_band in itertools.combinations(range(3), 2):
            for first_degree in range(2, size + 1):
                for second_degree in range(2, size + 1):
                    intersection = [
                        (index, 1)
                        for index, kind in enumerate(kinds)
                        if kind[first_band] == first_degree
                        and kind[second_band] == second_degree
                    ]
                    if not intersection:
                        continue
                    if args.joint_selector_lp:
                        first_max = size // first_degree
                        second_max = size // second_degree
                        for first_count in range(first_max + 1):
                            add_row(
                                [
                                    (
                                        joint_selector_column[
                                            first_band,
                                            first_degree,
                                            second_band,
                                            second_degree,
                                            first_count,
                                            second_count,
                                        ],
                                        1,
                                    )
                                    for second_count in range(second_max + 1)
                                ]
                                + [
                                    (
                                        selector_column[
                                            first_band, first_degree, first_count
                                        ],
                                        -1,
                                    )
                                ],
                                0,
                                0,
                            )
                        for second_count in range(second_max + 1):
                            add_row(
                                [
                                    (
                                        joint_selector_column[
                                            first_band,
                                            first_degree,
                                            second_band,
                                            second_degree,
                                            first_count,
                                            second_count,
                                        ],
                                        1,
                                    )
                                    for first_count in range(first_max + 1)
                                ]
                                + [
                                    (
                                        selector_column[
                                            second_band,
                                            second_degree,
                                            second_count,
                                        ],
                                        -1,
                                    )
                                ],
                                0,
                                0,
                            )
                        add_row(
                            intersection
                            + [
                                (
                                    joint_selector_column[
                                        first_band,
                                        first_degree,
                                        second_band,
                                        second_degree,
                                        first_count,
                                        second_count,
                                    ],
                                    -first_count * second_count,
                                )
                                for first_count in range(first_max + 1)
                                for second_count in range(second_max + 1)
                            ],
                            -highspy.kHighsInf,
                            0,
                        )
                        continue
                    # If the selected first-band cluster count is u, exact
                    # pair capacity gives z <= u*y_second.  The one-hot
                    # selector makes this implication linear.  M is the
                    # maximum possible product and disables unselected rows.
                    first_max = size // first_degree
                    second_max = size // second_degree
                    big_m = first_max * second_max
                    for first_count in range(first_max + 1):
                        add_row(
                            intersection
                            + [
                                (
                                    cluster_column[second_band, second_degree],
                                    -first_count,
                                ),
                                (
                                    selector_column[
                                        first_band, first_degree, first_count
                                    ],
                                    big_m,
                                ),
                            ],
                            -highspy.kHighsInf,
                            big_m,
                        )

        model.addRows(
            len(lower),
            np.asarray(lower),
            np.asarray(upper),
            len(indices),
            np.asarray(starts, dtype=np.int32),
            np.asarray(indices, dtype=np.int32),
            np.asarray(values),
        )
        model.changeObjectiveSense(highspy.ObjSense.kMaximize)
        if args.write_model_prefix is not None:
            args.write_model_prefix.parent.mkdir(parents=True, exist_ok=True)
            model_path = args.write_model_prefix.with_name(
                args.write_model_prefix.name
                + f"_{required_bands[0]}{required_bands[1]}.lp"
            )
            model.writeModel(str(model_path))
        model.run()
        if model.getModelStatus() != highspy.HighsModelStatus.kOptimal:
            raise SystemExit(f"tail degree MIP: {model.getModelStatus()}")
        objective = model.getObjectiveValue()
        if objective > best:
            best = objective
            best_required = required_bands
            solution = model.getSolution().col_value
            best_types = [
                (round(solution[index]), kinds[index])
                for index in range(kind_columns)
                if solution[index] > 0.5
            ]
            info = model.getInfo()
            best_solver_info = (
                info.mip_node_count,
                info.mip_dual_bound,
                info.mip_gap,
            )

    if args.integer_cost_scale:
        best /= args.integer_cost_scale

    embedding = (
        math.log2(BLOCKS)
        + (size - 1) * math.log2(LANES)
        - math.lgamma(size + 1) / math.log(2)
        - GRAPH_CODIMENSION
    )
    cutoff = args.support_cutoff * math.log2(args.target_pole / args.tilt)
    partition_overcount = 3 * math.log2(bell(size))
    total = embedding + cutoff + partition_overcount + best
    print("three-band connected-tail degree relaxation")
    print(
        f"blocks={size} tilt={args.tilt:.12f} "
        f"target_pole={args.target_pole:.12f}"
    )
    print(f"degree_types={len(kinds)} bell_cubed_log2={partition_overcount:.12f}")
    print(f"objective_degree_bins={degree_bins}")
    print(f"best_required_bands={best_required} best_types={best_types}")
    print(
        f"solver_nodes={best_solver_info[0]} "
        f"solver_dual_bound={best_solver_info[1]:.12f} "
        f"solver_gap={best_solver_info[2]:.12g}"
    )
    print(f"factor_objective_log2={best:.12f}")
    if args.integer_cost_scale:
        print(
            f"integer_cost_scale={args.integer_cost_scale} "
            f"cost_guard_bits={args.cost_guard_bits:.12g}"
        )
        if args.exact_moment_costs:
            print(
                f"moment_costs=EXACT_RATIONAL_OUTWARD_DYADIC "
                f"root_bits={args.root_bits}"
            )
    print(f"connected_tail_log2={total:.12f}")
    if args.exact_moment_costs and not args.relax_integrality:
        print("status=SAFE_COMBINATORIAL_RELAXATION_OUTWARD_INTEGER_MIP")
    else:
        print(
            "status=SAFE_COMBINATORIAL_RELAXATION_FLOATING_"
            + ("LP" if args.relax_integrality else "MIP")
        )


if __name__ == "__main__":
    main()
