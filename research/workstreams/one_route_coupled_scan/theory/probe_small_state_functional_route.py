"""Bounded FLOAT complete census using one common positive functional.

Fix g=(c,1)^T with c>0, common to all regions. For a nonnegative row v,
t=c*v0/(v*g) belongs to[0,1]. For positive raw local M, set
A=(M*g)0/c=M00+M01/c and B=(M*g)1=c*M10+M11. Then
E[(v*M*g)^alpha]/(v*g)^alpha=f(t)=sum_ell pi_ell[t*A_ell+(1-t)*B_ell]^alpha.
This is concave. Bounded derivative bisection finds a witness t; its global
tangent, not a sampled maximum, supplies lambda>=sup_[0,1] f.

The unconditioned128-batch norm bound is lambda^128. Divide ONCE by a
lower bound on the full4096-pair histogram probability. A common c avoids
norm-switch costs between regions. For alternating cells i,j, the full
route bound is max(c,1)^alpha *[(lambda_i^128/H_i)(lambda_j^128/H_j)]^128.
This scalar formula retains all intermediate-state sums inside alpha.
Complete cell-pair/run coverage and exact source counts are enforced.

All4845 raw occupancy profiles are included. This is ordinary floating
arithmetic, NOT a directed certificate. No files are written. Root owns
every run; use Python -B and external limits.
"""

import argparse
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import sys

sys.dont_write_bytecode = True
import probe_small_state_census_anchor as census


dense, np, bridge = census.dense, census.np, census.bridge
OUTPUT_LIMIT = 256 * 1024


def prepare(z, alpha, log2c, guard):
    table = dense.raw_local_table(z, guard)
    if table.shape != (4845, 2, 2) or np.any(~np.isfinite(table)):
        raise ArithmeticError("Incomplete positive raw local table")
    first_log = np.logaddexp2(table[:, 0, 0], table[:, 0, 1]-float(log2c))
    second_log = np.logaddexp2(table[:, 1, 0]+float(log2c), table[:, 1, 1])
    shift = float(max(first_log.max(), second_log.max()))
    with np.errstate(over="raise", invalid="raise", under="raise"):
        first, second = np.exp2(first_log-shift), np.exp2(second_log-shift)
    if np.any(first <= 0) or np.any(second <= 0):
        raise ArithmeticError("A positive functional coefficient vanished")
    guard()
    return dict(first=first, second=second, difference=first-second, alpha=float(alpha),
                shift=float(alpha)*shift, linear=np.empty(4845), powered=np.empty(4845),
                work=np.empty(4845)), dict(occupancies=4845, log2c=str(log2c),
                functional="g=(2^log2c,1); common to every cell and region",
                raw_A_log2_range=[float(first_log.min()), float(first_log.max())],
                raw_B_log2_range=[float(second_log.min()), float(second_log.max())],
                powered_log2_shift=float(alpha)*shift)


def value_derivative(prepared, weights, point):
    linear, powered, work = (prepared[key] for key in ("linear", "powered", "work"))
    # Use two positive terms, not B+t*(A-B): cancellation at t=1 could
    # erase a valid small A when the raw row asymmetry is large.
    np.multiply(prepared["first"], point, out=linear)
    np.multiply(prepared["second"], 1-point, out=work)
    np.add(linear, work, out=linear)
    if np.any(linear <= 0):
        raise ArithmeticError("A positive affine moment vanished")
    np.power(linear, prepared["alpha"], out=powered)
    value = float(weights @ powered)
    np.divide(prepared["difference"], linear, out=work)
    np.multiply(work, powered, out=work)
    derivative = prepared["alpha"]*float(weights @ work)
    if not math.isfinite(value) or value <= 0 or not math.isfinite(derivative):
        raise ArithmeticError("Nonfinite functional value or derivative")
    return value, derivative


def regional_cell(prepared, groups, cell, tilt, pseudocount, steps, guard):
    guard()
    rational = dense.proposal_for_cell(groups, cell, tilt, pseudocount)
    probabilities = np.asarray(tuple(map(float, rational)))
    if (np.any(probabilities <= 0) or np.any(~np.isfinite(probabilities))
            or abs(float(probabilities.sum())-1.) > 1e-12):
        raise ArithmeticError("Rational proposal conversion lost positive mass")
    with np.errstate(over="raise", invalid="raise", under="raise"):
        weights = bridge.profile_weights(probabilities)
    if (np.any(weights <= 0) or np.any(~np.isfinite(weights))
            or abs(float(weights.sum())-1.) > 1e-10):
        raise ArithmeticError("Complete positive striped proposal did not normalize")
    low, high, iterations = 0., 1., 0
    value, derivative = value_derivative(prepared, weights, low)
    if derivative <= 0:
        point = low
    else:
        value, derivative = value_derivative(prepared, weights, high)
        if derivative >= 0:
            point = high
        else:
            for iteration in range(steps):
                if (iteration & 3) == 0:
                    guard()
                point = (low+high)/2
                value, derivative = value_derivative(prepared, weights, point)
                if derivative > 0:
                    low = point
                else:
                    high = point
                iterations += 1
            point = (low+high)/2
            value, derivative = value_derivative(prepared, weights, point)
    # A tangent to a concave f majorizes it on the entire interval. Take
    # the larger tangent endpoint; do not accept f(point) as the supremum.
    tangent = max(value-point*derivative, value+(1-point)*derivative)
    if not math.isfinite(tangent) or tangent <= 0 or tangent < value:
        raise ArithmeticError("Invalid global functional tangent upper")
    log_lambda = math.log2(tangent)+prepared["shift"]
    endpoints = [bridge.histogram_log_probability(bridge.pair_histogram(groups, full),
                                                 probabilities) for full in cell]
    divisor = min(endpoints)
    region = 128*log_lambda-divisor
    if not math.isfinite(region):
        raise ArithmeticError("Nonfinite complete regional norm bound")
    guard()
    return region, dict(full_pair_low=cell[0], full_pair_high=cell[1],
                        covered_counts=cell[1]-cell[0]+1, proposal=list(map(str, rational)),
                        tangent_point=str(Fraction(point)), bisection_iterations=iterations,
                        derivative_bracket=[low, high], normalized_value=value,
                        normalized_derivative=derivative, normalized_tangent_upper=tangent,
                        log2_lambda_upper=log_lambda,
                        endpoint_histogram_log2=endpoints, divisor_log2=divisor,
                        regional_log2=region)


def scalar_census(groups, cells, regional, alpha, log2c, guard, progress, details):
    """Complete compatible-pair maxima for the scalar functional operator."""
    guard()
    count, maximum_runs = len(cells), min(groups, bridge.GROUPS-groups)
    regional = np.asarray(regional)
    if regional.shape != (count,) or np.any(~np.isfinite(regional)):
        raise ArithmeticError("Incomplete regional functional table")
    operators = (128*(regional[:, None]+regional[None, :])
                 + float(alpha)*max(float(log2c), 0.))
    lower = np.asarray([cell[0] for cell in cells], dtype=int)
    upper = np.asarray([cell[1] for cell in cells], dtype=int)
    first = np.maximum(1, groups-upper[:, None]-upper[None, :])
    last = np.minimum(maximum_runs, groups-lower[:, None]-lower[None, :])
    candidates = np.flatnonzero(first <= last)
    order = candidates[np.argsort(operators.ravel()[candidates])[::-1]]
    successor = list(range(maximum_runs+2))
    maxima, maximizing = np.full(maximum_runs+1, -math.inf), [None]*(maximum_runs+1)

    def next_unfilled(runs):
        while successor[runs] != runs:
            successor[runs] = successor[successor[runs]]
            runs = successor[runs]
        return runs

    assigned = 0
    progress.update(computed_ordered_cell_pairs=count*count,
                    compatible_ordered_cell_pairs=len(candidates), assigned_run_classes=0)
    for index, flat in enumerate(order):
        if (index & 63) == 0:
            guard()
        even, odd = divmod(int(flat), count)
        runs = next_unfilled(int(first[even, odd]))
        while runs <= int(last[even, odd]):
            maxima[runs], maximizing[runs] = operators[even, odd], (even, odd)
            assigned += 1
            successor[runs] = next_unfilled(runs+1)
            runs = successor[runs]
        progress["assigned_run_classes"] = assigned
        if assigned == maximum_runs:
            break
    if assigned != maximum_runs or np.any(~np.isfinite(maxima[1:])):
        raise ArithmeticError("Incomplete source-run coverage")
    multiplicity, count_sum, terms = bridge.GROUPS, 0, []
    for runs in range(1, maximum_runs+1):
        if (runs & 63) == 1:
            guard()
        if multiplicity <= 0:
            raise ArithmeticError("Nonpositive exact cycle-run count")
        even, odd = maximizing[runs]
        terms.append(dict(runs=runs, source_sets=str(multiplicity),
                          log2_source_sets=math.log2(multiplicity),
                          log2_route_operator=float(maxima[runs]),
                          log2_counted_route_operator=math.log2(multiplicity)+float(maxima[runs]),
                          even_cell_index=even, odd_cell_index=odd,
                          even_full_pair_cell=list(cells[even]), odd_full_pair_cell=list(cells[odd])))
        count_sum += multiplicity
        progress["counted_run_classes"] = runs
        if runs < maximum_runs:
            numerator = multiplicity*(groups-runs)*(bridge.GROUPS-groups-runs)
            multiplicity, remainder = divmod(numerator, runs*(runs+1))
            if remainder:
                raise ArithmeticError("Nonintegral exact cycle-run count")
    expected = math.comb(bridge.GROUPS, groups)
    if count_sum != expected:
        raise ArithmeticError("Cycle-run counts did not partition all source sets")
    total = bridge.log2_sum(term["log2_counted_route_operator"] for term in terms)
    dominant = sorted(terms, key=lambda term: term["log2_counted_route_operator"], reverse=True)[:8]
    progress.update(exact_source_count_partition=True, complete_run_cover=True)
    result = dict(raw_anchor_census_operator_sum_log2=total,
                  raw_anchor_average_operator_log2=total-math.log2(expected),
                  exact_source_sets=str(expected), counted_source_sets=str(count_sum),
                  run_classes=maximum_runs, worst_run_term=dominant[0], top_run_terms=dominant)
    if details:
        result["run_terms"] = terms
    guard()
    return result


def probe(groups=1024, z=Fraction(1, 2), alpha=Fraction(3, 40), end_groups=None,
          max_cells=128, tail_counts=32, tilt=Fraction(6, 5),
          pseudocount=Fraction(1, 65536), log2c=Fraction(80), bisection_steps=32,
          seconds=45., memory_mib=1024, details=False):
    end_groups = groups if end_groups is None else end_groups
    z, alpha, tilt, pseudocount, log2c = map(Fraction, (z, alpha, tilt, pseudocount, log2c))
    if (not isinstance(groups, int) or not isinstance(end_groups, int)
            or not 1 <= groups < bridge.GROUPS or not groups <= end_groups <= bridge.GROUPS
            or not isinstance(max_cells, int) or not 1 <= max_cells <= 128
            or not isinstance(tail_counts, int) or not 0 <= tail_counts <= 32
            or not 0 < z < 1 or not 0 < alpha <= 1 or tilt <= 0 or pseudocount <= 0
            or not isinstance(bisection_steps, int) or not 20 <= bisection_steps <= 40
            or not math.isfinite(float(log2c)) or not 0 < float(z) < 1
            or not 0 < float(alpha) <= 1):
        raise ValueError("Require proper anchor, valid rational markers, <=128 bulk+32 tail cells,"
                         "finite log2c and20..40 bisection steps")
    budget = dense.Budget(seconds, memory_mib)
    cells = census.cells_for_group(groups, max_cells, tail_counts)
    regional, reports, result = [], [], {}
    progress = dict(computed_ordered_cell_pairs=0, assigned_run_classes=0,
                    counted_run_classes=0, complete_run_cover=False)
    status, error, complete, completed_sources = "complete-functional-route-float", None, False, 0
    try:
        prepared, receipt = prepare(z, alpha, log2c, budget.guard)
        result["local_kernel"] = receipt
        for index, cell in enumerate(cells):
            value, receipt = regional_cell(prepared, groups, cell, tilt, pseudocount,
                                           bisection_steps, budget.guard)
            regional.append(value)
            reports.append(dict(index=index, **receipt))
        result.update(scalar_census(groups, cells, regional, alpha, log2c,
                                    budget.guard, progress, details))
        source, receipt = dense.source_interval_log(groups, end_groups, z, alpha, budget.guard)
        completed_sources = end_groups-groups+1
        bound = result["raw_anchor_average_operator_log2"]+source
        if not math.isfinite(bound):
            raise ArithmeticError("Nonfinite complete source-interval bound")
        prefix = float(alpha)*(bridge.source_factor(groups, "joint")["log2_prefactor"]
                               - bridge.LIMIT*math.log2(float(z)))
        result.update(log2_uncapped_bound=bound, interval_source_factor_log2=source,
                      source_interval=receipt,
                      anchor_log2_uncapped_bound=result["raw_anchor_census_operator_sum_log2"]+prefix)
        for term in result["top_run_terms"]:
            term["anchor_log2_bound"] = term["log2_counted_route_operator"]+prefix
        selected = sorted({term[key] for term in result["top_run_terms"]
                           for key in ("even_cell_index", "odd_cell_index")})
        result["dominant_cell_reports"] = [reports[index] for index in selected]
        budget.guard()
        complete = True
    except (TimeoutError, MemoryError, ArithmeticError, RuntimeError) as failure:
        status = str(failure) if isinstance(failure, (TimeoutError, MemoryError)) else "arithmetic-error"
        error = str(failure)
    result.update(schema="small-state-functional-route-float-v1", status=status, error=error,
                  anchor_groups=groups, end_groups=end_groups, z=str(z), alpha=str(alpha),
                  tilt=str(tilt), pseudocount=str(pseudocount), log2c=str(log2c),
                  bisection_steps=bisection_steps, maximum_bulk_cells=max_cells,
                  requested_tail_counts=tail_counts, requested_cells=len(cells), completed_cells=len(regional),
                  unprocessed_cells=len(cells)-len(regional), cell_partition=[list(cell) for cell in cells],
                  requested_run_classes=min(groups, bridge.GROUPS-groups),
                  unprocessed_run_classes=min(groups, bridge.GROUPS-groups)-progress["counted_run_classes"],
                  requested_source_sizes=end_groups-groups+1, completed_source_sizes=completed_sources,
                  unprocessed_source_sizes=end_groups-groups+1-completed_sources,
                  complete_profile_cell_cover=complete, complete_source_interval_cover=complete,
                  census_engine=progress, source_mode="joint",
                  functional="(v*(c,1)^T)^alpha; one common c for all regions",
                  lambda_rule="global tangent upper of concave IID scalar moment",
                  regional_rule="lambda^128 / one full4096-pair histogram probability",
                  pair_rule="max(c,1)^alpha*(regional_even*regional_odd)^128",
                  source_accounting="complete averaged raw anchor and actual-r joint source factors",
                  interval_rule="exact full-product raw inclusion monotonicity; fixed markers and functional",
                  empty_input_rule="raw[[1,1],[256^-16,V_empty]]; omega=1",
                  terminal="weighted norm comparison; no physical final refresh or reset charge",
                  directed=False, whole_code_certificate=False, generated_file_bytes=0,
                  output_cap_bytes=OUTPUT_LIMIT, **budget.receipt())
    if details:
        result["cells"] = reports
    paths = {Path(module.__file__).resolve() for module in (census, dense, bridge)}
    paths.add(Path(__file__).resolve())
    provider = sys.modules.get("hd_common_refresh_packets")
    if provider is not None:
        paths.add(Path(provider.__file__).resolve())
    result["dependencies"] = {path.name: hashlib.sha256(path.read_bytes()).hexdigest()
                              for path in sorted(paths)}
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--groups", type=int, default=1024)
    parser.add_argument("--z", type=Fraction, default=Fraction(1, 2))
    parser.add_argument("--alpha", type=Fraction, default=Fraction(3, 40))
    parser.add_argument("--end-groups", type=int)
    parser.add_argument("--max-cells", type=int, default=128)
    parser.add_argument("--tail-counts", type=int, default=32)
    parser.add_argument("--tilt", type=Fraction, default=Fraction(6, 5))
    parser.add_argument("--pseudocount", type=Fraction, default=Fraction(1, 65536))
    parser.add_argument("--log2c", type=Fraction, default=Fraction(80))
    parser.add_argument("--bisection-steps", type=int, default=32)
    parser.add_argument("--seconds", type=float, default=45.)
    parser.add_argument("--memory-mib", type=int, default=1024)
    parser.add_argument("--details", action="store_true")
    result = probe(**vars(parser.parse_args()))
    encoded = json.dumps(result, allow_nan=False, sort_keys=True)
    if len(encoded.encode("utf-8")) > OUTPUT_LIMIT:
        result.pop("cells", None)
        result.pop("run_terms", None)
        result["details_omitted_for_output_cap"] = True
        encoded = json.dumps(result, allow_nan=False, sort_keys=True)
    if len(encoded.encode("utf-8")) > OUTPUT_LIMIT:
        raise RuntimeError("Compact functional receipt exceeded its output cap")
    print(encoded, flush=True)
    if not result["complete_profile_cell_cover"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
