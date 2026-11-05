"""Bounded FLOAT marker proposals from an already-built rational pair table.

fit_markers(actual_degrees, anchor_actual, guard, progress=None) returns
31 bank rows (reference_center, reference_lo, reference_hi, u, v) and a
compact proposal report.  Each u,v is a positive Fraction with denominator
at most2^20.  The caller supplies all2145 actual-occupancy H/J coefficients,
already built coherently at its chosen z,alpha.  This helper does not build
local caps, recompute alpha, change source factors, or evaluate a bound.

The complete degree triangle is relabeled once by(f,s)->(64-f-s,s).
No caps are evaluated on complementary inputs; all six channels and their
temporal order are retained.  For hole anchor k0=8192-anchor_actual, a
reference center x targets a=k0*x/2048 and s=k0-2a.  Only the two fitting
endpoints use half-slot proxies.  The frozen finite-chain objective has
numerator e0 H^16383 J and target vector256*(a,s).  Its marker-independent
source/count/normalization constants are set to0 for proposal selection.

Every fit uses the frozen at-most12-evaluation coordinate procedure.
Arithmetic failure retains a positive rational seed and never omits a
center.  The caller must independently replay the returned rational bank
with directed arithmetic and each actual target's own denominator.  FLOAT
logs, gradients, and fit gains are not certificate inputs.  No job, provider
import, file write, or parameter sweep is performed at module import.
"""

from collections.abc import Mapping
from fractions import Fraction as F
import hashlib
import importlib.util
import math
from pathlib import Path
import sys


HOLE_PIN = "a430c6c72bdc4c5b19ccf3d99da49fbc0e1385ee6361bc020469f8bc1806b217"
DENOMINATOR_LIMIT = 1 << 20
REFERENCE_WEIGHT = 2048
MAX_PROPOSALS = 12
BASE_CENTERS = (0,3,11,25,44,69,98,131,169,211,256,304,354,406,458,512,
                566,618,670,720,768,813,855,893,926,955,980,999,1013,1021,1024)


def _load(name, pin, label):
    path = Path(__file__).with_name(name)
    if hashlib.sha256(path.read_bytes()).hexdigest() != pin:
        raise RuntimeError("frozen fitter dependency mismatch: "+name)
    spec = importlib.util.spec_from_file_location(label, path)
    module = importlib.util.module_from_spec(spec)
    saved = sys.dont_write_bytecode
    sys.dont_write_bytecode = True
    try:
        spec.loader.exec_module(module)
    finally:
        sys.dont_write_bytecode = saved
    return module


def _dependencies(guard):
    """Load only pinned interfaces; none of their numerical mains is called."""
    guard()
    hole = _load("cheap_lane_ordinary_hole_directed.py", HOLE_PIN, "fit_hole_bank")
    screen = _load("cheap_lane_ordered_hole_screen.py", hole.HOLE_PIN, "fit_hole_interface")
    large = _load("cheap_lane_ordered_large_screen.py", screen.LARGE_PIN, "fit_large_interface")
    rho = _load("cheap_lane_rho_bucket_control.py", large.RHO_PIN, "fit_rho_interface")
    control = _load("cheap_lane_pair_profile_calibration.py", rho.CONTROL_PIN, "fit_pair_objective")
    moments = _load("cheap_lane_broad_profile_calibration.py", control.CALIBRATION_PIN, "fit_pair_moments")
    sizes = _load("cheap_lane_pair_sizes_screen.py", large.SIZES_PIN, "fit_marker_seeds")
    if control.MAX_PROPOSALS != MAX_PROPOSALS:
        raise ArithmeticError("frozen proposal budget changed")
    bank = tuple(hole.BANK)
    if tuple(row[0] for row in bank) != BASE_CENTERS:
        raise ArithmeticError("reference center layout changed")
    if [f for _,lo,hi,_,_ in bank for f in range(lo,hi+1)] != list(range(1025)):
        raise ArithmeticError("reference integer cells do not partition0..1024")
    guard()
    return bank, control, moments, sizes, {
        "hole_bank": HOLE_PIN, "hole_screen": hole.HOLE_PIN,
        "large_interface": screen.LARGE_PIN, "rho_interface": large.RHO_PIN,
        "objective": rho.CONTROL_PIN, "moments": control.CALIBRATION_PIN,
        "marker_seed": large.SIZES_PIN}


def _fraction(value):
    if isinstance(value, bool) or not isinstance(value, (int, F)):
        raise TypeError("local coefficients must be exact int/Fraction values")
    return F(value)


def _prepare(table, guard, capacity=64):
    """Validate and copy a read-only complete six-channel rational triangle."""
    if not isinstance(table, Mapping):
        raise TypeError("local coefficients must be a read-only-compatible Mapping")
    expected = {(f,s) for f in range(capacity+1) for s in range(capacity+1-f)}
    if set(table) != expected:
        raise ValueError("complete local degree triangle required before reversal")
    result = {}
    fingerprint = hashlib.sha256()
    positives = [0]*6
    zeros = 0
    for f,s in sorted(expected):
        guard()
        values = table[f,s]
        if len(values) != 6:
            raise ValueError("each pair degree must have six MM/MT channels")
        row = tuple(_fraction(value) for value in values)
        if any(value < 0 for value in row):
            raise ValueError("negative local coefficient")
        result[f,s] = row
        fingerprint.update((str(f)+","+str(s)+":"+",".join(map(str,row))+"\n").encode("ascii"))
        for channel,value in enumerate(row):
            positives[channel] += int(value > 0)
            zeros += int(value == 0)
    if not all(positives):
        raise ValueError("each channel needs a positive coefficient")
    return result, fingerprint.hexdigest(), positives, zeros


def _reverse(table, capacity):
    return {(capacity-f-s,s): row for (f,s),row in table.items()}


def _log_fraction(value):
    """Proposal-only log; never convert a tiny positive Fraction to float."""
    if value < 0:
        raise ValueError("negative logarithm argument")
    if not value:
        return -math.inf
    result = math.log(value.numerator)-math.log(value.denominator)
    if not math.isfinite(result):
        raise ArithmeticError("nonfinite rational coefficient log")
    return result


def fit_target(holes, center):
    if type(holes) is not int or not 641 <= holes <= 4095:
        raise ValueError("hole anchor must be an integer in641..4095")
    if type(center) is not int or center not in BASE_CENTERS:
        raise ValueError("unknown normalized reference center")
    full = F(holes*center, REFERENCE_WEIGHT)
    if center == 0:
        full = F(1,2)
    elif center == 1024:
        full = F(holes-1,2)
    single = holes-2*full
    if full <= 0 or single <= 0:
        raise ArithmeticError("fitting proxy must be strictly interior")
    return full, single


def _positive_rational(value):
    if not isinstance(value, F) or value <= 0:
        raise ValueError("positive Fraction proposal required")
    rounded = value.limit_denominator(DENOMINATOR_LIMIT)
    return max(F(1,DENOMINATOR_LIMIT), rounded)


def fit_markers(rational_degrees, anchor_actual, guard, progress=None):
    """Return (bank, FLOAT-only report), with no local rebuild or source term.

    A MappingProxyType is accepted; the supplied mapping and rows are never
    modified.  The caller records z/alpha and owns deadline, memory, disk
    preflight and directed replay.  Timeout/MemoryError always propagate.
    """
    if type(anchor_actual) is not int or not 4097 <= anchor_actual <= 7551:
        raise ValueError("actual fitting anchor must be4097..7551")
    guard()
    template,control,moments,sizes,pins = _dependencies(guard)
    saved = sys.dont_write_bytecode
    sys.dont_write_bytecode = True
    try:
        import numpy as np
    finally:
        sys.dont_write_bytecode = saved
    actual,fingerprint,positives,zero_count = _prepare(rational_degrees,guard)
    holes_table = _reverse(actual,64)
    logs = np.full((65,65,6), -np.inf, dtype=float)
    for (f,s),row in holes_table.items():
        guard()
        logs[f,s] = tuple(_log_fraction(value) for value in row)
        if actual[64-f-s,s] != row:
            raise ArithmeticError("hole relabeling changed a channel")
    if np.isnan(logs).any() or np.isposinf(logs).any():
        raise ArithmeticError("invalid proposal log table")
    del actual,holes_table
    objective = control.Objective(logs,0.0,moments,np,guard)
    holes = 8192-anchor_actual
    bank,fit_reports,selected = [],[],[]
    failures = []
    total_evaluations = 0
    bank_fingerprint = hashlib.sha256()
    for center,lo,hi,_,_ in template:
        guard()
        full,single = fit_target(holes,center)
        objective.targets = 256*np.asarray([float(full),float(single)])
        seed = tuple(_positive_rational(value) for value in sizes.markers(holes,full))
        marker = seed
        history = []
        status = "positive rational seed retained"
        rounded = None
        failure = None
        try:
            best,history,status = control.coordinate_fit(objective,
                np.asarray([_log_fraction(value) for value in seed]),np,guard)
            marker = tuple(_positive_rational(F.from_float(math.exp(float(x)))) for x in best["x"])
            # One fixed witness check follows the bounded proposal procedure.
            rounded = objective.evaluate(np.asarray([_log_fraction(value) for value in marker]))
        except (ArithmeticError, ValueError) as error:
            # No center is omitted.  These are proposals, not proof terms.
            marker = seed
            failure = type(error).__name__+": "+str(error)
            status = "FLOAT fit failed; positive rational seed retained"
            failures.append({"center":center,"reason":failure})
        if len(history) > MAX_PROPOSALS:
            raise ArithmeticError("frozen coordinate fit exceeded12 proposals")
        if not all(value > 0 and value.denominator <= DENOMINATOR_LIMIT for value in marker):
            raise ArithmeticError("invalid exported rational witness")
        row = (center,lo,hi,*marker)
        bank.append(row)
        bank_fingerprint.update((",".join(map(str,row))+"\n").encode("ascii"))
        total_evaluations += len(history)
        report = {"center":center,"evaluations":len(history),"status":status}
        if rounded is not None and history:
            report["ordinary_fit_gain_bits_FLOAT"] = history[0]["log2_ordinary_FLOAT"]-rounded["objective"]/math.log(2)
        fit_reports.append(report)
        if center in (0,256,512,768,1024):
            detail = {"center":center,"fit_targets":[str(full),str(single)]}
            if rounded is not None:
                detail["rounded_region_means_FLOAT"] = [float(x/256) for x in rounded["means"]]
                detail["rounded_region_gradient_FLOAT"] = [float(x/256) for x in rounded["gradient"]]
            selected.append(detail)
        if progress:
            progress("hole_marker_fit",center=center,status=status,evaluations=len(history))
    guard()
    if len(bank) != 31 or tuple(row[:3] for row in bank) != tuple(row[:3] for row in template):
        raise ArithmeticError("fitting changed reference coverage")
    return tuple(bank), {
        "scope":"FLOAT-only witness selection from an already-directed rational local table",
        "certificate":False,"actual_anchor":anchor_actual,"hole_anchor":holes,
        "reference_hole_weight":REFERENCE_WEIGHT,"reference_centers":31,
        "local_cap_builds":0,"degree_reversals":1,"local_degree_rows":2145,
        "actual_rational_coefficients_sha256":fingerprint,
        "positive_coefficients_per_channel":positives,"zero_coefficients":zero_count,
        "frozen_coordinate_evaluation_limit":MAX_PROPOSALS,
        "successful_proposal_evaluations_recorded":total_evaluations,
        "one_rounded_witness_check_after_each_successful_fit":True,
        "marker_denominator_limit":DENOMINATOR_LIMIT,"positive_floor":"1/1048576",
        "objective_numerator":"e0 H^16383 J","objective_targets":"256*(full,single)",
        "marker_independent_objective_constant":0,"source_or_count_factor_applied":False,
        "new_local_alpha":False,"endpoint_proxies_only_for_fitting":True,
        "bank_sha256":bank_fingerprint.hexdigest(),"fits":fit_reports,
        "selected_fit_diagnostics":selected,"fit_failures":failures,"pins":pins}


def exact_tests(guard):
    """Tiny rational interface tests plus clearly labeled FLOAT log checks."""
    checks = 0

    def check(condition):
        nonlocal checks
        assert condition
        checks += 1

    toy = {(f,s):tuple(F(1+f+3*s+channel,17+channel) for channel in range(6))
           for f in range(3) for s in range(3-f)}
    toy[0,0] = (F(0),)+toy[0,0][1:]
    original = dict(toy)
    prepared,_,positive,zero_count = _prepare(toy,guard,capacity=2)
    reversed_toy = _reverse(prepared,2)
    check(_reverse(reversed_toy,2) == prepared)
    check(toy == original and all(x > 0 for x in positive) and zero_count == 1)
    u,v = F(2,3),F(5,7)
    for channel in range(6):
        guard()
        lhs = sum(row[channel]*u**f*v**s for (f,s),row in reversed_toy.items())
        rhs = u**2*sum(row[channel]*(1/u)**f*(v/u)**s for (f,s),row in prepared.items())
        check(lhs == rhs)
        for (f,s),row in prepared.items():
            check(reversed_toy[2-f-s,s][channel] == row[channel])
    for actual_anchor in (4097,5120,6144,7551):
        holes = 8192-actual_anchor
        for center in BASE_CENTERS:
            guard(); full,single = fit_target(holes,center)
            check(full > 0 and single > 0 and 2*full+single == holes)
            if center not in (0,1024):
                check(full == F(holes*center,2048))
            else:
                check(full == (F(1,2) if center == 0 else F(holes-1,2)))
            # Rational targets are proposals; no factorial of these values is taken.
            check(256*full+128*single == 128*holes)
    for value in (F(1,1<<80),F(1,DENOMINATOR_LIMIT),F(1,3),F(7,9),F(17),F(100001,17)):
        guard(); marker = _positive_rational(value)
        check(marker > 0 and marker.denominator <= DENOMINATOR_LIMIT)
        check(F(marker.numerator,marker.denominator) == marker)
    check(_positive_rational(F(1,1<<80)) == F(1,DENOMINATOR_LIMIT))
    check(_log_fraction(F(0)) == -math.inf)
    check(_log_fraction(F(1)) == 0.0)
    for exponent in (-4096,-128,1,128,4096):
        value = F(2)**exponent
        observed = _log_fraction(value)
        check(math.isfinite(observed))
        check(abs(observed-exponent*math.log(2)) < 1e-9)
    for bad in (F(-1),):
        try:
            _log_fraction(bad)
        except ValueError:
            checks += 1
        else:
            raise AssertionError("negative log input accepted")
    for bad in (640,4096):
        try:
            fit_target(bad,512)
        except ValueError:
            checks += 1
        else:
            raise AssertionError("invalid hole anchor accepted")
    return {"checks":checks,"rational_reversal_evaluation_identity":True,
            "all_six_channel_order_preserved":True,"odd_actual_anchors_tested":True,
            "rational_fitting_targets_not_actual_profiles":True,
            "log_conversion_tests_are_FLOAT_only":True,"numerical_bound_verified":False}
