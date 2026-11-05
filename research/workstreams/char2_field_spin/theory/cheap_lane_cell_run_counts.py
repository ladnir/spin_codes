"""Run-interval count envelopes for normalized ordinary marker cells.

The caller has a common local table, z and alpha.  Its conditional cell
log bounds already absorb beta^(alpha*r), through256 regional factors,
but NOT eta_r^alpha, indexed support counts, or z^(-alpha*L).

At actual size r, put k=r in actual coordinates and k=8192-r in hole
coordinates.  A marginal f belongs to the unique normalized cell with
x=2048*f/k: lower cuts are strict except at0, upper cuts are closed.
For assigned integer ranges [a0,a1],[b0,b1], all supports in that ordered
cell have c=k-a-b in [k-a1-b1,k-a0-b0], clipped to positive feasible runs.
The indexed count N_r(c) equals N_k(c) under complementation.  Summing it
over this interval is an exact RUN-INTERVAL mass and only an UPPER bound
on the count of that cell's profiles.  Different cells can overlap.

source_counts(...) streams one r's integer run prefix at a time and returns
496 masses for31 cells, weighted by the ACTUAL-r eta_r^alpha.  combine(...)
adds both terminal orders off diagonal and one on diagonal.  It also
returns the old global conditional maximum times the full support sum.
Both returned log bounds EXCLUDE the common output correction.  No extra
source, profile, marker or orientation factor belongs in the caller.

No provider imports, main, or work on import.  The caller owns all resource
guards and must supply valid directed conditional bounds for its geometry.
The helper validates counting coverage, not the local-cap theorem.
"""

from collections import Counter
from collections.abc import Mapping
from fractions import Fraction as F
import hashlib
from math import comb

POINTS = 8192
NORMALIZED = 2048
MAX_CELLS = 31


def _fraction(value):
    if isinstance(value, bool) or not isinstance(value, (int, F)):
        raise TypeError("exact integer or Fraction required")
    return F(value)


def _cells(cells):
    cells = tuple(cells)
    if not 1 <= len(cells) <= MAX_CELLS:
        raise ValueError("require1..31 normalized cells")
    out = []; previous = F(0); centers = set()
    for row in cells:
        if not isinstance(row, (tuple, list)) or len(row) != 5:
            raise TypeError("expected(center,lower,upper,u,v)")
        center, low, high, u, v = row
        if type(center) is not int or center in centers:
            raise ValueError("cell centers must be distinct integer labels")
        low, high, u, v = map(_fraction, (low, high, u, v))
        if low != previous or not 0 <= low < high <= 1024 or min(u, v) <= 0:
            raise ValueError("positive markers and contiguous normalized cuts required")
        out.append((center, low, high, u, v)); centers.add(center); previous = high
    if previous != 1024:
        raise ValueError("cells must cover the full continuous interval0..1024")
    return tuple(out)


def marginal_ranges(cells, k, guard):
    """Unique assigned integer ranges; empty cells are None, ties go lower.

    cells must have passed _cells.  Every returned integer range is inside
    the caller's CLOSED geometry cell, despite unique tie assignment here.
    """
    if type(k) is not int or not 1 <= k <= POINTS//2:
        raise ValueError("coordinate weight must be1..4096")
    result = []; next_full = 0
    for index, (_, low, high, _, _) in enumerate(cells):
        guard()
        lo = 0 if index == 0 else (k*low/NORMALIZED).__floor__()+1
        hi = min(k//2, (k*high/NORMALIZED).__floor__())
        if lo <= hi:
            if lo != next_full:
                raise ArithmeticError("integer marginal assignment has a gap or overlap")
            if not low <= F(NORMALIZED*lo, k) <= F(NORMALIZED*hi, k) <= high:
                raise ArithmeticError("assigned range leaves closed geometry")
            if index and not low < F(NORMALIZED*lo, k):
                raise ArithmeticError("normalized-cut tie was assigned twice")
            result.append((lo, hi)); next_full = hi+1
        else:
            result.append(None)
    if next_full != k//2+1:
        raise ArithmeticError("incomplete marginal assignment")
    return tuple(result)


def _run_prefix(counts, total, guard):
    prefix = [0]; running = 0
    if not counts or counts[0] != 0:
        raise ArithmeticError("run vector must start with zero")
    for c, value in enumerate(counts[1:], 1):
        if c % 64 == 1: guard()
        if type(value) is not int or value <= 0:
            raise ArithmeticError("positive exact indexed run counts required")
        running += value; prefix.append(running)
    if running != total:
        raise ArithmeticError("own-size run counts fail the binomial checksum")
    return prefix


def _mass_row(prefix, ranges, k, pairs, guard):
    """Exact unordered interval masses and ordered overlap diagnostics."""
    if len(prefix) != k+1:
        raise ArithmeticError("run prefix has the wrong coordinate weight")
    masses = []; intervals = []; difference = [0]*(k+2)
    ordered_mass = 0
    for index, (i, j) in enumerate(pairs):
        if index % 32 == 0: guard()
        left, right = ranges[i], ranges[j]
        interval = None; mass = 0
        if left is not None and right is not None:
            lo = max(1, k-left[1]-right[1])
            hi = min(k, k-left[0]-right[0])
            if lo <= hi:
                interval = (lo, hi)
                mass = prefix[hi]-prefix[lo-1]
                multiplicity = 1 if i == j else 2
                difference[lo] += multiplicity; difference[hi+1] -= multiplicity
                ordered_mass += multiplicity*mass
        masses.append(mass); intervals.append(interval)
    overlap = maximum = 0
    for c in range(1, k+1):
        if c % 64 == 1: guard()
        overlap += difference[c]
        if not 1 <= overlap <= len(ranges)**2:
            raise ArithmeticError("run intervals miss a class or exceed the cell count")
        maximum = max(maximum, overlap)
    if not prefix[-1] <= ordered_mass <= len(ranges)**2*prefix[-1]:
        raise ArithmeticError("duplicate count-mass control failed")
    return tuple(masses), tuple(intervals), maximum, ordered_mass


def source_counts(v, anchor, join, cells, alpha, first, last, coordinate, guard):
    """Return pairs/weights/full_weight/report; counts outside alpha.

    weights[t]=sum_r eta_r^alpha*sum_(c in assigned interval_t(r)) N_r(c).
    full_weight=sum_r eta_r^alpha*C(8192,r).  No beta or output factor.
    The supplied join.run_counts(points=...,weight=...) is the frozen exact
    recurrence; holes always query the complementary SMALL count weight.
    """
    cells = _cells(cells); alpha = _fraction(alpha)
    if not 0 < alpha <= 1:
        raise ValueError("alpha must be in(0,1]")
    if type(first) is not int or type(last) is not int or not 1 <= first <= last < POINTS:
        raise ValueError("nonempty inclusive actual-size band required")
    if coordinate not in ("actual", "hole"):
        raise ValueError("coordinate must be actual or hole")
    if coordinate == "actual" and last > POINTS//2:
        raise ValueError("actual coordinates require r<=4096")
    if coordinate == "hole" and first < POINTS//2:
        raise ValueError("hole coordinates require r>=4096")
    pairs = tuple((i, j) for i in range(len(cells)) for j in range(i, len(cells)))
    weights = [v.exact(F(0)) for _ in pairs]
    full_weight = v.exact(F(0)); exponent = v.exact(alpha)
    exact_cell_sums = [0]*len(pairs)
    indexed_sum = duplicate_sum = marginal_total = 0
    maximum_overlap = 0; maximum_at = None; factors = {}
    indexed_hash = hashlib.sha256(); assignment_hash = hashlib.sha256()
    interval_hash = hashlib.sha256(); source_hash = hashlib.sha256()
    for r in range(first, last+1):
        guard(); k = r if coordinate == "actual" else POINTS-r
        count_weight = min(r, POINTS-r)
        if k != count_weight:
            raise ArithmeticError("coordinate/count complementation mismatch")
        counts = join.run_counts(points=POINTS, weight=count_weight)
        total = comb(POINTS, r)
        if total != comb(POINTS, count_weight):
            raise ArithmeticError("actual/complement binomial identity failed")
        prefix = _run_prefix(counts, total, guard)
        if len(counts) != k+1:
            raise ArithmeticError("unexpected indexed run-count length")
        ranges = marginal_ranges(cells, k, guard)
        masses, intervals, overlap, duplicated = _mass_row(prefix, ranges, k, pairs, guard)
        h, c_h, eta = anchor.source_factor(r)
        c_h, eta = _fraction(c_h), _fraction(eta)
        if type(h) is not int or h < 0 or c_h <= 0 or eta <= 0:
            raise ArithmeticError("invalid actual-r source factor")
        factor = v.exact(eta)**exponent
        if not factor.is_finite() or not factor > 0:
            raise ArithmeticError("nonpositive directed eta factor")
        full_weight += total*factor
        for index, mass in enumerate(masses):
            if index % 32 == 0: guard()
            if mass: weights[index] += mass*factor
            exact_cell_sums[index] += mass
        indexed_sum += total; duplicate_sum += duplicated; marginal_total += k//2+1
        if overlap > maximum_overlap: maximum_overlap, maximum_at = overlap, r
        entry = factors.setdefault(h, {"first_actual":r, "last_actual":r, "sizes":0})
        entry["last_actual"] = r; entry["sizes"] += 1
        indexed_hash.update((str(r)+":"+str(total)+"\n").encode("ascii"))
        assignment_hash.update((str(r)+":"+str(k)+":"+repr(ranges)+"\n").encode("ascii"))
        interval_hash.update((str(r)+":"+repr(intervals)+"\n").encode("ascii"))
        source_hash.update((str(r)+":"+str(h)+":"+str(c_h)+":"+str(eta)+"\n").encode("ascii"))
        # Do not retain one large prefix/count vector per actual source size.
        del counts, prefix, ranges, masses, intervals
    guard()
    if not full_weight.is_finite() or not full_weight > 0:
        raise ArithmeticError("invalid full weighted source count")
    for value in weights:
        if not value.is_finite() or not value >= 0:
            raise ArithmeticError("invalid positive run-interval mass")
    mass_hash = hashlib.sha256()
    for pair, value in zip(pairs, exact_cell_sums):
        mass_hash.update((str(pair)+":"+str(value)+"\n").encode("ascii"))
    report = {"actual_first":first, "actual_last":last, "sizes":last-first+1,
        "coordinate":coordinate, "odd_sizes":sum(r%2 for r in range(first,last+1)),
        "cells":len(cells), "unordered_cells":len(pairs), "ordered_cells":len(cells)**2,
        "integer_marginals_checked":marginal_total, "own_r_binomial_checksum":True,
        "indexed_support_count":str(indexed_sum), "indexed_count_sha256":indexed_hash.hexdigest(),
        "assignment_sha256":assignment_hash.hexdigest(), "run_interval_sha256":interval_hash.hexdigest(),
        "unweighted_cell_mass_sha256":mass_hash.hexdigest(), "actual_source_sha256":source_hash.hexdigest(),
        "source_h_ranges":factors, "tie_policy":"lower cell; other lower cuts strict",
        "cell_counts_are_upper_bounds":True, "run_interval_sums_are_exact":True,
        "ordered_run_interval_overlap_max":maximum_overlap, "overlap_max_actual":maximum_at,
        "duplicate_count_mass_le_cells_squared":duplicate_sum <= len(cells)**2*indexed_sum,
        "count_overlap_diagnostic_only":True, "no_uniform_weighted_loss_claim":True,
        "counts_outside_alpha":True, "eta_uses_actual_r":True,
        "beta_factor_applied_here":False, "output_factor_applied_here":False,
        "one_actual_size_prefix_retained":True}
    return {"pairs":pairs, "weights":tuple(weights), "full_weight":full_weight,
        "centers":tuple(row[0] for row in cells), "report":report}


def combine(v, counts, cell_bounds, guard):
    """Combine complete directed log bounds; common output factor stays out.

    cell_bounds[(i,j)] is (("forward", logupper), ("reverse", logupper));
    diagonal cells have ONLY forward.  Empty geometry may be omitted, but
    every positive run-interval mass requires its corresponding bounds.
    The old maximum includes ALL supplied geometry cells, even zero-mass
    ones, matching the frozen max-conditional-times-full-count control.
    """
    if not isinstance(cell_bounds, Mapping):
        raise TypeError("cell bounds must be a mapping")
    pairs = counts["pairs"]; weights = counts["weights"]
    if len(pairs) != len(weights) or len(set(pairs)) != len(pairs):
        raise ValueError("invalid count-vector layout")
    if set(cell_bounds)-set(pairs):
        raise ValueError("unknown or reversed unordered cell index")
    maximum = None; old_worst = None; weighted = v.exact(F(0)); leaders = []
    used = ordered = zero_mass = 0
    logtwo = v.arb(2).log()
    for (i,j), weight in zip(pairs, weights):
        guard()
        bounds = cell_bounds.get((i,j))
        if not weight.is_finite() or not weight >= 0:
            raise ArithmeticError("invalid count weight")
        if bounds is None:
            if weight != 0:
                raise ArithmeticError("positive-count cell has no geometry bound")
            continue
        bounds = tuple(bounds)
        expected = ("forward",) if i == j else ("forward", "reverse")
        if tuple(row[0] for row in bounds) != expected:
            raise ValueError("diagonal once; off-diagonal needs BOTH terminal orders")
        subtotal = v.exact(F(0))
        for direction, logbound in bounds:
            if not logbound.is_finite():
                raise ArithmeticError("finite directed conditional log bound required")
            if maximum is None or logbound.upper() > old_worst[0]:
                old_worst = (logbound.upper(), (i,j), direction)
            maximum = logbound if maximum is None else maximum.max(logbound)
            if weight != 0: subtotal += logbound.exp()
            ordered += 1
        if weight == 0:
            zero_mass += 1; continue
        term = weight*subtotal
        if not term.is_finite() or not term > 0:
            raise ArithmeticError("invalid positive weighted cell contribution")
        weighted += term; used += 1
        leaders.append((term.upper(), {"cell_indices":[i,j],
            "centers":[counts["centers"][i],counts["centers"][j]],
            "terminal_orders":len(bounds), "log2_weighted_cell_interval":str(term.log()/logtwo)}))
    if maximum is None or not weighted > 0 or not counts["full_weight"] > 0:
        raise ArithmeticError("empty or nonpositive completed envelope")
    old = counts["full_weight"].log()+maximum
    new = weighted.log()
    report = {"supplied_geometry_cells":len(cell_bounds), "supplied_terminal_orders":ordered,
        "positive_count_cells":used, "supplied_zero_count_cells":zero_mass,
        "old_max_includes_zero_mass_geometry":True,
        "old_worst_cell_indices":list(old_worst[1]), "old_worst_order":old_worst[2],
        "largest_run_weighted_cells":[row for _,row in sorted(leaders,key=lambda x:x[0],reverse=True)[:8]],
        "both_terminal_orders_added":True, "counts_are_run_interval_upper_bounds":True,
        "beta_and_output_not_applied_again":True, "output_correction_still_required":True,
        "minimum_is_of_completed_valid_bounds":True,
        "log2_run_minus_old_interval":str((new-old)/logtwo)}
    return {"old_log_bound":old, "run_log_bound":new,
        "minimum_log_bound":old.min(new), "maximum_conditional_log_bound":maximum,
        "report":report}


def exact_tests(join, guard):
    """Small exhaustive cycles and Fraction unions; no proof providers loaded."""
    checks = masks = weighted_cases = 0
    cells = _cells(((128,F(0),F(256),F(1),F(1)),
                    (384,F(256),F(512),F(1),F(1)),
                    (768,F(512),F(1024),F(1),F(1))))
    pairs = tuple((i,j) for i in range(3) for j in range(i,3))
    assert marginal_ranges(cells,8,guard) == ((0,1),(2,2),(3,4)); checks += 1
    regions = tuple(((F(i+2,11),F(1,7),F(3,13),F(i+5,17)),
                     (F(i+7,19),F(2*i+3,23))) for i in range(3))
    def direct_chain(left, right, bits):
        row = (F(1),F(0)); a, _ = left; b, terminal = right
        for matrix in (a,b)*((1<<bits)-1)+(a,):
            row = (row[0]*matrix[0]+row[1]*matrix[2],
                   row[0]*matrix[1]+row[1]*matrix[3])
        return row[0]*terminal[0]+row[1]*terminal[1]
    for i,j in pairs:
        for bits in (1,2,3):
            guard()
            forward,reverse = join.two_orientations(regions[i],regions[j],bits=bits,diagonal=i==j)
            assert forward == direct_chain(regions[i],regions[j],bits)
            assert reverse == direct_chain(regions[j],regions[i],bits); checks += 2
    for points in (6,8,10):
        actual = Counter(); profiles = Counter()
        for mask in range(1, (1<<points)-1):
            guard(); bits = tuple((mask>>j)&1 for j in range(points)); r = sum(bits)
            c = sum(bits[j] and not bits[j-1] for j in range(points))
            a = sum(bits[j]*bits[(j+1)%points] for j in range(0,points,2))
            b = sum(bits[j]*bits[(j+1)%points] for j in range(1,points,2))
            assert c == r-a-b; checks += 1
            ha, hb = points//2-r+a, points//2-r+b
            assert c == points-r-ha-hb; checks += 1
            assert r-2*a == points-r-2*ha and r-2*b == points-r-2*hb; checks += 1
            actual[r,c] += 1; profiles[r,a,b] += 1; masks += 1
        total_true = total_upper = total_old = F(0)
        for r in range(1,points):
            guard(); k = min(r,points-r)
            run = join.run_counts(points=points,weight=k)
            assert tuple(run) == tuple(actual[r,c] for c in range(k+1)); checks += k+1
            prefix = _run_prefix(run,comb(points,r),guard)
            assert prefix == [sum(run[:c+1]) for c in range(k+1)]; checks += k+1
            ranges = marginal_ranges(cells,k,guard)
            masses, intervals, overlap, duplicate = _mass_row(prefix,ranges,k,pairs,guard)
            assigned = Counter(); index = {}
            for i, interval in enumerate(ranges):
                if interval is not None:
                    for f in range(interval[0],interval[1]+1): index[f] = i
            for (rr,a,b), number in profiles.items():
                if rr != r: continue
                aa,bb = (a,b) if r <= points//2 else (points//2-r+a,points//2-r+b)
                i,j = index[aa],index[bb]; assigned[i,j] += number
                position = pairs.index((min(i,j),max(i,j)))
                lo,hi = intervals[position]
                assert lo <= r-a-b <= hi; checks += 1
            assert sum(assigned.values()) == comb(points,r); checks += 1
            exact = upper = F(0); largest = F(0)
            for (i,j),mass in zip(pairs,masses):
                assert assigned[i,j] <= mass and assigned[j,i] <= mass; checks += 2
                forward,reverse = join.two_orientations(regions[i],regions[j],bits=1+r%3,diagonal=i==j)
                exact += assigned[i,j]*forward
                upper += mass*forward; largest = max(largest,forward)
                if i != j:
                    exact += assigned[j,i]*reverse
                    upper += mass*reverse; largest = max(largest,reverse)
            old = comb(points,r)*largest
            assert exact <= upper and exact <= old; checks += 2
            eta = F(r+2,17)
            total_true += eta*exact; total_upper += eta*upper; total_old += eta*old
            assert 1 <= overlap <= 9 and duplicate <= 9*comb(points,r); checks += 2
            weighted_cases += 1
        assert total_true <= min(total_upper,total_old); checks += 1
    # Empty marginal cells and odd coordinate endpoints remain in coverage.
    for k in (1,2,3,7,8,9):
        ranges = marginal_ranges(cells,k,guard)
        flat = [f for interval in ranges if interval for f in range(interval[0],interval[1]+1)]
        assert flat == list(range(k//2+1)); checks += 1
    return {"exact_checks":checks, "exhaustive_nonconstant_masks":masks,
        "weighted_union_cases":weighted_cases, "actual_and_hole_profile_mapping":True,
        "indexed_run_counts_checked_by_enumeration":True,
        "unique_cut_ties_and_empty_cells_checked":True,
        "fraction_noncommuting_terminal_chains":True,
        "cell_counts_tested_as_upper_bounds_not_equalities":True}
