"""Input/output moments with explicit nonzero-tile feedback constraints.

Free tile moments use exact split MDS enumerators. The quiet upper bound
counts each pattern of nonzero input tiles, retains its feedback rank, and
removes the exact all-zero solution of that pattern. It may overcount
proper subpatterns, hence is an upper bound rather than an exact enumerator.
"""

from fractions import Fraction as F
from functools import lru_cache
from math import comb, exp2, factorial, lgamma, log, log2, prod

from hd_randomized_outer_proof import prescribed_mds_count
from hd_randomized_outer_quiet_tile import profile_moments as original_profiles
from hd_randomized_outer_quiet_tile import marker_screen
from hd_randomized_outer_state_screen import NEG, add


@lru_cache(None)
def split_enumerator(a, width=8, field=65536):
    """Exact integer counts, indexed by input weight and output weight."""
    if not 0 <= a <= width:
        raise ValueError('Invalid potential input count')
    return tuple(tuple(comb(a, i)*comb(width, w)
                       *prescribed_mds_count(i+w, a+width, a, field)
                       for w in range(width+1)) for i in range(a+1))


def join(left, right, state=4, distinct=True):
    if distinct and state == 4:
        if left == right == 1:
            return 2
        if left == right == 3:
            return 4
    return max(left, right)


def local_integers(mark, input_mark, width=8, state=4, field=65536):
    """Numerators with common denominator z_den^width v_den^a q^a.

    inactive includes every selected output-zero mask. active arrays contain
    only masks admitting a nonzero input. q_weight counts residual solutions;
    unit_weight removes their known zero solution after the feedback join.
    """
    z, v = F(mark), F(input_mark)
    if not 0 < z <= 1 or v < 1 or not 1 <= state <= width:
        raise ValueError('Invalid marks or state dimension')
    result = []
    for a in range(width+1):
        q_weight, unit_weight = [0]*(state+1), [0]*(state+1)
        for zeros in range(a):
            dimension = a-zeros
            rank = min(state, dimension)
            coefficient = (comb(width, zeros)*z.numerator**(width-zeros)
                           *(z.denominator-z.numerator)**zeros*v.numerator**a)
            q_weight[rank] += coefficient*field**dimension
            unit_weight[rank] += coefficient
        inactive = z.denominator**width*v.denominator**a
        free = sum(count*v.numerator**i*v.denominator**(a-i)
                   *z.numerator**w*z.denominator**(width-w)
                   for i, row in enumerate(split_enumerator(a, width, field))
                   for w, count in enumerate(row))
        result.append((inactive, tuple(q_weight), tuple(unit_weight), free))
    return tuple(result)


def marked_pattern(profile, mark, input_mark, *, width=8, state=4,
                   field=65536, distinct=True):
    """One exact-rational upper moment; used for independent small tests."""
    local = local_integers(mark, input_mark, width, state, field)
    q_dp, unit_dp = [1]+[0]*state, [1]+[0]*state
    free = 1
    for a in profile:
        inactive, q_local, unit_local, tile_free = local[a]
        following_q, following_unit = [x*inactive for x in q_dp], [x*inactive for x in unit_dp]
        for r in range(state+1):
            for s in range(1, state+1):
                successor = join(r, s, state, distinct)
                following_q[successor] += q_dp[r]*q_local[s]
                following_unit[successor] += unit_dp[r]*unit_local[s]
        q_dp, unit_dp = following_q, following_unit
        free *= tile_free
    j = sum(profile)
    common = F(mark).denominator**(width*len(profile))*F(input_mark).denominator**j
    nonzero = sum(q_dp[r]*field**(state-r)-unit_dp[r]*field**state for r in range(state+1))
    if nonzero < 0:
        raise ArithmeticError('Nonzero-pattern count became negative')
    quiet = F(common*field**state+nonzero, common*field**(j+state))
    terminal = F(free, common*field**j)
    return min(quiet, terminal), terminal


@lru_cache(maxsize=12)
def profile_moments(mark, input_mark=F(1), field=65536):
    """All sorted FFT8 profiles, with shared-prefix positive integer DPs."""
    z, v = F(mark), F(input_mark)
    local = local_integers(z, v, field=field)
    base = {profile:(quiet, free) for profile, j, mult, quiet, free
            in original_profiles(z, field, distinct_hyperplanes=True)}
    result, counts = [], [0]*9

    def visit(profile, start, j, multiplicity, q_dp, unit_dp, free):
        if len(profile) == 8:
            multiplicity *= factorial(8)//prod(factorial(x) for x in counts)
            common = z.denominator**64*v.denominator**j
            nonzero = sum(q_dp[r]*field**(4-r)-unit_dp[r]*field**4 for r in range(5))
            if nonzero < 0:
                raise ArithmeticError('Nonzero-pattern count became negative')
            quiet = F(common*field**4+nonzero, common*field**(j+4))
            terminal = F(free, common*field**j)
            zero = F(1, field**j)
            old_quiet, old_free = base[profile]
            coarse = zero+v**j*(old_quiet-zero)
            quiet = min(quiet, coarse, terminal)
            if quiet < zero or terminal < zero:
                raise ArithmeticError('Lost the exact zero atom')
            result.append((profile, j, multiplicity, quiet, terminal))
            return
        for a in range(start, 9):
            inactive, q_local, unit_local, tile_free = local[a]
            following_q = [x*inactive for x in q_dp]
            following_unit = [x*inactive for x in unit_dp]
            for r in range(5):
                for s in range(1, 5):
                    successor = join(r, s)
                    following_q[successor] += q_dp[r]*q_local[s]
                    following_unit[successor] += unit_dp[r]*unit_local[s]
            counts[a] += 1
            visit(profile+(a,), a, j+a, multiplicity*comb(8, a),
                  following_q, following_unit, free*tile_free)
            counts[a] -= 1

    visit((), 0, 0, 1, [1, 0, 0, 0, 0], [1, 0, 0, 0, 0], 1)
    return tuple(result)


def fractional_logs(mark, alpha, input_mark=F(1)):
    alpha = F(alpha)
    if not 0 < alpha <= 1:
        raise ValueError('Invalid fractional exponent')
    quiet, free = [NEG]*65, [NEG]*65
    for profile, j, multiplicity, k, t in profile_moments(F(mark), F(input_mark)):
        weight = log2(multiplicity)
        quiet[j] = add(quiet[j], weight+float(alpha)*(log2(k.numerator)-log2(k.denominator)))
        free[j] = add(free[j], weight+float(alpha)*(log2(t.numerator)-log2(t.denominator)))
    return (tuple(x-log2(comb(64, j)) for j, x in enumerate(quiet)),
            tuple(x-log2(comb(64, j)) for j, x in enumerate(free)))


def screen(rows, mark, alpha, input_mark):
    return marker_screen(rows, mark, alpha, *fractional_logs(mark, alpha, input_mark),
                         input_mark=input_mark)


def fractional_enclosures(mark, alpha, input_mark=F(1)):
    from flint import arb
    from hd_randomized_outer_quiet_exact import _arb_fraction
    alpha = F(alpha)
    if not 0 < alpha <= 1:
        raise ValueError('Invalid fractional exponent')
    exponent = _arb_fraction(alpha)
    quiet, free = [arb(0) for _ in range(65)], [arb(0) for _ in range(65)]
    for profile, j, multiplicity, k, t in profile_moments(F(mark), F(input_mark)):
        quiet[j] += multiplicity*_arb_fraction(k)**exponent
        free[j] += multiplicity*_arb_fraction(t)**exponent
    return (tuple(x/comb(64, j) for j, x in enumerate(quiet)),
            tuple(x/comb(64, j) for j, x in enumerate(free)))


def directed_endpoints(rows, mark, alpha, input_mark, marker=None, coefficient_cap=None, *, enclosures=None):
    """Replay several row counts with one fixed positive marker witness."""
    from flint import arb
    from hd_randomized_outer_quiet_exact import _arb_fraction
    rows = tuple(rows)
    z, a, v = map(_arb_fraction, map(F, (mark, alpha, input_mark)))
    if marker is None:
        proposal = screen(rows[len(rows)//2], mark, alpha, input_mark)
        marker = F(1/(1+exp2(-proposal['marker_log_odds']))).limit_denominator(1 << 40)
    marker = F(marker)
    if not 0 < marker < 1:
        raise ValueError('Invalid marker')
    m = _arb_fraction(marker)
    quiet, free = fractional_enclosures(mark, alpha, input_mark) if enclosures is None else enclosures
    ordinary = sum((comb(64, j)*m**j*(1-m)**(64-j)*quiet[j] for j in range(65)), arb(0))
    terminal = sum((comb(64, j)*m**j*(1-m)**(64-j)*free[j] for j in range(65)), arb(0))
    logarithm = lambda value: value.log()/arb(2).log()
    coefficient_correction = arb(0)
    if coefficient_cap is not None:
        from flint import arb_poly
        cap = _arb_fraction(F(coefficient_cap))
        if not 0 < cap <= 1:
            raise ValueError('Invalid coefficient cap')
        u = arb_poly([comb(64, j)*m**j*(1-m)**(64-j)*quiet[j]/ordinary for j in range(65)])
        final_poly = arb_poly([comb(64, j)*m**j*(1-m)**(64-j)*free[j]/terminal for j in range(65)])
        first, last = u**128, u**127*final_poly
        if not all(value < cap for poly in (first, last) for value in poly.coeffs()):
            raise ArithmeticError('The proposed uniform coefficient cap was not proved')
        coefficient_correction = 256*logarithm(cap)
    beta = 4096-128*logarithm(arb(65535))
    result = []
    for r in rows:
        if not 1 <= r <= 8192:
            raise ValueError('Invalid row count')
        choose = logarithm(arb(comb(8192, r)))
        mass = choose+r*logarithm(m)+(8192-r)*logarithm(1-m)
        bound = (choose+a*(r*beta-logarithm(arb(65535))-838860*logarithm(z)-129*r*logarithm(v))
                 +32767*logarithm(ordinary)+logarithm(terminal)-256*mass+coefficient_correction)
        result.append(dict(active_rows=r, log2_upper=str(bound), below_minus_40=bool(bound < -40)))
    return dict(mark=str(mark), alpha=str(alpha), input_mark=str(input_mark), marker=str(marker),
                coefficient_cap=None if coefficient_cap is None else str(coefficient_cap),
                endpoints=result, directed=True, whole_code_certificate=False)


def directed_interval(first, last, mark, alpha, input_mark, marker, coefficient_cap=None, *, enclosures=None):
    from flint import arb
    if not 1 <= first <= last <= 8192:
        raise ValueError('Invalid interval')
    receipt = directed_endpoints((first, last), mark, alpha, input_mark, marker, coefficient_cap,enclosures=enclosures)
    left, right = (arb(e['log2_upper']) for e in receipt['endpoints'])
    # Keep an enclosing hull if endpoint intervals overlap.
    maximum = left if left > right else right if right > left else left.union(right)
    union = maximum+arb(last-first+1).log()/arb(2).log()
    receipt.update(active_rows_first=first, active_rows_last=last,
                   union_log2_upper=str(union), below_minus_40=bool(union < -40),
                   interior_justification='convex affine(r)-255 log2 binom(8192,r)')
    return receipt


def interval_proposal(center, mark, alpha, input_mark, threshold=-80, coefficient_cap=None, *, moments=None):
    """FLOAT-only interval selection; convexity is certified in the replay."""
    from hd_randomized_outer_quiet_exact import weighted_logsum_mean
    quiet, free = fractional_logs(mark, alpha, input_mark) if moments is None else moments
    proposal = marker_screen(center, mark, alpha, quiet, free, input_mark=input_mark)
    odds = proposal['marker_log_odds']
    marker = F(1/(1+exp2(-odds))).limit_denominator(1 << 40)
    ordinary, _ = weighted_logsum_mean(tuple(log2(comb(64, j))+quiet[j] for j in range(65)), odds)
    terminal, _ = weighted_logsum_mean(tuple(log2(comb(64, j))+free[j] for j in range(65)), odds)
    beta = 4096-128*log2(65535)
    accepted = []
    for r in range(1, 8193):
        choose = (lgamma(8193)-lgamma(r+1)-lgamma(8193-r))/log(2)
        bound = (float(alpha)*(r*beta-log2(65535)-838860*log2(mark)-129*r*log2(input_mark))
                 +32767*ordinary+terminal-256*r*odds-255*choose)
        if coefficient_cap is not None:
            bound += 256*log2(coefficient_cap)
        if bound < threshold:
            accepted.append(r)
    if not accepted:
        return None
    return dict(first=min(accepted), last=max(accepted), marker=str(marker),
                mark=str(mark), alpha=str(alpha), input_mark=str(input_mark),
                coefficient_cap=None if coefficient_cap is None else str(coefficient_cap))


def coverage_receipt():
    """Directed union of the nine new intervals, excluding prior receipts."""
    from flint import arb, ctx
    witnesses = (
        (1508,1644,6,16,'146211182168/577720884983',F(3,256)),
        (1544,1802,7,16,'172815637333/646891988433',F(3,256)),
        (1658,1998,9,16,'30758323013/106196160857',None),
        (1933,2466,16,16,'325686709449/951889216970',None),
        (2395,3070,32,16,'143989657028/347265555841',None),
        (2779,3709,64,16,'193629185397/402500592083',None),
        (3609,5215,256,8,'723745044563/1033346037698',None),
        (4318,6585,1024,8,'775311481896/966340120213',None),
        (5064,7574,4096,8,'441447453054/502077788525',None),
    )
    with ctx.workprec(192):
        receipts = [directed_interval(first,last,F(1,z),F(1,a),F(1),F(p),cap)
                    for first,last,z,a,p,cap in witnesses]
        probability = sum((arb(2)**arb(item['union_log2_upper']) for item in receipts),arb(0))
        union = probability.log()/arb(2).log()
        return dict(intervals=receipts, union_log2_upper=str(union),
                    below_minus_40=bool(union < -40), directed=True,
                    covered_first=1508, covered_last=7574, whole_code_certificate=False)


if __name__ == '__main__':
    import argparse
    import json
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--search', action='store_true')
    parser.add_argument('--receipt', action='store_true')
    args = parser.parse_args()
    if args.receipt:
        print(json.dumps(coverage_receipt()),flush=True)
        raise SystemExit(0)
    if args.search:
        rows = (1024, 1280, 1536, 2048, 3072, 4096, 5120, 5800)
        best = {r:None for r in rows}
        for shift in (3, 4, 6, 8, 10, 12, 14, 16):
            z = F(1, 2**shift)
            for v in (F(1), F(4), F(16), F(64), F(256), F(1024)):
                for alpha in (F(1, 64), F(1, 32), F(1, 16), F(1, 8)):
                    moments = fractional_logs(z, alpha, v)
                    for r in rows:
                        result = marker_screen(r, z, alpha, *moments, input_mark=v)
                        if best[r] is None or result['log2_upper_screen'] < best[r]['log2_upper_screen']:
                            best[r] = result
            print(json.dumps(dict(completed_shift=shift, best=list(best.values()))), flush=True)
        raise SystemExit(0)
    for rows, z in ((1536, F(1, 8)), (4096, F(1, 4096))):
        for v in (F(1), F(257, 256), F(17, 16), F(5, 4), F(2), F(4)):
            print(json.dumps(screen(rows, z, F(1, 32), v)), flush=True)
