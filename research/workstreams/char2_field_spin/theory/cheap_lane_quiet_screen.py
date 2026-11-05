"""Bounded screens and directed cover for the unchanged all-quiet subcode.

For fixed active groups and routed pair histograms (empty,single,full),
independent uniform reference symbols give the exact shortened-kernel
enumerator k_l(z). A half packet has sixteen pairs feeding two RS[16,4]
blocks. Fair swaps choose u of its s singletons for the first block.
The positive polynomial coefficient for (s,f) is
  multinomial(16-s-f,s,f) sum_u binom(s,u)/2**s
       * (k_(f+u)(z)*k_(f+s-u)(z))**alpha.
A region has 256 such halves. Its coefficient at the prescribed histogram,
divided by multinomial(4096; histogram), is the exact route mean of the
fractional quiet factor. Positive evaluation at any markers upper-bounds
that coefficient. Two alternating histograms each occur 128 times.

The sharper established Toeplitz outer density is 256**160/255**32 per
active group. Projectivization costs 1/255, before the fractional power.
The final packet is not required to have zero feedback. Since k_l(z) is
at least 256**-4 times the unconstrained moment, removing its four feedback
constraints costs at most 128*alpha bits. This is included, never omitted.

The profile/pair modes are FLOAT screens, not certificates. The --cover
mode proposes a complete cover, reevaluates every accepted vertex with Arb,
and independently checks all 8192 integer group-count slices. A successful
cover bounds only the true quiet subcode, never the whole code. This script
uses no generated inputs and writes no files. Optimizer limits are explicit.
Use --cover --include-witnesses to emit one standalone compact JSON receipt.
Use --verify-receipt PATH to recompute its coverage and directed bounds
without an optimizer or trust in the receipt's reported ceilings or hashes.
"""
import argparse
from hashlib import sha256
import json
from fractions import Fraction
from functools import lru_cache
from math import comb, isfinite, lgamma, log, log2
import re
from time import perf_counter
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp

from hd_shifted_pairs_count import profile_count_upper, one_pairing_count

Q, PAIRS, N = 256, 4096, 1 << 21
LOG_DENSITY = 1280 - 32 * log2(255)
LN2 = log(2)
ROWS = np.array([(s, f, u) for s in range(17) for f in range(17-s)
                 for u in range(s+1)], dtype=int)
S, F, U = ROWS.T
COUNTS = np.array([log2(comb(16, int(s))*comb(16-int(s), int(f))*comb(int(s), int(u))) - s
                   for s, f, u in ROWS])
DEG = np.arange(17)
profile_count_upper = lru_cache(256)(profile_count_upper)


def quiet_exact(l, z):
    z = Fraction(z)
    return sum(Fraction(comb(l, c), Q**(c+min(4, l-c))) * z**(l-c) * (1-z)**c
               for c in range(l+1))


def quiet_float(z):
    return np.array([sum(comb(l, c)*z**(l-c)*(1-z)**c*Q**(-c-min(4, l-c))
                         for c in range(l+1)) for l in range(17)])


def log_multinomial(hist):
    return (lgamma(sum(hist)+1)-sum(lgamma(x+1) for x in hist))/LN2


def half_log(logk, alpha, odds):
    return float(logsumexp((COUNTS + alpha*(logk[F+U]+logk[F+S-U])
                           + S*odds[0] + F*odds[1])*LN2)/LN2)


def regional(logk, alpha, groups, full, odds):
    singles, empty = groups-2*full, PAIRS-groups+full
    if min(singles, empty, full) < 0:
        raise ValueError('Infeasible pair histogram')
    return (256*half_log(logk, alpha, odds) - singles*odds[0] - full*odds[1]
            - log_multinomial((empty, singles, full)))


def bound(profile, delta, witness):
    groups, first, second = profile
    count = profile_count_upper(PAIRS, *profile)
    if count == 0:
        raise ValueError('Empty support profile')
    z, alpha, s0, f0, s1, f1 = witness
    logk = np.log2(quiet_float(z))
    regions = (regional(logk, alpha, groups, first, (s0, f0))
               + regional(logk, alpha, groups, second, (s1, f1)))
    limit = int(Fraction(delta)*N)
    return (log2(count) + alpha*(groups*LOG_DENSITY-log2(255)-limit*log2(z)+128)
            + 128*regions)


def optimize(profile, delta, maxiter):
    groups, first, second = profile
    def markers(full):
        empty, singles = PAIRS-groups+full, groups-2*full
        return (log2(max(singles, .00001)/max(empty, .00001)),
                log2(max(full, .00001)/max(empty, .00001)))
    best = None
    for z, alpha in ((.7, .2), (.96, .5), (.99, 1)):
        start = (z, alpha)+markers(first)+markers(second)
        result = minimize(lambda x: bound(profile, delta, x), start, method='L-BFGS-B',
            bounds=((.001, .999999), (.00001, 1), (-50, 50), (-50, 50), (-50, 50), (-50, 50)),
            options=dict(maxiter=maxiter, maxfun=maxiter*10, ftol=1e-13, gtol=1e-6))
        row = dict(profile=list(profile), distance=str(delta), log2_upper_screen=float(result.fun),
                   witness=list(map(float, result.x)), evaluations=int(result.nfev),
                   optimizer_success=bool(result.success), directed=False)
        if best is None or result.fun < best['log2_upper_screen']:
            best = row
    return best


def pairing_bound(pair, delta, witness):
    """Cauchy--Schwarz summand for all source supports, at one histogram.

    For support S, let R_i(S) be its regional quiet factor. The two
    matchings have the same marginal support count. Therefore
      sum_S R_0(S)**128 R_1(S)**128 <= sum_S R_0(S)**256.
    Group the right side by the even-matching full-pair count f. Its
    multiplicity is multinomial(4096;e,s,f)*2**s. This does not assume
    independent even/odd support geometries. A sum over every feasible
    (g,f), with the outer/projective/terminal factors, bounds the quiet
    bad event. For a fixed witness the log cap below is convex in (g,f).
    """
    groups, full = pair
    z, alpha, sx, fx = witness
    empty, singles = PAIRS-groups+full, groups-2*full
    if min(empty, singles, full) < 0:
        raise ValueError('Infeasible histogram')
    return (alpha*(groups*LOG_DENSITY-log2(255)-int(Fraction(delta)*N)*log2(z)+128)
            + singles - 255*log_multinomial((empty, singles, full))
            + 65536*half_log(np.log2(quiet_float(z)), alpha, (sx, fx))
            - 256*(singles*sx+full*fx))


def optimize_pair(pair, delta, maxiter):
    groups, full = pair
    empty, singles = PAIRS-groups+full, groups-2*full
    markers = (log2(max(singles, .00001)/max(empty, .00001)),
               log2(max(full, .00001)/max(empty, .00001)))
    best = None
    for z, alpha in ((.7, .2), (.96, .5), (.99, 1)):
        result = minimize(lambda x: pairing_bound(pair, delta, x), (z, alpha)+markers,
            method='L-BFGS-B', bounds=((.001, 1), (.00001, 1), (-50, 50), (-50, 50)),
            options=dict(maxiter=maxiter, maxfun=maxiter*10, ftol=1e-13, gtol=1e-6))
        row = dict(pair=list(pair), distance=str(delta), log2_upper_screen=float(result.fun),
                   witness=list(map(float, result.x)), evaluations=int(result.nfev),
                   optimizer_success=bool(result.success), directed=False)
        if best is None or result.fun < best['log2_upper_screen']:
            best = row
    return best


def polygon(glo, ghi, flo, fhi):
    """Clip with exact rational vertices, not floating intersections."""
    glo, ghi, flo, fhi = map(Fraction, (glo, ghi, flo, fhi))
    points = [(glo, flo), (ghi, flo), (ghi, fhi), (glo, fhi)]
    for a, b, c in ((0, -1, 0), (1, -1, -4096), (-1, 2, 0)):
        result = []
        if not points:
            break
        for left, right in zip(points, points[1:]+points[:1]):
            vl, vr = a*left[0]+b*left[1]+c, a*right[0]+b*right[1]+c
            if vl <= 0:
                result.append(left)
            if (vl < 0 < vr) or (vr < 0 < vl):
                ratio = vl/(vl-vr)
                result.append(tuple(x+ratio*(y-x) for x, y in zip(left, right)))
        points = result
    return list(dict.fromkeys(points))


def directed_vertices(vertices, delta, witness, precision=192):
    """Arb reevaluation; decimal witness strings define exact rationals."""
    from flint import arb, ctx
    ctx.prec = precision
    def rational(value):
        fraction = Fraction(str(value))
        return arb(fraction.numerator)/fraction.denominator
    z, alpha, sx, fx = map(rational, witness)
    ln2 = arb(2).log()
    quiet = [sum(arb(comb(l, c))*z**(l-c)*(1-z)**c/Q**(c+min(4, l-c))
                 for c in range(l+1)) for l in range(17)]
    coefficient = arb(0)
    for s, f, u in ROWS:
        s, f, u = int(s), int(f), int(u)
        count = comb(16, s)*comb(16-s, f)*comb(s, u)
        coefficient += (arb(count)/2**s * (quiet[f+u]*quiet[f+s-u])**alpha
                        * ((s*sx+f*fx)*ln2).exp())
    common = 65536*coefficient.log()/ln2
    density = 1280-32*arb(255).log()/ln2
    fixed = alpha*(-arb(255).log()/ln2-int(Fraction(delta)*N)*z.log()/ln2+128)
    values = []
    for groups, full in vertices:
        groups, full = rational(groups), rational(full)
        empty, singles = PAIRS-groups+full, groups-2*full
        norm = (arb(PAIRS+1).lgamma()-(empty+1).lgamma()
                -(singles+1).lgamma()-(full+1).lgamma())/ln2
        values.append(fixed+alpha*groups*density+singles-255*norm+common
                      -256*(singles*sx+full*fx))
    return values


def cover(delta, maxseconds=40, target=-100, precision=192, include_witnesses=False):
    """Bounded directed proposal; z/alpha are common within each g band.

    This does NOT combine independently optimized (g,f) fractional moments.
    A failed g band is split into disjoint integer bands and restarted.
    Only positive markers vary inside a successful band's f rectangles.
    All accepted polygons receive directed vertex reevaluation.
    """
    from flint import arb
    started = perf_counter()
    ends = [1, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 1536, 2048,
            2560, 3072, 3584, 3840, 4096, 4352, 4608, 5120, 6144, 8192]
    bands = [(1, 1)]+[(a+1, b) for a, b in zip(ends, ends[1:])]
    accepted, attempts = [], 0
    while bands and perf_counter()-started < maxseconds:
        glo, ghi = bands.pop(0)
        centre = (glo+ghi)/2
        witness = optimize_pair((centre, centre/2), delta, 40)['witness']
        z, alpha = witness[:2]
        # The strict receipt interface uses z<1, including sparse bands.
        z = min(z, 1-2**-40)
        work, leaves, failed = [(max(0, glo-PAIRS), min(PAIRS, ghi/2))], [], False
        while work:
            if perf_counter()-started >= maxseconds or len(leaves)+len(work) > 128:
                failed = True
                break
            flo, fhi = work.pop()
            vertices = polygon(glo, ghi, flo, fhi)
            if not vertices:
                continue
            groups = sum(x for x, _ in vertices)/len(vertices)
            full = sum(y for _, y in vertices)/len(vertices)
            empty, singles = PAIRS-groups+full, groups-2*full
            odds = (log2(max(singles, .00001)/max(empty, .00001)),
                    log2(max(full, .00001)/max(empty, .00001)))
            result = minimize(lambda p: pairing_bound((groups, full), delta, (z, alpha, *p)),
                odds, method='L-BFGS-B', bounds=((-50, 50), (-50, 50)),
                options=dict(maxiter=50, maxfun=180, ftol=1e-12))
            candidate = [z, alpha, *map(float, result.x)]
            attempts += 1
            floating = max(pairing_bound(v, delta, candidate) for v in vertices)
            if floating <= target-1e-5:
                directed = directed_vertices(vertices, delta, candidate, precision)
                if not all(value < arb(target) for value in directed):
                    raise ArithmeticError('Directed check did not prove target')
                leaves.append(dict(f_range=[flo, fhi], witness=candidate,
                                   vertex_upper=max(float(v.upper()) for v in directed)))
            elif fhi-flo > .25:
                mid = (flo+fhi)/2
                work.extend(((flo, mid), (mid, fhi)))
            else:
                failed = True
                break
        if failed:
            if glo == ghi or perf_counter()-started >= maxseconds:
                bands.insert(0, (glo, ghi))
                break
            mid = (glo+ghi)//2
            bands[:0] = [(glo, mid), (mid+1, ghi)]
        else:
            accepted.append(dict(g_range=[glo, ghi], mark=z, alpha=alpha,
                leaf_count=len(leaves), maximum_vertex_upper=max(row['vertex_upper'] for row in leaves),
                leaves=leaves))
    # Every feasible integer histogram is determined by (e,s,f) summing to
    # 4096. Removing the all-empty histogram gives binom(4098,2)-1 terms.
    terms = comb(PAIRS+2, 2)-1
    # Independent interval coverage check, on every integer group count.
    seen_groups = set()
    for band in accepted:
        glo, ghi = band['g_range']
        for groups in range(glo, ghi+1):
            if groups in seen_groups:
                raise ArithmeticError('Overlapping integer group bands')
            seen_groups.add(groups)
            low, high = max(0, groups-PAIRS), min(PAIRS, groups/2)
            reached = low
            for leaf in sorted(band['leaves'], key=lambda row: row['f_range']):
                a, b = leaf['f_range']
                a, b = max(a, low), min(b, high)
                if a > b:
                    continue
                if a > reached:
                    raise ArithmeticError('Gap in full-pair coverage')
                reached = max(reached, b)
            if reached < high:
                raise ArithmeticError('Incomplete full-pair coverage')
    if not bands and seen_groups != set(range(1, 2*PAIRS+1)):
        raise ArithmeticError('Incomplete integer group coverage')
    # Persist exact rational literals, not machine-dependent float objects.
    for band in accepted:
        band['mark'], band['alpha'] = (str(Fraction(str(band[name]))) for name in ('mark', 'alpha'))
        for leaf in band['leaves']:
            leaf['f_range'] = list(map(lambda x: str(Fraction(x)), leaf['f_range']))
            leaf['witness'] = list(map(lambda x: str(Fraction(str(x))), leaf['witness']))
    witness_hash = sha256(json.dumps(accepted, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
    for band in accepted:
        if not include_witnesses:
            del band['leaves']
    return dict(receipt_schema='cheap-lane-quiet-cover-v1',
                scope='all-quiet scalar subcode only', directed=True,
                distance=str(Fraction(delta)), forbidden_weight=int(Fraction(delta)*N),
                complete_g_coverage=not bands, target_log2_per_summand=target,
                integer_histograms=terms, log2_union_ceiling=target+log2(terms) if not bands else None,
                exact_union_numerator=terms if not bands else None,
                exact_union_denominator_power2=-target if not bands else None,
                guaranteed_margin_bits=-target-terms.bit_length() if not bands else None,
                precision_bits=precision, independently_checked_integer_groups=len(seen_groups),
                accepted_rectangle_count=sum(row['leaf_count'] for row in accepted),
                witness_sha256=witness_hash,
                source_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
                z_alpha_common_per_g_band=True, terminal_loss_bits='128*alpha',
                accepted_bands=accepted, pending_g_bands=bands, rectangle_attempts=attempts,
                elapsed_seconds=perf_counter()-started, writes=0,
                whole_code_certificate=False)


def receipt_fraction(value):
    """Read a bounded finite rational literal; never evaluate an expression."""
    if type(value) not in (str, int, float) or (type(value) is float and not isfinite(value)):
        raise ValueError('Expected a finite rational literal')
    if type(value) is str:
        if len(value) > 128 or not re.fullmatch(r'[+-]?\d+(?:/\d+|\.\d+)?', value):
            raise ValueError('Expected a bounded integer, decimal, or fraction literal')
    result = Fraction(str(value))
    if max(result.numerator.bit_length(), result.denominator.bit_length()) > 4096:
        raise ValueError('Rational literal is too large')
    return result


def validate_receipt(value):
    """Validate scope and reconstruct geometry; claimed bounds are ignored."""
    if (type(value) is not dict or value.get('receipt_schema') != 'cheap-lane-quiet-cover-v1'
            or value.get('scope') != 'all-quiet scalar subcode only'
            or value.get('whole_code_certificate') is not False
            or receipt_fraction(value.get('distance')) != Fraction(2, 5)
            or type(value.get('forbidden_weight')) is not int or value['forbidden_weight'] != 838860
            or type(value.get('target_log2_per_summand')) is not int
            or value['target_log2_per_summand'] != -100):
        raise ValueError('Receipt is not the fixed quiet-only 40% / -100-bit cover')
    bands = value.get('accepted_bands')
    if type(bands) is not list or not 1 <= len(bands) <= 512:
        raise ValueError('Require 1..512 group bands')
    result, total_leaves = [], 0
    for band in bands:
        if type(band) is not dict:
            raise ValueError('Invalid group band')
        interval = band.get('g_range')
        if (type(interval) is not list or len(interval) != 2
                or any(type(x) is not int for x in interval)
                or not 1 <= interval[0] <= interval[1] <= 8192):
            raise ValueError('Invalid integer group-count interval')
        z, alpha = map(receipt_fraction, (band.get('mark'), band.get('alpha')))
        if not 0 < z < 1 or not 0 < alpha <= 1:
            raise ValueError('Require 0<z<1 and 0<alpha<=1')
        leaves = band.get('leaves')
        if type(leaves) is not list or not leaves:
            raise ValueError('Missing leaf witnesses')
        total_leaves += len(leaves)
        if total_leaves > 512:
            raise ValueError('At most 512 leaf witnesses are allowed')
        checked = []
        for leaf in leaves:
            if type(leaf) is not dict:
                raise ValueError('Invalid leaf')
            frange, witness = leaf.get('f_range'), leaf.get('witness')
            if type(frange) is not list or len(frange) != 2 or type(witness) is not list or len(witness) != 4:
                raise ValueError('Require two full-pair endpoints and four witness entries')
            low, high = map(receipt_fraction, frange)
            witness = list(map(receipt_fraction, witness))
            if not 0 <= low <= high <= 4096 or witness[:2] != [z, alpha]:
                raise ValueError('Invalid full-pair interval or noncommon mark/exponent')
            # Both log markers must be representable finite numbers. The
            # exact rationals, not these float conversions, enter Arb.
            if any(not isfinite(float(x)) for x in witness[2:]):
                raise ValueError('Nonfinite marker')
            vertices = polygon(*interval, low, high)
            if not vertices:
                raise ValueError('Leaf has no feasible histogram')
            checked.append(dict(f_range=(low, high), witness=witness, vertices=vertices))
        result.append(dict(g_range=tuple(interval), leaves=checked))
    result.sort(key=lambda row: row['g_range'])
    previous = 0
    for band in result:
        glo, ghi = band['g_range']
        if glo != previous+1:
            raise ValueError('Overlapping or missing integer group-count band')
        previous = ghi
        for groups in range(glo, ghi+1):
            low, high = Fraction(max(0, groups-PAIRS)), min(Fraction(PAIRS), Fraction(groups, 2))
            reached = low
            for leaf in sorted(band['leaves'], key=lambda row: row['f_range']):
                a, b = leaf['f_range']
                a, b = max(a, low), min(b, high)
                if a > b:
                    continue
                if a > reached:
                    raise ValueError('Gap in rational full-pair coverage')
                reached = max(reached, b)
            if reached < high:
                raise ValueError('Incomplete rational full-pair coverage')
    if previous != 8192:
        raise ValueError('Missing terminal group counts')
    return result


def verify_value(value, precision=256, maxseconds=40):
    """Recompute all directed vertices without invoking any optimizer."""
    from flint import arb
    started = perf_counter()
    bands = validate_receipt(value)
    worst, vertices, leaves = float('-inf'), 0, 0
    for band in bands:
        for leaf in band['leaves']:
            if perf_counter()-started >= maxseconds:
                raise TimeoutError('Bounded receipt verification exceeded 40 seconds')
            bounds = directed_vertices(leaf['vertices'], '2/5', leaf['witness'], precision)
            if not all(bound.is_finite() and bound < arb(-100) for bound in bounds):
                raise ArithmeticError('A recomputed vertex does not prove the required -100-bit cap')
            worst = max(worst, *(float(bound.upper()) for bound in bounds))
            vertices += len(bounds)
            leaves += 1
    if perf_counter()-started >= maxseconds:
        raise TimeoutError('Bounded receipt verification exceeded 40 seconds')
    return dict(verification='passed', scope='all-quiet scalar subcode only',
        distance='2/5', forbidden_weight=838860, target_log2_per_summand=-100,
        independently_checked_integer_groups=8192, group_bands=len(bands),
        accepted_rectangle_count=leaves, recomputed_vertices=vertices,
        maximum_vertex_upper_diagnostic=worst, precision_bits=precision,
        exact_union_numerator=8394752, exact_union_denominator_power2=100,
        guaranteed_margin_bits=76, optimizer_invocations=0, writes=0,
        source_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
        elapsed_seconds=perf_counter()-started, whole_code_certificate=False)


def verify_receipt(path, precision=256):
    """Read at most 256 KiB and verify; hashes never replace arithmetic."""
    path = Path(path)
    if path.stat().st_size > 256*1024:
        raise ValueError('Receipt exceeds 256 KiB')
    with path.open('rb') as handle:
        data = handle.read(256*1024+1)
    if len(data) > 256*1024:
        raise ValueError('Receipt exceeds 256 KiB')
    def reject_constant(value):
        raise ValueError('Nonfinite JSON constant: '+value)
    value = json.loads(data, parse_constant=reject_constant)
    result = verify_value(value, precision)
    result['receipt_sha256'] = sha256(data).hexdigest()
    return result


def checks():
    for z in (Fraction(1, 4), Fraction(1, 2), Fraction(9, 10), Fraction(1)):
        for l in range(17):
            value = quiet_exact(l, z)
            unconstrained = (z+(1-z)/Q)**l
            assert unconstrained/Q**4 <= value <= unconstrained
            assert value == Fraction(1, Q**l) if l <= 4 else value > 0
            assert abs(float(value)-quiet_float(float(z))[l]) < 1e-14
    # alpha=0 must recover the multinomial polynomial (1+x+y)^16.
    for odds in ((0., 0.), (-3., -1.), (2., -2.)):
        got = half_log(np.zeros(17), .5, odds)
        expected = 16*log2(1+2**odds[0]+2**odds[1])
        assert abs(got-expected) < 1e-12
    assert set(polygon(1, 8192, 0, 4096)) == {
        (Fraction(1), Fraction(0)), (Fraction(4096), Fraction(0)),
        (Fraction(8192), Fraction(4096)), (Fraction(1), Fraction(1, 2))}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--profile', nargs=3, type=int, action='append')
    parser.add_argument('--pair', nargs=2, type=int, action='append')
    parser.add_argument('--distance', default='2/5')
    parser.add_argument('--maxiter', type=int, default=80)
    parser.add_argument('--cover', action='store_true')
    parser.add_argument('--precision', type=int, default=192)
    parser.add_argument('--include-witnesses', action='store_true')
    parser.add_argument('--verify-receipt', type=Path)
    args = parser.parse_args()
    if not 1 <= args.maxiter <= 150:
        parser.error('Require 1 <= maxiter <= 150')
    if len(args.profile or [])+len(args.pair or []) > 12:
        parser.error('At most twelve profiles per bounded invocation')
    if not 128 <= args.precision <= 512:
        parser.error('Require 128 <= precision <= 512')
    if args.verify_receipt and (args.cover or args.profile or args.pair or args.include_witnesses):
        parser.error('Receipt verification cannot be combined with proposal modes')
    if args.include_witnesses and not args.cover:
        parser.error('--include-witnesses requires --cover')
    if args.include_witnesses and (args.profile or args.pair):
        parser.error('Standalone receipt output cannot be combined with profile screens')
    if args.verify_receipt:
        print(json.dumps(verify_receipt(args.verify_receipt, args.precision)), flush=True)
        raise SystemExit(0)
    started = perf_counter()
    checks()
    if args.cover and args.include_witnesses:
        print(json.dumps(cover(args.distance, precision=args.precision,
                               include_witnesses=True)), flush=True)
        raise SystemExit(0)
    print(json.dumps(dict(checks='pass', scope='quiet-only profile FLOAT screen',
                          terminal_loss_bits='128*alpha', writes=0)), flush=True)
    for profile in args.profile or []:
        print(json.dumps(optimize(tuple(profile), args.distance, args.maxiter)), flush=True)
    for pair in args.pair or []:
        print(json.dumps(optimize_pair(tuple(pair), args.distance, args.maxiter)), flush=True)
    if args.cover:
        print(json.dumps(cover(args.distance, precision=args.precision,
                               include_witnesses=args.include_witnesses)), flush=True)
    print(json.dumps(dict(elapsed_seconds=perf_counter()-started)), flush=True)
