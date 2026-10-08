"""Quiet fractional moments conditioned on the eight tile occupancies.

The default run compares two fixed middle witnesses. The original quiet
residual rank rule is retained. A separate control isolates the earlier
Jensen step from applying the unweighted syndrome cap before averaging.
"""

from fractions import Fraction as F
from functools import lru_cache
from math import comb, exp2, factorial, log2

from hd_randomized_outer_quiet_exact import (
    _arb_fraction, optimized_marker, regional_coefficients, weighted_logsum_mean,
)
from hd_randomized_outer_state_screen import NEG, add


def join(previous, local, distinct_hyperplanes=False):
    if distinct_hyperplanes and previous == local == 3:
        return 4
    return 2 if previous == local == 1 else max(previous, local)


@lru_cache(None)
def local_rank_coefficients(mark, field=65536):
    """Scaled exact conditional-tile coefficients; no binom(8,a) factor."""
    mark = F(mark)
    if not 0 < mark <= 1:
        raise ValueError('Need 0<z<=1')
    result = [[F(0) for unused in range(5)] for unused in range(9)]
    for a in range(9):
        for v in range(9):
            residual = max(a-v, 0)
            result[a][min(4, residual)] += (comb(8, v)*mark**(8-v)*(1-mark)**v
                                          *field**residual)
    return tuple(tuple(row) for row in result)


@lru_cache(None)
def profile_moments(mark, field=65536, distinct_hyperplanes=False):
    """Exact rational upper K and exact T for every sorted tile profile.

    Integer prefix coefficients share denominator mark.denominator**64.
    This avoids cancellation when subtracting the known all-zero input atom.
    """
    mark = F(mark)
    if not 0 < mark <= 1 or field < 2:
        raise ValueError('Need 0<z<=1 and field>=2')
    numerator, denominator = mark.numerator, mark.denominator
    local = [[0]*5 for unused in range(9)]
    for a in range(9):
        for v in range(9):
            residual = max(a-v, 0)
            local[a][min(4, residual)] += (comb(8, v)*numerator**(8-v)
                                         *(denominator-numerator)**v*field**residual)
    tile_free = tuple(sum(row) for row in local)
    common = denominator**64
    output, counts = [], [0]*9
    rank_weights = tuple(field**(4-r) for r in range(5))

    def visit(profile, start, active, multiplicity, coefficients, free_numerator):
        if len(profile) == 8:
            duplicate_factor = 1
            for count in counts:
                duplicate_factor *= factorial(count)
            ordered_multiplicity = factorial(8)//duplicate_factor*multiplicity
            raw = F(sum(x*y for x, y in zip(coefficients, rank_weights)),
                    common*field**(active+4))
            quiet = min(raw, F(1, field**min(active, 4)))
            free = F(free_numerator, common*field**active)
            zero = F(1, field**active)
            if not zero <= quiet <= free <= 1:
                raise ArithmeticError('Invalid profile moment or zero atom')
            output.append((profile, active, ordered_multiplicity, quiet, free))
            return
        for a in range(start, 9):
            following = [0]*5
            for previous, value in enumerate(coefficients):
                for residual, coefficient in enumerate(local[a]):
                    following[join(previous, residual, distinct_hyperplanes)] += value*coefficient
            counts[a] += 1
            visit(profile+(a,), a, active+a, multiplicity*comb(8, a), following,
                  free_numerator*tile_free[a])
            counts[a] -= 1

    visit((), 0, 0, 1, (1, 0, 0, 0, 0), 1)
    return tuple(output)


def profile_logs(mark, field=65536):
    """Enumerate all12870 sorted tile profiles with their exact multiplicities.

    Each entry is (profile,j,integer_multiplicity,log2 Khat_profile,log2 T_profile).
    Prefix sharing evaluates only24309 nonempty prefixes of sorted profiles.
    """
    logarithm = lambda value: log2(value.numerator)-log2(value.denominator)
    return tuple((profile, active, multiplicity, logarithm(quiet), logarithm(free))
                 for profile, active, multiplicity, quiet, free in profile_moments(F(mark), field))


def fractional_profile_logs(mark, alpha, field=65536, input_mark=F(1)):
    """log E[Khat_profile**alpha], log E[T_profile**alpha], by packet size."""
    alpha, input_mark = F(alpha), F(input_mark)
    if not 0 < alpha <= 1 or input_mark < 1:
        raise ValueError('Need 0<alpha<=1 and input mark>=1')
    quiet, free = [NEG]*65, [NEG]*65
    mass = [0]*65
    profiles = profile_moments(F(mark), field)
    for profile, active, multiplicity, raw_quiet, raw_free in profiles:
        zero = F(1, field**active)
        weighted_quiet = zero+input_mark**active*(raw_quiet-zero)
        weighted_free = zero+input_mark**active*(raw_free-zero)
        log_quiet = log2(weighted_quiet.numerator)-log2(weighted_quiet.denominator)
        log_free = log2(weighted_free.numerator)-log2(weighted_free.denominator)
        weight = log2(multiplicity)
        quiet[active] = add(quiet[active], weight+float(alpha)*log_quiet)
        free[active] = add(free[active], weight+float(alpha)*log_free)
        mass[active] += multiplicity
    for active in range(65):
        if mass[active] != comb(64, active):
            raise ArithmeticError('Tile-profile multiplicities do not normalize')
        normalizer = log2(comb(64, active))
        quiet[active] -= normalizer
        free[active] -= normalizer
    return tuple(quiet), tuple(free)


def fractional_profile_enclosures(mark, alpha, field=65536, input_mark=F(1)):
    """Raw Arb bounds E_profile[Khat**alpha], E_profile[T**alpha].

    Both arrays are ALREADY fractionated. Optional input support marking
    removes the exact zero atom in rational arithmetic before conversion.
    """
    from flint import arb
    alpha, input_mark = F(alpha), F(input_mark)
    if not 0 < alpha <= 1 or input_mark < 1:
        raise ValueError('Need 0<alpha<=1 and input mark>=1')
    a = _arb_fraction(alpha)
    quiet, free = [arb(0) for unused in range(65)], [arb(0) for unused in range(65)]
    mass = [0]*65
    for profile, active, multiplicity, k, t in profile_moments(F(mark), field):
        zero = F(1, field**active)
        k = zero+input_mark**active*(k-zero)
        t = zero+input_mark**active*(t-zero)
        quiet[active] += multiplicity*_arb_fraction(k)**a
        free[active] += multiplicity*_arb_fraction(t)**a
        mass[active] += multiplicity
    for active in range(65):
        if mass[active] != comb(64, active):
            raise ArithmeticError('Tile-profile multiplicities do not normalize')
        quiet[active] /= mass[active]
        free[active] /= mass[active]
    return tuple(quiet), tuple(free)


def marker_screen(active_rows, mark, alpha, quiet, free, input_mark=F(1)):
    """Optimized marker for already-fractional local moments; FLOAT only."""
    mark, alpha, input_mark = F(mark), F(alpha), F(input_mark)
    combinations = tuple(log2(comb(64, j)) for j in range(65))
    ordinary = tuple(c+k for c, k in zip(combinations, quiet))
    terminal = tuple(c+k for c, k in zip(combinations, free))
    lower, upper = -128.0, 128.0
    for unused in range(80):
        middle = (lower+upper)/2
        _, mean = weighted_logsum_mean(ordinary, middle)
        _, last_mean = weighted_logsum_mean(terminal, middle)
        if 32767*mean+last_mean < 256*active_rows:
            lower = middle
        else:
            upper = middle
    odds = (lower+upper)/2
    first, _ = weighted_logsum_mean(ordinary, odds)
    last, _ = weighted_logsum_mean(terminal, odds)
    beta_log = 4096-128*log2(65535)
    bound = (float(alpha)*(active_rows*beta_log-log2(65535)-838860*log2(mark))
             +32767*first+last-256*active_rows*odds-255*log2(comb(8192, active_rows))
             -129*active_rows*float(alpha)*log2(input_mark))
    return dict(active_rows=active_rows, mark=str(mark), alpha=str(alpha),
                input_mark=str(input_mark),
                marker_log_odds=odds, log2_upper_screen=bound,
                directed=False, whole_code_certificate=False)


def compare_fixed(active_rows, mark, alpha):
    quiet, free = fractional_profile_logs(mark, alpha)
    result = marker_screen(active_rows, mark, alpha, quiet, free)
    original = optimized_marker(active_rows, mark, alpha)
    capped_quiet, capped_free = fractional_profile_logs(mark, F(1))
    cap_control = marker_screen(active_rows, mark, alpha,
                               tuple(float(alpha)*x for x in capped_quiet),
                               tuple(float(alpha)*x for x in capped_free))
    result['old_optimized_marker_log2'] = original['marker_screen_log2']
    result['gain_bits'] = original['marker_screen_log2']-result['log2_upper_screen']
    result['early_cap_control_log2'] = cap_control['log2_upper_screen']
    result['early_cap_gain_bits'] = original['marker_screen_log2']-cap_control['log2_upper_screen']
    result['conditioning_gain_bits'] = cap_control['log2_upper_screen']-result['log2_upper_screen']
    result['profiles'] = 12870
    return result


def directed_endpoint(active_rows, mark, alpha, input_mark=F(1), exact_regional=False):
    """Directed replay at the floating optimized marker, optionally exact regional."""
    from flint import arb
    mark, alpha, input_mark = F(mark), F(alpha), F(input_mark)
    logs = fractional_profile_logs(mark, alpha, input_mark=input_mark)
    proposal = marker_screen(active_rows, mark, alpha, *logs, input_mark=input_mark)
    marker = F(1/(1+exp2(-proposal['marker_log_odds']))).limit_denominator(1 << 40)
    quiet, free = fractional_profile_enclosures(mark, alpha, input_mark=input_mark)
    m, a = _arb_fraction(marker), _arb_fraction(alpha)
    logarb = lambda value: value.log()/arb(2).log()
    ordinary = sum((comb(64, j)*m**j*(1-m)**(64-j)*quiet[j] for j in range(65)), arb(0))
    terminal = sum((comb(64, j)*m**j*(1-m)**(64-j)*free[j] for j in range(65)), arb(0))
    choose = logarb(arb(comb(8192, active_rows)))
    beta = 4096-128*logarb(arb(65535))
    prefix = (choose+a*(active_rows*beta-logarb(arb(65535))-838860*logarb(_arb_fraction(mark))
                       -129*active_rows*logarb(_arb_fraction(input_mark))))
    mass = arb(comb(8192, active_rows))*m**active_rows*(1-m)**(8192-active_rows)
    bound = prefix+32767*logarb(ordinary)+logarb(terminal)-256*logarb(mass)
    result = dict(active_rows=active_rows, mark=str(mark), alpha=str(alpha),
                  input_mark=str(input_mark), marker=str(marker), marker_log2_upper=str(bound),
                  directed=True, whole_code_certificate=False)
    if exact_regional:
        region = regional_coefficients(quiet, free, F(1), marker, 128, active_rows)
        exact = prefix+255*logarb(region['quiet'])+logarb(region['last'])
        result['exact_regional_log2_upper'] = str(exact)
    return result


if __name__ == '__main__':
    import json
    for rows, mark in ((1536, F(1, 4)), (4096, F(1, 4096))):
        print(json.dumps(compare_fixed(rows, mark, F(1, 32)), indent=2), flush=True)
