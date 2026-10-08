"""Bounded quiet proof experiment retaining the global outer support floor.

The actual randomized MDS[256,128] outer has input weight at least 129r for
r active rows. A positive input mark enforces this floor after the artificial
uniform comparison. Arb subtraction preserves the small nonzero-input terms.
"""

from fractions import Fraction as F
from functools import lru_cache
from math import comb, exp2, log2
import json

from hd_randomized_outer_joint_quiet import joint_moment_enclosures
from hd_randomized_outer_proof import Parameters
from hd_randomized_outer_quiet_exact import regional_coefficients, weighted_logsum_mean


def _arb_fraction(value):
    from flint import arb
    value = F(value)
    return arb(value.numerator)/value.denominator


@lru_cache(None)
def base_nonzero_logs(mark, field=65536, precision=384):
    """Arb subtraction followed by float witness data, never a certificate."""
    from flint import arb, ctx
    previous = ctx.prec
    ctx.prec = precision
    try:
        quiet, free = joint_moment_enclosures(F(mark), field)
        result = []
        for moments, quiet_kind in ((quiet, True), (free, False)):
            row = []
            for j, value in enumerate(moments):
                if j == 0 or (quiet_kind and j <= 4):
                    row.append(float('-inf'))
                    continue
                nonzero = value-arb(1)/field**j
                if not nonzero > 0:
                    raise ArithmeticError('Nonzero-input mass unresolved; increase subtraction precision')
                row.append(float(nonzero.log()/arb(2).log()))
            result.append(tuple(row))
        return tuple(result)
    finally:
        ctx.prec = previous


def weighted_moment_logs(mark, log_input, field=65536):
    if log_input < 0:
        raise ValueError('The input mark must be at least one')
    nonzero = base_nonzero_logs(F(mark), field)
    result = []
    for moments in nonzero:
        row = []
        for j, log_nonzero in enumerate(moments):
            zero = -j*log2(field)
            other = j*float(log_input)+log_nonzero
            maximum, minimum = max(zero, other), min(zero, other)
            row.append(maximum+log2(1+exp2(minimum-maximum)))
        result.append(tuple(row))
    return tuple(result)


def weighted_moment_enclosures(mark, log_input, field=65536):
    """K_j(v,z)<=q^-j+v^j(Khat_j(z)-q^-j), and similarly for T_j."""
    from flint import arb
    log_input = F(log_input)
    if log_input < 0:
        raise ValueError('The input mark must be at least one')
    quiet, free = joint_moment_enclosures(F(mark), field)
    v = arb(2)**_arb_fraction(log_input)
    result = []
    for moments, quiet_kind in ((quiet, True), (free, False)):
        row = []
        for j, value in enumerate(moments):
            zero = arb(1)/field**j
            if j == 0 or (quiet_kind and j <= 4):
                row.append(zero)
                continue
            nonzero = value-zero
            if not nonzero > 0:
                raise ArithmeticError('Nonzero-input mass unresolved; increase Arb precision')
            row.append(zero+v**j*nonzero)
        result.append(tuple(row))
    return tuple(result)


def optimized_marker(active_rows, mark, log_input, alpha, parameters=Parameters()):
    p = parameters
    if not 0 < active_rows < p.rows or not 0 < alpha <= 1:
        raise ValueError('Need interior occupancy and 0<alpha<=1')
    quiet, free = weighted_moment_logs(mark, log_input, p.field)
    choose_local = tuple(log2(comb(64, j)) for j in range(65))
    a = float(alpha)
    ordinary = tuple(c+a*k for c, k in zip(choose_local, quiet))
    terminal = tuple(c+a*k for c, k in zip(choose_local, free))
    lower, upper = -128.0, 128.0
    for unused in range(80):
        middle = (lower+upper)/2
        _, mean = weighted_logsum_mean(ordinary, middle)
        _, final_mean = weighted_logsum_mean(terminal, middle)
        if (p.packets-1)*mean+final_mean < p.outer_n*active_rows:
            lower = middle
        else:
            upper = middle
    log_odds = (lower+upper)/2
    log_ordinary, _ = weighted_logsum_mean(ordinary, log_odds)
    log_terminal, _ = weighted_logsum_mean(terminal, log_odds)
    choose = log2(comb(p.rows, active_rows))
    beta_log = p.outer_n*log2(p.field)-(p.outer_n-p.outer_k)*log2(p.field-1)
    floor = (p.outer_n-p.outer_k+1)*active_rows
    value = (a*(active_rows*beta_log-log2(p.field-1)-p.limit*log2(float(mark))
                -floor*float(log_input))
             +(p.packets-1)*log_ordinary+log_terminal
             -p.outer_n*active_rows*log_odds-(p.outer_n-1)*choose)
    marker = F(1/(1+exp2(-log_odds))).limit_denominator(1 << 40)
    if not 0 < marker < 1:
        raise ArithmeticError('Saddle marker rounded to a boundary')
    return dict(active_rows=active_rows, mark=str(F(mark)), log2_input_mark=str(F(log_input)),
                alpha=str(F(alpha)), marker=str(marker), marker_screen_log2=value,
                directed=False)


def screen_selected(occupancies=(1024, 2048, 4096)):
    marks = tuple(F(1, 1 << shift) for shift in (0, 2, 4, 8, 12, 14, 16))
    input_logs = tuple(map(F, (0, F(1, 2), 1, 2, 4, 8, 12, 16)))
    alphas = tuple(F(n, 32) for n in (1, 3, 8, 12, 16, 32))
    best = {r: None for r in occupancies}
    for mark in marks:
        for log_input in input_logs:
            for alpha in alphas:
                for active in occupancies:
                    witness = optimized_marker(active, mark, log_input, alpha)
                    if best[active] is None or witness['marker_screen_log2'] < best[active]['marker_screen_log2']:
                        best[active] = witness
    return tuple(best[r] for r in occupancies)


def directed_endpoint(active_rows, mark, log_input, alpha, marker=None, precision=384):
    from flint import arb, ctx
    previous = ctx.prec
    ctx.prec = precision
    try:
        p = Parameters()
        mark, log_input, alpha = map(F, (mark, log_input, alpha))
        proposal = optimized_marker(active_rows, mark, log_input, alpha)
        marker = F(proposal['marker']) if marker is None else F(marker)
        quiet, free = weighted_moment_enclosures(mark, log_input)
        region = regional_coefficients(quiet, free, alpha, marker, 128, active_rows)
        logarb = lambda value: value.log()/arb(2).log()
        choose = logarb(arb(comb(p.rows, active_rows)))
        beta_log = p.outer_n*logarb(arb(p.field))-(p.outer_n-p.outer_k)*logarb(arb(p.field-1))
        prefactor = choose+_arb_fraction(alpha)*(active_rows*beta_log-logarb(arb(p.field-1))
            -p.limit*logarb(_arb_fraction(mark))-129*active_rows*_arb_fraction(log_input))
        exact = prefactor+255*logarb(region['quiet'])+logarb(region['last'])
        marker_bound = (prefactor+(p.packets-1)*logarb(region['ordinary'])
                        +logarb(region['terminal'])-256*logarb(region['binomial_mass']))
        return dict(active_rows=active_rows, mark=str(mark), log2_input_mark=str(log_input),
                    alpha=str(alpha), marker=str(marker), input_support_floor=129*active_rows,
                    marker_log2_upper=str(marker_bound), exact_coefficient_log2_upper=str(exact),
                    below_2_to_minus_40=bool(exact < -40), directed=True,
                    covers_only_this_occupancy=True, whole_code_certificate=False)
    finally:
        ctx.prec = previous


def directed_input_derivative(active_rows, mark, alpha, precision=384):
    """Exact-coefficient derivative at log2(v)=0 for the proposed surrogate.

    The log bound is convex in log(v). A positive derivative certifies that
    this input mark cannot improve the specified fixed z,alpha witness.
    """
    from flint import arb, arb_poly, ctx
    previous = ctx.prec
    ctx.prec = precision
    try:
        mark, alpha = F(mark), F(alpha)
        proposal = optimized_marker(active_rows, mark, F(0), alpha)
        marker = F(proposal['marker'])
        m, aa = _arb_fraction(marker), _arb_fraction(alpha)
        quiet, free = joint_moment_enclosures(mark)
        polynomials, score_polynomials = [], []
        for moments, quiet_kind in ((quiet, True), (free, False)):
            coefficients, scored = [], []
            for j, value in enumerate(moments):
                coefficient = comb(64, j)*m**j*(1-m)**(64-j)*value**aa
                if j == 0 or (quiet_kind and j <= 4):
                    score = arb(0)
                else:
                    nonzero = value-arb(1)/65536**j
                    if not nonzero > 0:
                        raise ArithmeticError('Nonzero-input mass unresolved in derivative')
                    score = j*nonzero/value
                coefficients.append(coefficient)
                scored.append(coefficient*score)
            normalization = sum(coefficients, arb(0))
            polynomials.append(arb_poly([x/normalization for x in coefficients]))
            score_polynomials.append(arb_poly([x/normalization for x in scored]))
        ordinary, terminal = polynomials
        ordinary_score, terminal_score = score_polynomials
        prefix = ordinary**126
        prefix127 = prefix*ordinary
        ordinary_mass = (prefix127*ordinary)[active_rows]
        terminal_mass = (prefix127*terminal)[active_rows]
        if not ordinary_mass > 0 or not terminal_mass > 0:
            raise ArithmeticError('Coefficient lost positivity in derivative')
        score_region = 128*(prefix127*ordinary_score)[active_rows]/ordinary_mass
        score_last = ((127*prefix*ordinary_score*terminal+prefix127*terminal_score)[active_rows]
                      /terminal_mass)
        derivative = aa*(255*score_region+score_last-129*active_rows)
        return dict(active_rows=active_rows, mark=str(mark), alpha=str(alpha),
                    derivative_at_log2_v_zero=str(derivative),
                    no_input_mark_improves_this_witness=bool(derivative >= 0),
                    justification='convexity in log(input mark)', fixed_witness_only=True,
                    directed=True, whole_code_certificate=False)
    finally:
        ctx.prec = previous


if __name__ == '__main__':
    for proposal in screen_selected():
        print(json.dumps(directed_endpoint(proposal['active_rows'], proposal['mark'],
                                            proposal['log2_input_mark'], proposal['alpha'],
                                            proposal['marker'])), flush=True)
        print(json.dumps(directed_input_derivative(proposal['active_rows'], proposal['mark'],
                                                  proposal['alpha'])), flush=True)
