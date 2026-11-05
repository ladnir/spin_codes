"""General positive quartic weights; bounded(2,1,4) temporal-capacity control.

The encoder is unchanged. Network capacities(source,temporal,input) only
change the proof comparison. For old W=w,E=u-w,Z=4-u and non-H counts
r,s,t, put j=r+s+t, gain=c*j-(a-b)*(r+s)-b*r,
 loss=max(0,4*a+b*u-gain), cost=(a-b)*(u-r-s)+b*(w-r).
The local term is C(w,r)C(u-w,s)C(4-u,t)*x^j*y^cost
 times[p+g*alpha^4*beta^u]. Forced zero-gain/loss cases are explicit.

Only this15-coefficient table changes. The validated quartic interface
uses the existing positive power engine and old multinomial division.
Research envelope callbacks are explicit; no globals are monkeypatched.
The endpoint is(1+x)^(n-ell)*(y^(a-b)+x)^(ell-w)*(y^a+x)^w.

The sole production control is(a,b,c)=(2,1,4), y=3/2, threshold184116,
endpoint-matched central gauge(15/13,6/5), the same nine+-0.5 modes,
and quota x=1,9/4,81/16. This is a FLOAT geometry-tail diagnostic,
not the required tail proof or a modification of the construction.
No jobs or large arrays occur on import. Parent owns execution.
"""

import argparse
from fractions import Fraction
from functools import lru_cache
from itertools import combinations, product
from math import comb
from operator import index
import os
import shutil
import threading
import time

for _pt_thread_var in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[_pt_thread_var] = '1'

import numpy as np

PT_LD = np.longdouble
PT_TERMS = tuple((w, e) for w in range(5) for e in range(5-w))


def _pt_dependency(name):
    if name in globals():
        return globals()[name]
    if name.startswith('_ae_') or name.startswith('adaptive_envelope'):
        import quartet_adaptive_envelope as module
    elif name.startswith('_ac_') or name == 'adaptive_cover_local_coefficients':
        import quartet_adaptive_cover_checks as module
    else:
        import quartet_adaptive_population as module
    return getattr(module, name)


def _pt_capacities(capacities):
    a, b, c = map(index, capacities)
    if a < 1 or not 0 <= b <= a or c < a+b:
        raise ValueError('Require source>=1,0<=temporal<=source,input>=source+temporal')
    return a, b, c


@lru_cache(maxsize=None)
def general_flow_weights(y, capacities, u, w, j, r_A, r_W):
    y = Fraction(y)
    a, b, c = _pt_capacities(capacities)
    if (y <= 1 or not 0 <= w <= u <= 4 or not 0 <= r_W <= w
            or not 0 <= r_A-r_W <= u-w or not 0 <= j-r_A <= 4-u):
        raise ValueError('Invalid local labels or marker')
    gain = c*j-(a-b)*r_A-b*r_W
    loss = max(0, 4*a+b*u-gain)
    if j == 0 or gain == 0:
        return Fraction(1), Fraction(0), Fraction(1), loss, gain
    if loss == 0:
        return Fraction(0), Fraction(1), Fraction(0), loss, gain
    A, B = y**(-loss), y**(-gain)
    absent, selected = (1-A)/(1-A*B), (1-B)/(1-A*B)
    return B*absent, selected, absent, loss, gain


def _pt_probability_logs(log_y, loss, gain):
    if gain == 0:
        return PT_LD(0), PT_LD('-inf')
    if loss == 0:
        return PT_LD('-inf'), PT_LD(0)
    ln_y = log_y*np.log(PT_LD(2))
    denominator = -np.expm1(-(loss+gain)*ln_y)
    first, second = -np.expm1(-loss*ln_y), -np.expm1(-gain*ln_y)
    if min(denominator, first, second) <= 0 or not np.isfinite(denominator+first+second):
        raise ArithmeticError('General adaptive probability range failure')
    return (-gain*log_y+np.log2(first)-np.log2(denominator),
            np.log2(second)-np.log2(denominator))


def general_flow_quartic(y, x, alpha, beta, capacities=(2, 1, 4)):
    a, b, c = _pt_capacities(capacities)
    y, x, alpha, beta = map(PT_LD, (y, x, alpha, beta))
    if y <= 1 or min(x, alpha, beta) <= 0 or not all(np.isfinite(v) for v in (y, x, alpha, beta)):
        raise ValueError('Require finite y>1 and positive x,alpha,beta')
    log_y, log_x, log_alpha, log_beta = map(np.log2, (y, x, alpha, beta))
    local_logs = []
    for w, e in PT_TERMS:
        u = w+e
        terms = []
        for r in range(w+1):
            for s in range(e+1):
                for t in range(5-u):
                    j = r+s+t
                    gain = c*j-(a-b)*(r+s)-b*r
                    loss = max(0, 4*a+b*u-gain)
                    cost = (a-b)*(u-r-s)+b*(w-r)
                    log_g, log_p = _pt_probability_logs(log_y, loss, gain)
                    common = (np.log2(PT_LD(comb(w, r)*comb(e, s)*comb(4-u, t)))
                              +j*log_x+cost*log_y)
                    if np.isfinite(log_p):
                        terms.append(common+log_p)
                    if np.isfinite(log_g):
                        terms.append(common+log_g+4*log_alpha+u*log_beta)
        local_logs.append(np.log2(PT_LD(comb(4, u)*comb(u, w)))+_pt_dependency('_ap_logadd')(terms))
    logs = np.asarray(local_logs, dtype=PT_LD)
    scale = logs.max()
    coefficients = np.exp2(logs-scale)
    if np.any(coefficients <= 0) or not np.all(np.isfinite(coefficients)):
        raise ArithmeticError('General positive quartic normalization failed')
    return dict(coefficients=coefficients, scale=scale, terms=PT_TERMS,
                log2_coefficients=logs, y=y, x=x, alpha=alpha, beta=beta,
                log2_y=log_y, log2_x=log_x, log2_alpha=log_alpha, log2_beta=log_beta,
                capacities=(a, b, c))


def general_flow_endpoint_logs(prepared, y, x, capacities=(2, 1, 4)):
    source, temporal, _ = _pt_capacities(capacities)
    y, x = PT_LD(y), PT_LD(x)
    ell, omega, width = prepared['a'], prepared['w'], prepared['width']
    return ((width-ell)*np.log2(1+x)+(ell-omega)*np.log2(y**(source-temporal)+x)
            +omega*np.log2(y**source+x))


def general_flow_population_action(y, x, alpha, beta, width=512, progress=False, capacities=(2, 1, 4)):
    quartic = general_flow_quartic(y, x, alpha, beta, capacities)
    result = _pt_dependency('adaptive_population_action')(y, x, alpha, beta, width, progress, quartic)
    result['network_capacities'] = _pt_capacities(capacities)
    if quartic['alpha'] == quartic['beta'] == 1:
        expected = general_flow_endpoint_logs(result, y, x, capacities)
        error = np.max(np.abs(result['log2_action']-expected))
        if error > PT_LD('1e-8'):
            raise ArithmeticError('General all-row unit-potential identity failed')
        result['unit_potential_all_row_log2_error'] = error
    return result


def partial_temporal_action(y, x, alpha, beta, width=512):
    return general_flow_population_action(y, x, alpha, beta, width, capacities=(2, 1, 4))


def partial_temporal_endpoint(prepared, y, x):
    return general_flow_endpoint_logs(prepared, y, x, (2, 1, 4))


def partial_temporal_control(iterations=2, progress=False):
    if not 1 <= index(iterations) <= 2:
        raise ValueError('Use one or two bounded witnesses')
    started = time.monotonic()
    prepared = _pt_dependency('adaptive_envelope_prepare')(
        512, PT_LD(15)/13, PT_LD(6)/5, (-1/2, 0, 1/2))
    cache = _pt_dependency('adaptive_envelope_cache')(
        prepared, PT_LD(3)/2, (1, PT_LD(9)/4, PT_LD(81)/16), 16, 4096, progress,
        action_callback=partial_temporal_action, endpoint_callback=partial_temporal_endpoint)
    result = _pt_dependency('adaptive_envelope_cached_bound')(cache, 184116, 256, iterations, progress)
    result.update(network_capacities=(2, 1, 4), central_gauge=('15/13', '6/5'),
                  all_death_geometry_target=True, geometry_tail_only=True,
                  encoder_changed=False, seconds=time.monotonic()-started)
    return result


def _pt_exact_local(y, capacities, u, w, alpha, beta):
    source, temporal, _ = _pt_capacities(capacities)
    coefficients = [Fraction(0)]*5
    for r in range(w+1):
        for s in range(u-w+1):
            for t in range(5-u):
                j = r+s+t
                g, p, _, _, _ = general_flow_weights(y, capacities, u, w, j, r+s, r)
                cost = (source-temporal)*(u-r-s)+temporal*(w-r)
                coefficients[j] += (comb(w, r)*comb(u-w, s)*comb(4-u, t)*y**cost
                                    *(p+g*alpha**4*beta**u))
    return coefficients


@lru_cache(maxsize=None)
def _pt_physical_quartet(y, x, alpha, beta, capacities, A, W):
    source, temporal, _ = capacities
    value = Fraction(0)
    for J in range(16):
        g, p, _, _, _ = general_flow_weights(y, capacities, A.bit_count(), W.bit_count(),
                                            J.bit_count(), (J & A).bit_count(), (J & W).bit_count())
        cost = (source-temporal)*((15 ^ J) & A).bit_count()+temporal*((15 ^ J) & W).bit_count()
        for V, probability in ((0, p), (15, g)):
            value += probability*x**J.bit_count()*y**cost*alpha**V.bit_count()*beta**(V & A).bit_count()
    return value


def _pt_physical_action(width, y, x, alpha, beta, capacities, ell, omega):
    A, W = (1 << ell)-1, (1 << omega)-1
    partitions = _pt_dependency('_ap_partitions')(width)
    total = Fraction(0)
    for partition in partitions:
        value = Fraction(1)
        for quartet in partition:
            local_A = sum(((A >> coordinate) & 1) << slot for slot, coordinate in enumerate(quartet))
            local_W = sum(((W >> coordinate) & 1) << slot for slot, coordinate in enumerate(quartet))
            value *= _pt_physical_quartet(y, x, alpha, beta, capacities, local_A, local_W)
        total += value
    return total/len(partitions)


def _pt_copy_score(C, lookup, capacities):
    source, temporal, input_capacity = capacities
    short, long = list(C), list(C)
    for t, mask in enumerate(C):
        while mask:
            bit = mask & -mask
            quartet = lookup[t][bit.bit_length()-1]
            if t >= 1:
                short[t-1] |= quartet
                long[t-1] |= quartet
            if t >= 2:
                long[t-2] |= quartet
            mask -= bit
    return (input_capacity*sum(mask.bit_count() for mask in C)
            -(source-temporal)*sum(mask.bit_count() for mask in short)
            -temporal*sum(mask.bit_count() for mask in long))


def _pt_decision_check(J):
    capacities, y = (2, 1, 4), Fraction(3, 2)
    source, temporal, input_capacity = capacities
    updates, full = len(J), 15
    quartets, lookup = _pt_dependency('_ac_geometry')((((0, 1, 2, 3),),)*updates, 4)
    candidates = [t for t in range(1, updates) if J[t]]
    maximum = max(_pt_copy_score(C, lookup, capacities) for C in _pt_dependency('_ac_input_subsets')(J))
    probability_total = moment = weighted = Fraction(0)
    canonical_maximum = -4*source*updates
    offset = source*4*updates-input_capacity*sum(mask.bit_count() for mask in J)
    for bits in range(1 << len(candidates)):
        U, C = [0]*updates, [J[0]]+[0]*(updates-1)
        for i, t in enumerate(candidates):
            if bits >> i & 1:
                U[t], C[t] = full, J[t]
        V = [full ^ mask for mask in U]+[full, full]
        score = _pt_copy_score(C, lookup, capacities)
        canonical_maximum = max(canonical_maximum, score)
        probability = weight = Fraction(1)
        cost = absent_gain = 0
        for t in range(updates):
            A, W = V[t+1], V[t+1] & V[t+2]
            cost += (source-temporal)*((full ^ J[t]) & A).bit_count()+temporal*((full ^ J[t]) & W).bit_count()
            if t == 0:
                continue
            g, p, absent, loss, gain = general_flow_weights(y, capacities, A.bit_count(), W.bit_count(),
                                                           J[t].bit_count(), (J[t] & A).bit_count(), (J[t] & W).bit_count())
            assert g+p == 1 and absent+y**(-loss)*p == 1
            probability *= p if U[t] else g
            weight *= p if U[t] else absent
            if not U[t]:
                absent_gain += gain
        assert score+absent_gain == cost-offset
        probability_total += probability
        moment += probability*y**cost
        weighted += weight*y**score
    assert canonical_maximum == maximum and probability_total == 1
    assert weighted*y**offset == moment and y**(maximum+offset) <= moment
    return 1 << len(candidates)


def partial_temporal_checks():
    counts = dict(local_coefficients=0, physical_actions=0, baseline_equivalence=0,
                  unit_identities=0, decision_instances=0, decision_branches=0, callback_checks=0)
    for capacities, y in (((1, 1, 2), Fraction(2)), ((2, 1, 4), Fraction(3, 2))):
        for width in (4, 8):
            for x, alpha, beta in ((Fraction(1), Fraction(1), Fraction(1)),
                                   (Fraction(9, 4), Fraction(15, 13), Fraction(6, 5)),
                                   (Fraction(81, 16), Fraction(3, 2), Fraction(2, 3))):
                quartic = {}
                for w, e in PT_TERMS:
                    u = w+e
                    local = _pt_exact_local(y, capacities, u, w, alpha, beta)
                    # Explicit local labels give the same polynomial degree by degree.
                    A, W = (1 << u)-1, (1 << w)-1
                    physical = [Fraction(0)]*5
                    for J in range(16):
                        g, p, _, _, _ = general_flow_weights(y, capacities, u, w, J.bit_count(),
                                                            (J & A).bit_count(), (J & W).bit_count())
                        cost = (capacities[0]-capacities[1])*((15 ^ J) & A).bit_count()+capacities[1]*((15 ^ J) & W).bit_count()
                        physical[J.bit_count()] += y**cost*(p+g*alpha**4*beta**u)
                    assert physical == local
                    counts['local_coefficients'] += 5
                    if capacities == (1, 1, 2):
                        old = _pt_dependency('adaptive_cover_local_coefficients')(y, alpha, beta, u, w)
                        assert local == [v*alpha**u*beta**w for v in old]
                    quartic[(w, e)] = comb(4, u)*comb(u, w)*sum(v*x**j for j, v in enumerate(local))
                power = _pt_dependency('_ap_exact_power')(quartic, width//4)
                found = general_flow_population_action(float(y), float(x), float(alpha), float(beta), width, capacities=capacities)
                old_action = (_pt_dependency('adaptive_population_action')(float(y), float(x), float(alpha), float(beta), width)
                              if capacities == (1, 1, 2) else None)
                for i, (aa, ww) in enumerate(zip(found['a'], found['w'])):
                    ell, omega = int(aa), int(ww)
                    exact = power[(omega, ell-omega)]/(comb(width, ell)*comb(ell, omega))
                    assert exact == _pt_physical_action(width, y, x, alpha, beta, capacities, ell, omega)
                    log_expected = np.log2(PT_LD(exact.numerator)/exact.denominator)
                    assert abs(found['log2_action'][i]-log_expected) < PT_LD('1e-10')
                    counts['physical_actions'] += 1
                    if old_action is not None:
                        assert abs(old_action['log2_action'][i]-found['log2_action'][i]) < PT_LD('1e-10')
                        counts['baseline_equivalence'] += 1
                    if alpha == beta == 1:
                        endpoint = ((1+x)**(width-ell)*(y**(capacities[0]-capacities[1])+x)**(ell-omega)
                                    *(y**capacities[0]+x)**omega)
                        assert exact == endpoint
                        counts['unit_identities'] += 1
    for updates in (2, 3):
        for positions in combinations(range(4*updates), 2*updates):
            J = [0]*updates
            for position in positions:
                J[position//4] |= 1 << (position % 4)
            counts['decision_branches'] += _pt_decision_check(tuple(J))
            counts['decision_instances'] += 1
    for J in ((0, 0), (15, 15), (1, 0, 3)):
        counts['decision_branches'] += _pt_decision_check(J)
        counts['decision_instances'] += 1
    # Exercise both endpoint consumers with the explicit research callback.
    seen = []
    def checked_endpoint(prepared, y, x):
        seen.append(x)
        return partial_temporal_endpoint(prepared, y, x)
    prepared = _pt_dependency('adaptive_envelope_prepare')(4, PT_LD(15)/13, PT_LD(6)/5)
    cache = _pt_dependency('adaptive_envelope_cache')(
        prepared, PT_LD(3)/2, (1, PT_LD(9)/4), 2, 4, False,
        action_callback=partial_temporal_action, endpoint_callback=checked_endpoint)
    _pt_dependency('adaptive_envelope_endpoint')(cache)
    assert len(seen) == 2
    _pt_dependency('_ae_single_mode_controls')(cache, 2, 2)
    assert len(seen) == 4
    counts['callback_checks'] += 4
    # A mismatched supplied table must be rejected before the power engine.
    bad = general_flow_quartic(PT_LD(3)/2, 1, 1, 1)
    try:
        _pt_dependency('adaptive_population_action')(2, 1, 1, 1, 4, quartic=bad)
    except ValueError:
        counts['callback_checks'] += 1
    else:
        raise AssertionError('Mismatched quartic metadata was accepted')
    assert (Fraction(3, 2)+Fraction(9, 4))/(1+Fraction(9, 4)) == Fraction(15, 13)
    assert (Fraction(9, 4)+Fraction(9, 4))/(Fraction(3, 2)+Fraction(9, 4)) == Fraction(6, 5)
    print(dict(**counts, exact_width8_partitions=35, encoder_changed=False,
               shared_positive_power_engine=True, runtime_globals_monkeypatched=False,
               full_code_certificate=False), flush=True)


def partial_temporal_main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--checks', action='store_true')
    mode.add_argument('--control', action='store_true')
    parser.add_argument('--iterations', type=int, choices=(1, 2), default=2)
    args = parser.parse_args()
    import resource
    if shutil.disk_usage('/tmp').free < 20*(1 << 30):
        raise RuntimeError('20GiB /tmp free-space guard')
    resource.setrlimit(resource.RLIMIT_AS, (2*(1 << 30), 2*(1 << 30)))
    timer = threading.Timer(175, lambda: os._exit(124))
    timer.daemon = True
    timer.start()
    try:
        if args.checks:
            partial_temporal_checks()
        else:
            print(partial_temporal_control(args.iterations, True), flush=True)
    finally:
        timer.cancel()


if __name__ == '__main__':
    partial_temporal_main()
