"""Bounded diagnostic for separable causal-flow exponential potentials.

Five-state spectral control only; floating point, not a certificate.
No code change, output files, or large state enumeration. Standard library.
"""
import math
import shutil
import time
import sys
from fractions import Fraction as F
from pathlib import Path
from random_pairing_demand_flow import causal_pair

P = 1974 / 8192
LAW = ((1-P)**2, 2*P*(1-P), P*P)
STEPS = 4095
PAIRS = 256
LOSS_THRESHOLD = 371901
WITNESS_ROWS = (
    ((611302,336623,52075),(595844,346845,57311),(578284,357838,63878),
     (560519,370472,69009),(541404,386585,72011)),
    ((605409,340581,54010),(589289,351240,59471),(571406,362460,66134),
     (553718,375015,71267),(535180,390593,74227)),
)


def perron(theta, eta):
    entries = []
    for a in range(-2, 3):
        row = []
        for j, probability in enumerate(LAW):
            na, nb, _, loss = causal_pair(a, a, j)
            assert na == nb
            row.append((na+2, probability*math.exp(theta*loss+eta*(j-2*P))))
        entries.append(row)
    vector = [1.]*5
    for iteration in range(2000):
        result = [sum(value*vector[b] for b, value in row) for row in entries]
        ratios = [x/y for x, y in zip(result, vector)]
        low, high = min(ratios), max(ratios)
        scale = max(result)
        vector = [x/scale for x in result]
        if high-low < 2e-13*high:
            return (high+low)/2, vector, high-low
    raise AssertionError('Perron iteration did not converge')


def optimistic_objective(point):
    theta, eta = point
    if theta < 0 or theta > 2 or abs(eta) > 2:
        return float('inf')
    value, _, _ = perron(theta, eta)
    # Missing first input batch is retained when conditioning all route tags.
    initial = 512*math.log((1-P)*math.exp(-eta*P)+P*math.exp(eta*(1-P)))
    return STEPS*PAIRS*math.log(value)-theta*LOSS_THRESHOLD+initial


def search():
    point = [0.05, 0.05]
    best = optimistic_objective(point)
    evaluations = 1
    for step in (0.08, 0.04, 0.02, 0.01, 0.005, 0.0025, 0.00125,
                 0.000625, 0.0003125, 0.00015625, 0.000078125):
        for _ in range(80):
            choices = []
            for dx, dy in ((1,0),(-1,0),(0,1),(0,-1),(1,1),(1,-1),(-1,1),(-1,-1)):
                trial = [point[0]+step*dx, point[1]+step*dy]
                value = optimistic_objective(trial); evaluations += 1
                choices.append((value, trial))
            value, trial = min(choices)
            if value >= best:
                break
            best, point = value, trial
    return point, best, evaluations


def stationary(rows):
    matrix = [[F(0) for _ in range(5)] for _ in range(5)]
    for a in range(5):
        for j in range(3):
            na, nb, _, _ = causal_pair(a-2,a-2,j)
            assert na == nb
            matrix[a][na+2] += rows[a][j]
    equations = [[matrix[j][i]-int(i == j) for j in range(5)]+[F(0)] for i in range(4)]
    equations.append([F(1)]*5+[F(1)])
    for j in range(5):
        pivot = next(i for i in range(j,5) if equations[i][j])
        equations[j], equations[pivot] = equations[pivot], equations[j]
        value = equations[j][j]; equations[j] = [x/value for x in equations[j]]
        for i in range(5):
            if i != j:
                value = equations[i][j]
                equations[i] = [x-value*y for x,y in zip(equations[i],equations[j])]
    mu = [equations[i][-1] for i in range(5)]
    assert sum(mu) == 1 and min(mu) > 0
    assert all(sum(mu[a]*matrix[a][b] for a in range(5)) == mu[b] for b in range(5))
    flow = [[mu[a]*rows[a][j] for j in range(3)] for a in range(5)]
    return flow


def rational_tilt(theta,eta):
    _, vector, _ = perron(theta,eta)
    rows = []
    for a in range(-2,3):
        raw = []
        for j, probability in enumerate(LAW):
            na, _, _, loss = causal_pair(a,a,j)
            raw.append(probability*math.exp(theta*loss+eta*(j-2*P))*vector[na+2])
        mass = sum(raw); denominator = 1000000
        numerators = [round(x/mass*denominator) for x in raw[:2]]
        numerators.append(denominator-sum(numerators))
        assert min(numerators) > 0
        rows.append([F(x,denominator) for x in numerators])
    return rows, stationary(rows)


def entropy_witness():
    # Two exact stationary flows, mixed to preserve the original mean tag count.
    pair_mean = F(1974,4096)
    cases = [rational_tilt(.082, eta) for eta in (.16,.18)]
    means = [sum(j*x for row in flow for j,x in enumerate(row)) for _,flow in cases]
    weight = (pair_mean-means[1])/(means[0]-means[1])
    assert 0 <= weight <= 1, [float(x) for x in means]
    flow = [[weight*cases[0][1][a][j]+(1-weight)*cases[1][1][a][j]
             for j in range(3)] for a in range(5)]
    mu = [sum(row) for row in flow]
    rows = [[x/mu[a] for x in row] for a,row in enumerate(flow)]
    assert stationary(rows) == flow
    mean_tags = sum(j*x for row in flow for j,x in enumerate(row))
    loss = 2*flow[0][0]
    assert mean_tags == pair_mean
    original = [(1-F(1974,8192))**2,2*F(1974,8192)*(1-F(1974,8192)),F(1974,8192)**2]
    chi = sum(mu[a]*(rows[a][j]-original[j])**2/original[j] for a in range(5) for j in range(3))
    entropy_float = sum(float(flow[a][j])*math.log(float(rows[a][j]/original[j]))
                        for a in range(5) for j in range(3))
    print(dict(witness_rows=[[[int(x*1000000) for x in row] for row in case[0]] for case in cases],
               mixture_weight=float(weight), witness_pair_loss=float(loss),
               required_pair_loss=LOSS_THRESHOLD/(STEPS*PAIRS), exact_mean_tags=True,
               chi_square_cost_bits=float(chi*STEPS*PAIRS)/math.log(2),
               entropy_cost_bits_diagnostic=entropy_float*STEPS*PAIRS/math.log(2)))
    return loss, chi


def verify_barrier():
    """Exact rational witness; no dependence on floating-point search."""
    flows = []
    for integers in WITNESS_ROWS:
        assert all(sum(row) == 1000000 and min(row) > 0 for row in integers)
        rows = [[F(x,1000000) for x in row] for row in integers]
        flows.append(stationary(rows))
    means = [sum(j*x for row in flow for j,x in enumerate(row)) for flow in flows]
    target = F(1974,4096)
    weight = (target-means[1])/(means[0]-means[1])
    assert 0 < weight < 1
    flow = [[weight*flows[0][a][j]+(1-weight)*flows[1][a][j]
             for j in range(3)] for a in range(5)]
    mu = [sum(row) for row in flow]
    rows = [[x/mu[a] for x in row] for a,row in enumerate(flow)]
    assert stationary(rows) == flow
    assert sum(j*x for row in flow for j,x in enumerate(row)) == target
    assert 2*flow[0][0]*STEPS*PAIRS >= LOSS_THRESHOLD
    p = F(1974,8192)
    original = ((1-p)**2,2*p*(1-p),p*p)
    chi = sum(mu[a]*(rows[a][j]-original[j])**2/original[j]
              for a in range(5) for j in range(3))
    # D(Q||P) <= chi^2(Q,P), and ln(2)>69/100.
    # First two positive terms of 2*atanh(1/3) already give 56/81.
    assert F(56,81) > F(69,100)
    assert chi*STEPS*PAIRS < 5500*F(69,100)
    assert math.comb(8192,1974) > 1 << 6519
    print(dict(exact_stationarity=True, exact_tag_mean=True,
               loss_threshold_met=True, entropy_budget_less_than_bits=5500,
               rowset_union_more_than_bits=6519,
               separable_envelope_union_exponent_greater_than=1019,
               whole_code_counterexample=False, setup_tail_certificate=False))


def main():
    assert shutil.disk_usage(Path.cwd()).free >= 20*(1 << 30), 'disk guard'
    began = time.monotonic()
    if '--verify' in sys.argv:
        verify_barrier()
        print(dict(seconds=time.monotonic()-began))
        return
    theta_eta, exponent, evaluations = search()
    theta, eta = theta_eta
    value, vector, width = perron(theta, eta)
    weights = [math.sqrt(x) for x in vector]
    ratios = []
    for a in range(-2,3):
        for b in range(a,3):
            total = 0.
            for j, probability in enumerate(LAW):
                na, nb, _, loss = causal_pair(a,b,j)
                total += probability*math.exp(theta*loss+eta*(j-2*P))*weights[na+2]*weights[nb+2]
            ratios.append((total/(weights[a+2]*weights[b+2]), a, b))
    log_choose = math.lgamma(8193)-math.lgamma(1975)-math.lgamma(6219)
    log_region_probability = log_choose+1974*math.log(P)+6218*math.log1p(-P)
    overhead = log_choose-256*log_region_probability
    upper_ratio, a, b = max(ratios)
    if '--pair-ratios' in sys.argv:
        print(dict(normalized_pair_weights=[(a,b,ratio/value) for ratio,a,b in ratios]))
    boundary = 512*math.log(weights[2]/min(weights))
    trial_bound = exponent+STEPS*PAIRS*math.log(upper_ratio/value)+boundary+overhead
    print(dict(theta=theta, eta=eta, spectral_radius=value, perron_ratio_width=width,
               optimistic_fixed_pair_exponent_bits=exponent/math.log(2),
               rowset_union_bits=log_choose/math.log(2),
               regional_conditioning_bits=-256*log_region_probability/math.log(2),
               optimistic_with_overhead_bits=(exponent+overhead)/math.log(2),
               candidate_any_pair_bound_bits=trial_bound/math.log(2),
               worst_pair=(a,b), candidate_weights=weights, evaluations=evaluations,
               seconds=time.monotonic()-began, floating_screen_not_certificate=True))
    entropy_witness()
    assert time.monotonic()-began < 25


if __name__ == '__main__':
    main()
