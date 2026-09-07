"""Exact ledger and outward kernel check for the t256_s14 zero-state obstruction.

Reuses the even-row convolution argument in ZERO_STATE_FIRST_MOMENT_OBSTRUCTION.md.
The conclusion concerns E[bad messages], not Pr[any bad message].
No frozen producer or old certificate is changed.
"""
import argparse
import math
from pathlib import Path

from flint import arb, arb_poly, ctx

import certificate_search_core as core
import exact_tail_lower_bound as tail_check

base = core.base
OUTPUT = base.HERE / 'generated/t256_s14_zero_state_obstruction_v1.json'
TAIL = base.HERE / 'generated/joint_tail_lower_80/lower.json'


def calculation(precision):
    spec = core.instance('t256_s14', 20)
    ctx.prec = precision
    t, s, _, kernel = core.inputs.load(spec['configuration'])
    rows, q, lo, hi = spec['rows'], 2620, 64, 1536
    core.require((t, s, rows) == (256, 14, 8192), 'Wrong instance')
    core.require(q % 2 == 0 and 80*q <= spec['cutoff'], 'Invalid row family')
    core.require(sum(kernel.values()) == 2**(t-s), 'Kernel dimension mismatch')
    core.require(all(j % 2 == 0 for j, n in kernel.items() if n), 'Odd kernel word')
    polynomial = arb_poly([kernel.get(j, 0) for j in range(t+1)])**(rows//t)
    minimum = min((polynomial[j]/math.comb(rows, j)).lower()
                  for j in range(lo, hi+1, 2))
    region_exponent = s*(rows//t)-1
    core.require(minimum > arb(2)**(-region_exponent), 'Kernel lower bound failed')

    def tail(a, p):
        divergence = a*(a/p).log() + (1-a)*((1-a)/(1-p)).log()
        return (-q*divergence).exp()

    bad = 256*(tail(arb(lo-1)/q, arb(38)/256) +
               tail(arb(hi+1)/q, arb(80)/256))
    gamma = arb(2)**(-255)-bad
    core.require(gamma > arb(2)**(-256), 'Parity-and-range lower bound failed')
    tail_receipt = base.read(TAIL)
    core.require(tail_receipt['status'] == 'EXACT_BCH_TAIL_LOWER' and
                 tail_receipt['weights'] == list(range(38, 81, 2)) and
                 tail_receipt['rational_primal_dual_checks_passed'], 'Wrong BCH tail')
    for root, field in ((base.HERE, 'local_sha256'), (base.BCH, 'outer_sha256')):
        for name, digest in tail_receipt[field].items():
            core.require(base.sha(root/name) == digest, 'Changed tail dependency: '+name)
    a = tail_receipt['lower']
    denominator_exponent = 256*(region_exponent+1)
    numerator = math.comb(rows, q)*a**q
    claim_bits = 150700
    core.require(numerator > 1 << (denominator_exponent+claim_bits),
                 'Exact first-moment lower bound failed')
    # Store the small exact formula, not its roughly 80 KB expanded integer.
    return dict(instance=spec, occupation=q, row_weight_interval=[38, 80],
                region_weight_interval=[lo, hi], input_weight_upper=80*q,
                tail_cardinality_lower=a, region_probability_strict_lower_power=-region_exponent,
                parity_and_range_probability_strict_lower_power=-256,
                first_moment_strict_lower_formula='binom(rows, occupation) * tail_cardinality_lower^occupation / 2^denominator_exponent',
                denominator_exponent=denominator_exponent,
                first_moment_exceeds_power_of_two=claim_bits,
                setup_failure_probability_lower_bound=False)


def sources():
    hashes = core.inputs.source_hashes()
    hashes.update(core.outer_dependencies())
    for p in (Path(__file__), Path(core.__file__), Path(tail_check.__file__), TAIL,
              base.HERE/'ZERO_STATE_FIRST_MOMENT_OBSTRUCTION.md'):
        hashes[p.relative_to(base.ROOT).as_posix()] = base.sha(p)
    tail = base.read(TAIL)
    for root, field in ((base.HERE, 'local_sha256'), (base.BCH, 'outer_sha256')):
        for name, digest in tail[field].items():
            hashes[(root/name).relative_to(base.ROOT).as_posix()] = digest
    return hashes


def run(verify=False):
    if verify:
        saved = base.read(OUTPUT)
        core.require(saved['status'] == 'ACTUAL_FIRST_MOMENT_OBSTRUCTION_NOT_FAILURE_PROBABILITY',
                     'Wrong obstruction status')
        core.authenticate(saved, base.ROOT)
    else:
        core.require(not OUTPUT.exists(), 'Use --verify for an existing certificate')
    # Recheck the exact LP solution, not just its success flag or file hash.
    tail_check.run(verify=True)
    result = calculation(512 if verify else 256)
    if verify:
        core.require(result == saved['result'], 'Obstruction replay mismatch')
    else:
        base.write_new(OUTPUT, dict(
            status='ACTUAL_FIRST_MOMENT_OBSTRUCTION_NOT_FAILURE_PROBABILITY',
            result=result, source_sha256=sources()))
    print(('512-bit replay passed. ' if verify else '256-bit producer passed. ')+
          'E[bad messages] > 2^150700 at K=2^20, Q=2620.', flush=True)
    print('This rules out unconditional first-moment closure, not low setup failure probability.', flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--verify', action='store_true')
    run(p.parse_args().verify)
