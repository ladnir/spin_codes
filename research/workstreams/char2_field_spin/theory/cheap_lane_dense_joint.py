"""Directed 35% joint bound for actual support sizes7552..8192.

Same encoder, ideal setup and input-restricted kernel as the frozen
cheap_lane_dense_gpadding.py. Only the initial source comparison changes:
average missing outer maps together while retaining nonzero projections.
With h=32 and r<=8192, c_h>=1/33, so the new factor is
 beta^r *33/255 * (1+sum_{j=1}^32 C(8192,j)/255).
No q^256 or g-padding cost is charged for the missing groups. The joint
source lemma is applied before the existing local kernel majorization.
For every route the allowed-input total is256r; hence K gives
 A^r*q^8*rho^32767, A=beta*g^256<1. Summing at most2^8192 supports with
 r>=7552 gives the expression below. Terminal is output-only.

This is not a whole-code certificate. Bounded10s,20GiB free,stdout only.
"""

import argparse
import hashlib
import json
from math import comb
from pathlib import Path
import shutil
import sys
import time

sys.dont_write_bytecode = True


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--precision', type=int, choices=(256, 384), default=256)
    args = parser.parse_args()
    source = Path(__file__).resolve()
    if shutil.disk_usage(source.anchor).free < 20*(1 << 30):
        raise RuntimeError('requires20GiB free')
    started = time.monotonic()
    from flint import arb, ctx
    ctx.prec = args.precision
    q = arb(256)
    groups, first, h, packets, records = 8192, 7552, 32, 32768, 1 << 21
    limit = 7*records//20
    assert h <= 256 and h < first and groups < 256*(h+1)
    assert 256*(h+1)-groups == 256
    beta = q**160/arb(255)**32
    z = arb(3)/1250
    g = z+(1-z)/q
    group_factor = beta*g**256
    assert 0 < group_factor < 1
    missing = 1+sum(arb(comb(groups, j)) for j in range(1, h+1))/255
    rho = (1+1/q)**16+65536*q**(-8)
    bound = (arb(2)**groups*33/255*group_factor**first*missing
             *q**8*rho**(packets-1)*z**(-limit))
    bits = bound.log()/arb(2).log()
    assert bound.is_finite() and bound > 0
    if time.monotonic()-started > 10:
        raise TimeoutError('10-second dense-joint budget exceeded')
    print(json.dumps({
        'scope': 'all actual source supports of size7552..8192, all state paths',
        'whole_code_certificate': False, 'verifier_completed': True,
        'directed': True, 'precision_bits': args.precision,
        'parameters': {'distance_limit': limit, 'z': '3/1250',
                       'source_min': first, 'source_max': groups,
                       'maximum_missing_groups': h, 'uniform_c_lower': '1/33'},
        'source_factor': 'beta^r*33/255*(1+sum_{j=1}^32 C(8192,j)/255)',
        'log2_bound_interval': str(bits),
        'proved_below_2_pow_minus_40': bool(bound < arb(2)**-40),
        'proved_below_2_pow_minus_768': bool(bound < arb(2)**-768),
        'tests': {'positive_group_factor_less_than_one': True,
                  'joint_missing_h_at_most_q': True,
                  'uniform_c_lower_exact': True, 'charged_refreshes': packets-1},
        'generated_bytes': 0, 'seconds': time.monotonic()-started,
        'source_sha256': hashlib.sha256(source.read_bytes()).hexdigest()
    }, indent=2))


if __name__ == '__main__':
    main()
