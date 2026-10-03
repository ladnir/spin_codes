"""Optimistic return ablations: diagnostics only, never probability bounds."""
import argparse
import json
from math import comb, log

import numpy as np

import refresh_gate as gate


def run(q=119, tilt=.4):
    data, _ = gate.maps.prepare()
    z = np.exp(-tilt)
    weighted, emission, _ = gate.prior.moments(data, z)
    local, record = gate.operators(data, weighted, emission, z, mode='shared',
        census=gate.exact_returns(data), holder=gate.prepare_holder(data))
    beta = gate.prior.rs_uniform_envelope.UniformInputEnvelope(16, 8, 8, 2).beta
    beta_log = log(beta.numerator)-log(beta.denominator)
    outcomes = {}
    for name, keep_through in [('all_bounds', 8), ('exact_j1_j2_only', 2), ('exact_j1_only', 1), ('no_returns', 0)]:
        candidate = local.copy()
        if keep_through < 8:
            candidate[keep_through+1:, 1, 0] = 0
            candidate[keep_through+1:, 1, 1] = np.array(record['uniform_emission'])[keep_through+1:]
        regional, _ = gate.prior.placement(candidate, q, epochs=64, windows=8)
        moment = gate.prior.logarithmic.log_power_matrix(gate.prior.uniform_mixture(regional, q), 32)
        outcomes[name] = -(log(comb(512, q))+q*beta_log+tilt*13107+moment)/log(2)
    total = np.array(record['uniform_emission'])
    ratios = np.array(record['return_upper'])*65535/total
    return dict(q=q, tilt=tilt, margins=outcomes,
        shared_return_over_independent_return=ratios.tolist(),
        bounds_valid_only_for_all_bounds=True, optimistic_ablations_are_not_certificates=True,
        local_diagnostic=record)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--q', type=int, default=119)
    parser.add_argument('--tilt', type=float, default=.4)
    args = parser.parse_args()
    print(json.dumps(run(args.q, args.tilt), indent=2))
