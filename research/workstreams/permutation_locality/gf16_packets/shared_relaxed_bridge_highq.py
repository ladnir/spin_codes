"""Selected high-occupancy bridge covers with log-scaled proposals.

Each exact regional coefficient is normalized before binary64 conversion.
Bernoulli masses are combined in log space. These changes affect only witness
search; the existing support-cover verifier still replays exact coefficients.
"""
import argparse
from collections import Counter
from contextlib import contextmanager
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
from scipy.special import gammaln
from flint import arb, arb_mat, ctx
import sparse_cover
import occupancy_cdf_cover as cover
import occupancy_birth_classes
from shared_support import shared_counts
from shared_relaxed_bridge_sparse import stable_moment


class ScaledRegion:
    def __init__(self, exact):
        self.logs = []
        self.values = []
        for matrix in exact:
            size = matrix.nrows()
            peak = max(matrix[i,j].upper() for i in range(size) for j in range(size))
            if not peak > 0:
                raise ArithmeticError('positive exact regional coefficient required')
            self.logs.append(float(peak.log()))
            self.values.append(np.array([[float(matrix[i,j]/peak)
                                          for j in range(size)] for i in range(size)]))
        self.logs = np.asarray(self.logs)
        self.values = np.asarray(self.values)

    def __len__(self):
        return len(self.logs)


def log_masses(probabilities):
    """Log Poisson-binomial masses, grouping identical probabilities."""
    masses = np.array([0.])
    for probability, count in Counter(probabilities).items():
        if not 0 < probability < 1:
            raise ValueError('strictly interior Bernoulli probabilities required')
        support = np.arange(count+1)
        factor = (gammaln(count+1)-gammaln(support+1)-gammaln(count-support+1)
                  +support*np.log(probability)+(count-support)*np.log1p(-probability))
        if len(masses) == 1:
            masses = factor
        else:
            result = np.full(len(masses)+count, -np.inf)
            for index, value in enumerate(factor):
                np.logaddexp(result[index:index+len(masses)], value+masses,
                             out=result[index:index+len(masses)])
            masses = result
    return masses


def scaled_matrix(region, probabilities):
    if not isinstance(region, ScaledRegion):
        raise TypeError('the normalized bridge adapter requires ScaledRegion')
    count = len(probabilities)+1
    weights = log_masses(probabilities)+region.logs[:count]
    peak = float(np.max(weights))
    matrix = np.einsum('i,ijk->jk', np.exp(weights-peak), region.values[:count])
    return matrix, peak


def scaled_moment(value, length=128, terminal=None):
    matrix, scale = value
    return stable_moment(matrix, length, terminal)+length*scale


@contextmanager
def proposal_adapter():
    original_matrix = cover.matrix_for_probabilities
    original_moment = cover.log_power_moment
    cover.matrix_for_probabilities = scaled_matrix
    cover.log_power_moment = scaled_moment
    try:
        yield
    finally:
        cover.matrix_for_probabilities = original_matrix
        cover.log_power_moment = original_moment


def self_test():
    ctx.prec = 256
    cases = 0
    for scales in ((0,-400,-800,-1200), (-1200,-900,-500,-100)):
        exact = [arb_mat([[arb('0.6'),arb('0.1')],[arb('0.3'),arb('0.4')]])
                 *arb(scale).exp() for scale in scales]
        region = ScaledRegion(exact)
        for probabilities in ([.1,.1,.1], [.01,.6,.6], [.99,.99,.01]):
            masses = [arb(1)]
            for probability in probabilities:
                p = arb(probability)
                nxt = [arb(0)]*(len(masses)+1)
                for i,value in enumerate(masses):
                    nxt[i] += value*(1-p)
                    nxt[i+1] += value*p
                masses = nxt
            matrix = sum((m*r for m,r in zip(masses,exact)),arb_mat(2,2))
            for length in (1, 256):
                power = matrix**length
                expected = float((power[0,0]+power[0,1]).log())
                actual = scaled_moment(scaled_matrix(region,probabilities),length,np.ones(2))
                assert abs(actual-expected) < 1e-7, (actual,expected)
                cases += 1
    print('Normalized proposal adapter:',cases,'Arb comparison cases passed.',flush=True)


def run(options):
    ctx.prec = options.precision
    witnesses = [Path(p) for p in options.count_witnesses]
    counts, _ = shared_counts(True, True, False, False, witnesses)
    args = sparse_cover.build_args(options.occupancy,options.tilts,options.precision,
                                   options.max_splits,32,None,2)
    args.exact_feedback = True
    args.analytic_gradient = False
    args.joint_return_through = 3
    args.lazy_density_through = 6
    args.proposal_buffer_bits = .125
    operators = occupancy_birth_classes.build_operators(args)
    operators = {choice:(exact,ScaledRegion(exact))
                 for choice,(exact,_) in operators.items()}
    terminal = np.ones(next(iter(operators.values()))[0][0].nrows())
    receipt = dict(schema='shared-gf16-highq-bridge-1',occupancy=options.occupancy,
                   updates=2,precision=options.precision,tilts=options.tilts,
                   requested_cutoffs=options.cutoffs,max_splits=options.max_splits,
                   target_bits=32,refined_counts=True,coupled_counts=True,
                   joint_return_through=3,lazy_density_through=6,
                   count_witnesses=[dict(path=str(p),sha256=sha256(p.read_bytes()).hexdigest())
                                    for p in witnesses],results=[],
                   note='All supports at each successful listed occupancy/cutoff only. '
                        'Normalized floating proposals; original exact outward cover verifier. '
                        'No whole-code claim.')
    with proposal_adapter():
        for cutoff in options.cutoffs:
            print('NORMALIZED BRIDGE COVER occupancy',options.occupancy,'cutoff',cutoff,flush=True)
            upper = sparse_cover.sparse.cover(args,operators,{'1':counts},terminal,cutoff=cutoff)
            receipt['results'].append(dict(cutoff=cutoff,
                upper=None if upper is None else [int(v) for v in upper.upper().man_exp()]))
            options.output.parent.mkdir(parents=True,exist_ok=True)
            options.output.write_text(json.dumps(receipt,indent=2)+'\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--self-test',action='store_true')
    parser.add_argument('--occupancy',type=int,default=256)
    parser.add_argument('--tilts',nargs='+',default=['.16','.192','.224','.256'])
    parser.add_argument('--cutoffs',type=int,nargs='+',default=[136314,146800])
    parser.add_argument('--precision',type=int,default=256)
    parser.add_argument('--max-splits',type=int,default=4)
    parser.add_argument('--count-witnesses',nargs='+')
    parser.add_argument('--output',type=Path)
    options = parser.parse_args()
    if options.self_test:
        self_test()
    else:
        if not options.count_witnesses or options.output is None:
            parser.error('count witnesses and output path are required')
        if not 1 <= options.occupancy <= 2048 or options.precision < 128 or options.max_splits < 0:
            parser.error('valid occupancy, precision, and split limit required')
        run(options)
