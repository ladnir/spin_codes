"""Diagnostic ablations of the row-conditioned transfer, never certificates.

Deleting or scaling a transition need not preserve an upper bound. These
experiments locate candidate sources of slack; they cannot certify distance
or exhibit a bad codeword. The production verifier does not import this file.
"""
import argparse
from math import comb, log
from fractions import Fraction as Q

import numpy as np
from row_verify import packet_law
from row_counts import row_gamma_function
from shape_inner import build, prepare_shared
from averaged_windows import AveragedHighInner
from cone_moment import log_moment
from bch_joint_support import authenticated_caps
from occupancy_sensitivity import float_placement
from occupancy_screen import matrix_for_probabilities
from occupancy_memory import Z, F, M, C, U
from mature_tail import L48, L56


def ablations():
    yield 'baseline', 0, [], 1.
    for name, pairs in (
        ('density cancellation', [(C,Z)]),
        ('fresh cancellation', [(F,Z)]),
        ('uniform cancellation', [(U+i,Z) for i in range(5)]),
        ('mature refresh cancellation', [(i,Z) for i in (M,L48,L56)]),
        ('direct zero return', [(Z,Z)]),
        ('density persistence', [(C,C)]),
        ('uniform density', [(U+i,C) for i in range(5)]),
        ('fresh density', [(F,C)]),
        ('mature lazy mass', [(i,M) for i in (M,L48,L56)]),
        ('tail arrival', [(i,j) for i in range(11) for j in (L48,L56)]),
        ('all nonempty returns', [(i,Z) for i in range(11)]),
    ):
        for factor in (.5,0.):
            yield name, 1, pairs, factor
    # Separate the exact joint-cancellation range from higher occupancies.
    for first in (4,9):
        yield 'higher occupancy zero returns', first, [(i,Z) for i in range(11)], 0.


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--occupancies',type=int,nargs='+',default=[64,80,128])
    p.add_argument('--tilt',default='.052')
    p.add_argument('--p',type=Q,default=Q(1,2))
    p.add_argument('--cut',type=int,default=8)
    p.add_argument('--full-feedback',type=int,default=6)
    p.add_argument('--window-histogram',type=int,default=8)
    p.add_argument('--precision',type=int,default=192)
    p.add_argument('--column-density',action='store_true')
    p.add_argument('--feedback-density',type=int,choices=range(0,7),default=0)
    p.add_argument('--snapshot',help='save diagnostic floats only, never consumed by a verifier')
    p.add_argument('--load-snapshot',help='reuse matching diagnostic floats; no proof claim')
    p.set_defaults(class_tail=True,joint_cancellation=True,penalty='1')
    args = p.parse_args()
    if not 0 < args.p <= 1 or any(not 1 <= q <= 2048 for q in args.occupancies):
        p.error('invalid probability or group count')
    caps = authenticated_caps()
    active,theta = packet_law(float(args.p),4)
    parameters = tuple(map(str,(args.tilt,args.p,args.cut,args.full_feedback,args.window_histogram,
                               args.precision,args.column_density,args.feedback_density)))
    if args.load_snapshot:
        with np.load(args.load_snapshot,allow_pickle=False) as saved:
            recorded = tuple(saved['parameters'])
            # Legacy snapshots predate both optional density refinements.
            if len(recorded) == 6:
                recorded += ('False','0')
            if recorded != parameters:
                raise ValueError('diagnostic parameters differ from snapshot')
            base = saved['base']
        if base.shape != (33,11,11) or not np.isfinite(base).all() or (base<0).any():
            raise ValueError('invalid diagnostic snapshot')
    else:
        model = AveragedHighInner(build(args,prepare_shared(args)))
        base = model.float_mix(theta)
    if args.snapshot:
        np.savez_compressed(args.snapshot,base=base,parameters=np.asarray(parameters))
    gamma = row_gamma_function(caps,64,192,exclude_zero=True)(float(args.p))
    maximum = max(args.occupancies)
    for name,first,pairs,factor in ablations():
        changed = base.copy()
        for i,j in pairs:
            changed[first:,i,j] *= factor
        regions = float_placement(changed,maximum)
        scores = {}
        for q in args.occupancies:
            matrix = matrix_for_probabilities(regions[:q+1],[active]*q)
            scores[q] = (log_moment(matrix)+float(args.tilt)*209715
                         +log(comb(2048,q))+4*q*gamma)/log(2)
        print('DIAGNOSTIC, NOT A BOUND:',name,'local j >=',first,'factor',factor,
              'log2 scores',scores,flush=True)


if __name__ == '__main__':
    main()
