"""Run the existing sparse verifier with a scaled floating proposal helper.

This process-local adapter only changes witness search. The support partition,
fresh counts, directed arithmetic, and final acceptance test remain unchanged.
The analytic-gradient helper already rescales before its first squaring; using
zero derivative directions supplies the same stable scalar moment here.
"""
from math import log
import runpy
import sys

import numpy as np
from flint import arb, arb_mat, ctx
import sparse_cover
import occupancy_cdf_cover as cover
from cdf_gradient import log_power_gradient


def stable_moment(matrix, length=128, terminal=None):
    matrix = np.asarray(matrix)
    if terminal is None:
        terminal = np.ones(len(matrix))
    derivatives = np.empty((0, *matrix.shape))
    return log_power_gradient(matrix, derivatives, length, terminal)[0]


def self_test():
    for magnitude in (1e-20, 1e-200, 1e-300):
        for length in (1, 3, 256):
            result = stable_moment(np.array([[magnitude]]), length)
            expected = length*log(magnitude)
            assert abs(result-expected) < 1e-8
    previous = ctx.prec
    ctx.prec = 256
    try:
        for scale in (1., 1e-200):
            matrix = np.array([[.2,.1],[.3,.4]])*scale
            for length in (3,256):
                for terminal in (np.ones(2),np.array([0.,1.])):
                    exact = arb_mat([[arb(float(x)) for x in row] for row in matrix])**length
                    expected = sum((exact[0,j]*int(terminal[j]) for j in range(2)),arb(0)).log()
                    assert abs(stable_moment(matrix,length,terminal)-float(expected)) < 1e-8
        for matrix in (np.zeros((2,2)),np.array([[float('nan')]])):
            try:
                stable_moment(matrix)
            except ArithmeticError:
                pass
            else:
                raise AssertionError('invalid zero/nonfinite proposal was accepted')
    finally:
        ctx.prec = previous
    print('Stable sparse proposal helper: 17 positive cases and 2 rejection cases passed.',flush=True)


def main():
    if sys.argv[1:] == ['--self-test']:
        self_test()
        return
    original = cover.log_power_moment
    cover.log_power_moment = stable_moment
    try:
        print('PROPOSAL ADAPTER: initial matrix rescaling; original outward verifier unchanged.',flush=True)
        runpy.run_module('shared_sparse',run_name='__main__')
    finally:
        cover.log_power_moment = original


if __name__=='__main__':
    main()
