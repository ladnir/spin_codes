"""Numerical convergence of empty/one-impulse IMT regions; not a limit proof."""
import argparse
import json
import math
from pathlib import Path

import numpy as np
import screen


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    args = p.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    model = screen.Model()
    e = model.engine
    n = e.n
    left = np.zeros((2,n));left[0,0]=1;left[1,2:]=model.counts/model.m
    right = np.zeros((n,2));right[0,0]=1;right[1:,1]=1
    rows = []
    for theta in (2.,8.):
        gamma = theta*(1<<18)/model.m
        gap = -math.expm1(-gamma)/gamma
        limit0 = np.diag([1.,math.exp(-gamma)])
        limit1 = np.array([[0.,gap],[gap/model.m,math.exp(-gamma)]])
        for exponent in (12,16,20,24,28):
            length = 1<<exponent
            epochs = np.exp(screen.candidate.model.independent.g.epochs(e.spectrum,e.kernel,e.columns,e.caps,e.low,theta/length,maximum=1))
            transfer = np.block([[epochs[0],epochs[1]],[np.zeros((n,n)),epochs[0]]])
            raised = np.linalg.matrix_power(transfer,length//128)
            actual0 = left @ raised[:n,:n] @ right
            actual1 = left @ (raised[:n,n:]/(length//128)) @ right
            error0 = float(np.max(np.abs(actual0-limit0)))
            relative1 = float(np.max(np.abs(actual1-limit1)/np.maximum(limit1,1/model.m)))
            rows.append(dict(theta=theta,region_length=length,empty_max_absolute_error=error0,
                             one_impulse_max_relative_error=relative1))
    result = dict(status='BINARY64_IMT_CONTINUUM_DIAGNOSTIC_NOT_PROOF',rows=rows,
                  limit_impulse_upper=[[0.,1.],[1/model.m,1.]],
                  source_sha256={p.relative_to(screen.ROOT).as_posix():screen.sha(p) for p in [Path(__file__),Path(screen.__file__)]},
                  limitations=['Only empty and one-impulse regions, with zero or stationary live initial distributions.',
                    'No uniform error bound over spacings or arbitrary entering states.'])
    with args.output.open('x',encoding='utf-8') as stream:
        json.dump(result,stream,indent=2);stream.write('\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
