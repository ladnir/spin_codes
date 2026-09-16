"""Bounded dense retuning for the new mixer, using direct syndrome smoothing.

Every result is binary64 discovery, not outward certification. Box geometry is
unchanged; optimization changes only the fixed proposal and Chernoff tilt.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import minimize,minimize_scalar
from scipy.special import gammaln

import general_occupancies as g

HERE=Path(__file__).resolve().parent


def optimize_box(box,dense,counts,spectrum,kernel,cutoff):
    witness=box['witness'];bank=dense['probability_banks'][witness['probability_bank']]
    ps=np.array(bank['probabilities']);costs=[0.]
    for band,p in zip(dense['bands'],ps[1:]):
        costs.append(math.log(counts[128]) if band==[128] else max(
            math.log(counts[w])-math.log(math.comb(128,w))-w*math.log(p)-(128-w)*math.log1p(-p) for w in band))
    corners=g.typed.vertices(box['lower'],box['upper'],32768)
    multinomial=gammaln(32769)-gammaln(corners+1).sum(axis=1)
    lattice=g.typed.lattice_log_count(box['lower'],box['upper'])
    evaluations=0
    def objective(x):
        nonlocal evaluations
        evaluations+=1
        logits=np.array([x[0],x[1],0.]);logs=logits-float(np.logaddexp.reduce(logits));proposal=np.exp(logs)
        theta=float(ps@proposal);lam=math.exp(x[2])
        moment=g.terminal(g.bernoulli_epoch(spectrum,kernel,theta,lam),32768)
        value=moment+cutoff*lam+max(corners@(np.array(costs)-128*logs)-127*multinomial)+lattice
        return float(value)/4194304
    proposal=np.array(witness['proposal']);start=np.log(proposal[:2]/proposal[2])
    candidates=[]
    grid=np.linspace(-10.,-.1,25)
    values=[objective(np.r_[start,z]) for z in grid]
    for i in sorted(range(len(grid)),key=lambda i:values[i])[:3]:
        z=grid[i]
        trial=minimize(objective,np.r_[start,z],method='L-BFGS-B',
                       bounds=[(-25.,25.),(-25.,25.),(-12.,0.)],
                       options=dict(maxiter=120,ftol=1e-13,gtol=1e-9))
        candidates.append((objective(trial.x),trial.x,trial.success))
    value,x,success=min(candidates,key=lambda entry:entry[0])
    logits=np.r_[x[:2],0.];proposal=np.exp(logits-float(np.logaddexp.reduce(logits)))
    return dict(log_bound=value*4194304,margin_bits=-value*4194304/math.log(2),
                witness=dict(probability_bank=witness['probability_bank'],proposal=proposal.tolist(),log_surprisal=float(x[2])),
                optimizer_reported_success=bool(success),evaluations=evaluations)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--limit',type=int,default=4)
    parser.add_argument('--map',type=Path)
    args=parser.parse_args()
    spectrum,kernel,sources=g.tv.fixed.load_inner()
    witnesses_path=g.tv.fixed.HERE/'SMALLER_OUTWARD_WITNESSES.json'
    old_path=HERE/'SPARSE_MIXER_R1_Q128_DENSE.json'
    if args.map:
        path=args.map.resolve();record=json.loads(path.read_text())
        rows=g.tv.fixed.maps.generators(record['columns'],19)
        exact=g.tv.fixed.maps.spectrum(rows);dual=g.tv.fixed.maps.dual_spectrum(exact,128,19)
        assert exact=={int(w):int(n) for w,n in record['spectrum'].items()}
        assert dual=={int(w):int(n) for w,n in record['kernel'].items()}
        spectrum={w:n for w,n in exact.items() if w};kernel=[dual.get(w,0) for w in range(129)]
        sources.append(path)
        old_path=HERE/'NO_CONSTANT_SPARSE_MIXER_R1_Q3_DENSE.json'
    witnesses=json.loads(witnesses_path.read_text());old=json.loads(old_path.read_text())
    counts={w:n for w,n in enumerate(g.tv.smaller_outer.spectrum()) if w and n};results=[]
    for row,prior in zip(witnesses['results'],old['results']):
        selected=sorted(range(len(prior['dense_box_log_bounds'])),key=lambda i:-prior['dense_box_log_bounds'][i])[:args.limit]
        for index in selected:
            box=row['dense']['selected_boxes'][index]
            answer=optimize_box(box,row['dense'],counts,spectrum,kernel,row['bad_weight'])
            answer.update(distance_target=row['distance_target'],box_index=index,lower=box['lower'],upper=box['upper'])
            results.append(answer);print(json.dumps(answer),flush=True)
    sources += [Path(__file__),Path(g.__file__),witnesses_path,old_path]
    output=dict(status='BINARY64_SELECTED_DENSE_BOXES_NOT_CERTIFICATE',results=results,
                source_sha256={p.relative_to(g.tv.ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in sources})
    (HERE/f'{"NO_CONSTANT_" if args.map else ""}DENSE_RETUNE_{args.limit}.json').write_text(json.dumps(output,indent=2)+'\n',encoding='utf-8')


if __name__=='__main__':main()
