"""Retune regional candidates even when an iid bound wins initially.

Selected-cell diagnostics only. The comparison is freshly regenerated;
the winning candidate is then evaluated with outward arithmetic.
"""
import argparse
from fractions import Fraction as Q
import hashlib
import json
from math import isfinite
from pathlib import Path
from flint import arb,ctx

import pairwise_cover as cover
from pairwise_assemble import validate_geometry
from pairwise_dense import probe_cell
import regional_count as regional
import variance_partition as variance

FACTORS=tuple(map(Q,('1/4','1/2','3/4','7/8','9/8','5/4','3/2','2')))
SCHEMA='pairwise-gf16-regional-retune-1'
STAGES=(
    ('direct',dict(regional_tilted_atom=True,regional_fine_tilts=True,
                   regional_tilted_variance=True,regional_direct_counts=True,regional_exact_zero=True)),
    ('lazy',dict(regional_lazy_density_through=6)),
    ('joint',dict(regional_joint_return_through=3)),
    ('classified',dict(regional_feedback_classes_from=3,regional_feedback_classes_through=32)),
    ('uniform',dict(regional_feedback_uniform_classes=True,regional_feedback_uniform_replace=True)),
)


def search(model,cell,baseline,keep=2,factors=FACTORS,checkpoint=None):
    if type(keep) is not int or not 1<=keep<=len(STAGES) or any(Q(f)<=0 for f in factors):
        raise ValueError('bounded regional shortlist and positive output-tilt factors required')
    best=baseline;trials=[];candidates=[];witness=baseline[1]
    if not isfinite(best[0]):raise ArithmeticError('finite initial proposal required')
    def record(stage,factor,candidate):
        nonlocal best
        score,checked=candidate
        if not isfinite(score):raise ArithmeticError('finite regional proposal required')
        if score<best[0]:best=candidate
        trials.append(dict(stage=stage,factor=str(factor),proposal=float(score)))
        print('PAIRWISE RETUNE',stage,'factor',str(factor),'proposal',score,flush=True)
        if checkpoint:checkpoint(dict(trials=trials,best_proposal=float(best[0]),best_witness=best[1]))
    # Reuse the checked MGF family when only the inner envelope changes.
    for stage,options in STAGES:
        candidate=regional.propose(model,cell,dict(witness,**options))
        record(stage,Q(1),candidate);candidates.append((candidate[0],stage,candidate[1]))
        witness=candidate[1]
    for _,stage,witness in sorted(candidates,key=lambda row:row[0])[:keep]:
        lam=Q(witness['parameters'][0])
        for factor in map(Q,factors):
            trial=dict(witness,parameters=[str(lam*factor),*witness['parameters'][1:]])
            record(stage,factor,regional.propose(model,cell,trial))
    return best,trials


def replay(saved,precision):
    """Fresh selected-cell check; saved numerical bounds are not inputs."""
    if (saved.get('schema')!=SCHEMA or saved.get('ensemble')!=cover.ENSEMBLE
            or not isinstance(saved.get('witness'),dict) or not isinstance(saved.get('cell'),list)
            or len(saved['cell'])!=2 or type(precision) is not int or precision<256):
        raise ValueError('matching selected-cell witness and precision >=256 required')
    cover.validate_parameters(saved.get('parameters'))
    cell=tuple(map(Q,saved['cell']))
    if not 0<cell[0]<=cell[1]<1:raise ValueError('ordered interior cell required')
    model=cover.build_model(precision,saved['parameters']);validate_geometry(model)
    if not model.root[0]<=cell[0]<=cell[1]<=model.root[1]:raise ValueError('cell outside the actual model domain')
    upper=model.outward(cell,saved['witness'])
    if ctx.prec!=precision or not upper>0:raise ArithmeticError('invalid fresh selected-cell bound')
    return dict(schema=SCHEMA+'-replay',ensemble=cover.ENSEMBLE,parameters=saved['parameters'],
        cell=list(map(str,cell)),precision=precision,partial_only=True,
        upper=[int(v) for v in upper.upper().man_exp()],
        note='Fresh model and outward evaluation of the stated cell only; no full certificate.')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('cover',type=Path)
    parser.add_argument('--mean')
    parser.add_argument('--replay',action='store_true',help='Replay a saved retune witness without proposal search')
    parser.add_argument('--radius',default='0')
    parser.add_argument('--keep',type=int,default=2)
    parser.add_argument('--precision',type=int,default=256)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    proposal_path=args.output.with_suffix('.proposals.json')
    if args.output.exists() or proposal_path.exists() or args.precision<256:
        parser.error('new output paths and precision >=256 required')
    raw=args.cover.read_bytes();saved=json.loads(raw)
    if args.replay:
        if args.mean is not None or Q(args.radius)!=0:parser.error('replay takes its exact scope from the saved witness')
        result=replay(saved,args.precision);result['source_sha256']=hashlib.sha256(raw).hexdigest()
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(json.dumps(result,indent=2)+'\n')
        m,e=result['upper'];print('PAIRWISE RETUNE fresh replay log2',arb(m).log()/arb(2).log()+e,flush=True)
        return
    if args.mean is None:parser.error('proposal search requires --mean')
    cover.validate_record(saved)
    cell=probe_cell(args.mean,args.radius)
    model=cover.build_model(args.precision,saved['parameters']);validate_geometry(model)
    if model.empty(cell) or not 0<cell[0]<=cell[1]<1:parser.error('nonempty interior cell required')
    base=variance.propose(model,cell,model.propose_with(cell,model.tilt))
    record=dict(schema=SCHEMA,ensemble=cover.ENSEMBLE,
        parameters=saved['parameters'],cell=list(map(str,cell)),precision=args.precision,
        source_sha256=hashlib.sha256(raw).hexdigest(),partial_only=True,
        baseline_proposal=float(base[0]),factors=list(map(str,FACTORS)),keep=args.keep)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    def save(proposals):
        proposal_path.write_text(json.dumps(dict(record,**proposals,screen_only=True),indent=2)+'\n')
    (score,witness),trials=search(model,cell,base,args.keep,checkpoint=save)
    upper=model.outward(cell,witness)
    if ctx.prec!=args.precision or not upper>0:raise ArithmeticError('invalid outward result')
    record.update(trials=trials,proposal=float(score),witness=witness,
        upper=[int(v) for v in upper.upper().man_exp()],
        note='Only the stated mean cell; a positive log2 upper bound is a failed diagnostic, not a counterexample.')
    args.output.write_text(json.dumps(record,indent=2)+'\n')
    print('PAIRWISE RETUNE outward log2',upper.log()/arb(2).log(),flush=True)


if __name__=='__main__':main()
