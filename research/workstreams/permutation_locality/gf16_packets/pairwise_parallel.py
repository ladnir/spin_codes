"""Parallel continuation of a pairwise dense cover, with one common model.

The parent regenerates and checks the exact comparison once. Workers use
those identical rational components, fresh actual inner maps, and separate
Arb contexts. Every retained leaf is recomputed; no acceptance label is
imported. The separate full assembler must still replay the final cover.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor
from contextlib import contextmanager,redirect_stdout
from fractions import Fraction as Q
import hashlib
import heapq
import json
from math import isfinite
from multiprocessing import get_context
import os
from pathlib import Path
import sys
from flint import arb,ctx

import pairwise_cover as cover
from pairwise_assemble import validate_geometry
from audit_complete import dyadic


def component_digest(components):
    rows=[[str(Q(c)),str(Q(p)),a] for c,p,a in components]
    return hashlib.sha256(json.dumps(rows,separators=(',',':')).encode()).hexdigest()


def worker_model(components,parameters,precision):
    cover.validate_parameters(parameters)
    rows=[(str(i),c,tuple(cover.sc.probabilities(p)),a) for i,(c,p,a) in enumerate(components)]
    data=cover.birth_classes.actual(parameters['updates']);ctx.prec=precision
    model=cover.pairwise_search.Model(rows,data,parameters['threshold'],parameters['minimum_groups'],
        Q(parameters['base_tilt']),inner=cover.birth_classes,variance_shuffle=True,
        variance_bins=parameters['variance_bins'],regional_count=True)
    if model.components!=components:raise ArithmeticError('worker comparison differs from the checked parent model')
    validate_geometry(model)
    return model


def initialize(components,parameters,precision,log_directory):
    global _model,_precision,_digest,_log
    _log=open(Path(log_directory)/f'worker-{os.getpid()}.log','x')
    with redirect_stdout(_log):_model=worker_model(components,parameters,precision)
    _precision=precision;_digest=component_digest(_model.components)


def evaluate(model,job,precision,digest):
    path,cell,target_bits,old=job
    if ctx.prec!=precision:raise ArithmeticError('worker precision changed')
    upper=None;witness=old
    if old is not None:
        value=model.outward(cell,old)
        if not value>0:raise ArithmeticError('positive rechecked bound required')
        if value<arb(2)**-target_bits:
            upper=value;score=float(value.log()/arb(2).log())
    if upper is None:
        # A sufficient candidate needs verification, not another 24 bits of
        # expensive proposal refinement. This never changes acceptance.
        model.proposal_stop_bits=target_bits+2
        score,witness=model.proposal(cell)
        if score<-(target_bits+2):
            upper=model.outward(cell,witness)
            if not 0<upper<arb(2)**-target_bits:
                raise ArithmeticError('fresh outward bound misses the target')
    if ctx.prec!=precision or not isfinite(score):raise ArithmeticError('invalid worker arithmetic state')
    return dict(path=path,cell=list(map(str,cell)),proposal=float(score),witness=witness,
        upper=None if upper is None else [int(v) for v in upper.upper().man_exp()],component_sha256=digest)


def evaluate_worker(job):
    with redirect_stdout(_log):
        result=evaluate(_model,job,_precision,_digest)
        print('PAIRWISE WORKER finished',job[0],result['proposal'],flush=True)
        return result


def cell_target(path,target_bits,aggregate_bits=None):
    """Allocate 2^-aggregate_bits by normalized binary-cell width."""
    if (not isinstance(path,str) or any(c not in '01' for c in path)
            or type(target_bits) is not int or target_bits<40
            or (aggregate_bits is not None and (type(aggregate_bits) is not int or aggregate_bits<40))):
        raise ValueError('binary path and margin >=40 required')
    return target_bits if aggregate_bits is None else aggregate_bits+len(path)


def search(model,record,mapping,workers,max_cells,max_depth,target_bits,checkpoint,aggregate_bits=None):
    """A full partition after every batch, including all unfinished cells."""
    if (type(workers) is not int or not 1<=workers<=4 or type(max_cells) is not int or max_cells<1
            or type(max_depth) is not int or not 1<=max_depth<=64
            or type(target_bits) is not int or target_bits<40):
        raise ValueError('bounded search budget, depth, and margin required')
    cell_target('',target_bits,aggregate_bits)
    cover.validate_record(record)
    cells=cover.sc.partition(model,record['leaves'],record['unresolved'])
    previous=record.get('visited',0)
    if type(previous) is not int or previous<0:raise ValueError('valid previous work count required')
    pending=[(len(path),path,cell,record['leaves'].get(path,{}).get('witness')) for path,cell in cells.items()]
    heapq.heapify(pending);leaves={};unresolved={};visited=0
    digest=component_digest(model.components)
    while pending and visited<max_cells:
        batch=[heapq.heappop(pending) for _ in range(min(workers,len(pending),max_cells-visited))]
        jobs=[(path,cell,cell_target(path,target_bits,aggregate_bits),old) for _,path,cell,old in batch]
        results=list(mapping(evaluate_worker,jobs))
        if len(results)!=len(jobs):raise ArithmeticError('missing worker result')
        for (depth,path,cell,_),result in zip(batch,results):
            if (result.get('path')!=path or result.get('cell')!=list(map(str,cell))
                    or result.get('component_sha256')!=digest or not isinstance(result.get('witness'),dict)
                    or not isfinite(result.get('proposal',float('nan')))):
                raise ArithmeticError('worker scope or comparison mismatch')
            visited+=1
            if result['upper'] is not None:
                if not 0<dyadic(result['upper'])<Q(2)**-cell_target(path,target_bits,aggregate_bits):
                    raise ArithmeticError('invalid worker upper endpoint')
                leaves[path]=dict(witness=result['witness'],proposal=result['proposal'])
            elif depth>=max_depth:
                unresolved[path]=dict(cell=list(map(str,cell)),proposal=result['proposal'],witness=result['witness'])
            else:
                mid=sum(cell)/2
                heapq.heappush(pending,(depth+1,path+'0',(cell[0],mid),None))
                heapq.heappush(pending,(depth+1,path+'1',(mid,cell[1]),None))
        result=dict(leaves=dict(leaves),unresolved={**unresolved,**{p:dict(cell=list(map(str,c))) for _,p,c,_ in pending}},
                    visited=previous+visited,component_sha256=digest)
        cover.sc.partition(model,result['leaves'],result['unresolved'])
        checkpoint(result)
        print('PAIRWISE PARALLEL visited',previous+visited,'accepted',len(leaves),
              'unresolved',len(result['unresolved']),flush=True)
    return result


@contextmanager
def spawn_path():
    original=list(sys.path);sys.path.insert(0,str(Path(__file__).resolve().parent))
    try:yield
    finally:sys.path[:]=original


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('resume',type=Path)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--workers',type=int,choices=range(1,5),default=4)
    parser.add_argument('--max-cells',type=int,default=500)
    parser.add_argument('--max-depth',type=int,default=24)
    parser.add_argument('--target-bits',type=int,default=64)
    parser.add_argument('--aggregate-budget-bits',type=int,
        help='Allocate this dense budget by cell width; overrides fixed per-cell target bits')
    parser.add_argument('--precision',type=int,default=256)
    args=parser.parse_args()
    logs=args.output.with_suffix('.workers');staging=args.output.with_suffix('.writing.json')
    if args.output.exists() or logs.exists() or staging.exists() or args.resume.resolve()==args.output.resolve():
        parser.error('use new output and worker-log paths; previous proofs are preserved')
    if args.max_cells<1 or not 1<=args.max_depth<=64 or args.target_bits<40 or args.precision<256:
        parser.error('valid work budget, depth, margin, and precision >=256 required')
    if args.aggregate_budget_bits is not None and args.aggregate_budget_bits<40:
        parser.error('dense aggregate budget must be at least 40 bits')
    raw=args.resume.read_bytes();record=json.loads(raw);cover.validate_record(record)
    model=cover.build_model(args.precision,record['parameters']);validate_geometry(model)
    logs.mkdir(parents=True)
    def save(result):
        result.update(schema=cover.SCHEMA,ensemble=cover.ENSEMBLE,parameters=record['parameters'],
            precision=args.precision,target_bits=args.target_bits,search_workers=args.workers,
            aggregate_budget_bits=args.aggregate_budget_bits,
            resume_sha256=hashlib.sha256(raw).hexdigest(),
            note='All retained leaves were freshly recomputed with one common exact comparison. '
                 'Work count includes rechecks. Full assembly and sparse regeneration still required.')
        staging.write_text(json.dumps(result,indent=2)+'\n');staging.replace(args.output)
    with spawn_path(),ProcessPoolExecutor(max_workers=args.workers,mp_context=get_context('spawn'),
            initializer=initialize,initargs=(model.components,record['parameters'],args.precision,logs)) as executor:
        result=search(model,record,executor.map,args.workers,args.max_cells,args.max_depth,args.target_bits,save,
            args.aggregate_budget_bits)
    print('PAIRWISE PARALLEL finished, accepted',len(result['leaves']),'unresolved',len(result['unresolved']),flush=True)


if __name__=='__main__':main()
