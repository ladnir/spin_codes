"""Whole-domain dense cover for pairwise4-gf16-r4, with fresh replay.

This covers occupancies 49..2048 only. The sparse prefix is a separate
proof obligation, and successful cells alone are never a full certificate.
"""
import argparse
import json
import hashlib
from fractions import Fraction as Q
from pathlib import Path
from flint import arb,ctx

import pairwise_dense
import scalar_cover as sc
import birth_classes
import pairwise_search

SCHEMA='pairwise-gf16-dense-cover-1'
ENSEMBLE='pairwise4-gf16-r4'
PARAMETERS=dict(step=8,zero_bits=80,cost_tilt='1/4',mass_bits=None,
                refined_counts=True,coupled_counts=True,pair_shells=True,
                central_bits=260,maximum_activity='7/8',transform_shells=True,
                threshold=209715,updates=4,minimum_groups=49,
                base_tilt='3/16',variance_bins=16,prune_mixture=False,union_cost_tilt=None)


def validate_parameters(parameters):
    if not isinstance(parameters,dict):
        raise ValueError('model parameters required')
    expected=dict(PARAMETERS,mass_bits=parameters.get('mass_bits'),
                  prune_mixture=parameters.get('prune_mixture',False),
                  variance_bins=parameters.get('variance_bins'),
                  union_cost_tilt=parameters.get('union_cost_tilt'))
    if 'prune_mixture' not in parameters:expected.pop('prune_mixture')
    if 'union_cost_tilt' not in parameters:expected.pop('union_cost_tilt')
    if (parameters!=expected or parameters['mass_bits'] not in (None,260)
            or parameters['variance_bins'] not in (16,32,64)
            or type(parameters.get('prune_mixture',False)) is not bool
            or parameters.get('union_cost_tilt') not in (None,'3/16')
            or (parameters.get('union_cost_tilt') is not None and not parameters.get('prune_mixture',False))):
        raise ValueError('matching pairwise proof model required')


def validate_record(record):
    """Reject another ensemble or claim before using saved witnesses."""
    if (record.get('schema')!=SCHEMA or record.get('ensemble')!=ENSEMBLE
            or not isinstance(record.get('leaves'),dict)
            or not isinstance(record.get('unresolved'),dict)):
        raise ValueError('matching pairwise ensemble, parameters and partition required')
    validate_parameters(record.get('parameters'))


def build_model(precision,parameters=PARAMETERS):
    if type(precision) is not int or precision<128:
        raise ValueError('at least 128-bit directed precision required')
    validate_parameters(parameters)
    p=parameters
    rows,_,_=pairwise_dense.components(p['step'],p['zero_bits'],Q(p['cost_tilt']),
        p['mass_bits'],p['refined_counts'],p['coupled_counts'],p['pair_shells'],
        p['central_bits'],Q(p['maximum_activity']),p['transform_shells'],p.get('prune_mixture',False),
        union_cost_tilt=p.get('union_cost_tilt'))
    data=birth_classes.actual(p['updates']);ctx.prec=precision
    return pairwise_search.Model(rows,data,p['threshold'],p['minimum_groups'],Q(p['base_tilt']),
        inner=birth_classes,variance_shuffle=True,variance_bins=p['variance_bins'],
        regional_count=True)


def replay(model,record,target_bits=40):
    validate_record(record)
    if type(target_bits) is not int or target_bits<40:
        raise ValueError('at least 40 aggregate bits required')
    cells=sc.partition(model,record['leaves'],record['unresolved'])
    total=arb(0)
    for index,(path,row) in enumerate(record['leaves'].items(),1):
        upper=model.outward(cells[path],row['witness'])
        if not upper>0:
            raise ArithmeticError('positive outward cell bound required')
        total=sc.kernel.up(total+upper)
        print('PAIRWISE DENSE replay cell',index,'of',len(record['leaves']),flush=True)
    complete=not record['unresolved'] and 0<total<arb(2)**-target_bits
    return dict(upper=[int(v) for v in total.upper().man_exp()],
                unresolved=len(record['unresolved']),dense_complete=bool(complete),
                note='Dense occupancies 49..2048 only; sparse prefix must be added separately.')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--max-cells',type=int,default=200)
    parser.add_argument('--max-depth',type=int,default=24)
    parser.add_argument('--target-bits',type=int,default=64)
    parser.add_argument('--aggregate-bits',type=int,default=48)
    parser.add_argument('--precision',type=int,default=256)
    parser.add_argument('--mass-bits',type=int,choices=(260,),help='Optional bound on each pair comparison component')
    parser.add_argument('--prune-mixture',action='store_true',default=None)
    parser.add_argument('--variance-bins',type=int,choices=(16,32,64))
    parser.add_argument('--union-cost-tilt',choices=('3/16',),help='Use the pair-union mass gradient when pruning the comparison')
    parser.add_argument('--output',type=Path,required=True)
    source=parser.add_mutually_exclusive_group()
    source.add_argument('--resume',type=Path)
    source.add_argument('--replay',type=Path)
    source.add_argument('--partition-from',type=Path,help='Reuse only subdivision paths with a new comparison; repropose and recheck every cell')
    parser.add_argument('--repropose-all',action='store_true')
    args=parser.parse_args()
    if (args.max_cells<1 or not 1<=args.max_depth<=64 or args.target_bits<40 or args.aggregate_bits<40
            or args.precision<128 or (args.repropose_all and not args.resume)):
        parser.error('valid work limits and directed precision required')
    path=args.resume or args.replay or args.partition_from
    if path and path.resolve()==args.output.resolve():
        parser.error('preserve the input checkpoint with a separate output path')
    saved=json.loads(path.read_text()) if path else None
    if saved is not None:validate_record(saved)
    new_model=saved is None or args.partition_from
    parameters=dict(PARAMETERS,mass_bits=args.mass_bits,prune_mixture=bool(args.prune_mixture),
                    variance_bins=args.variance_bins or 16,union_cost_tilt=args.union_cost_tilt) if new_model else saved['parameters']
    if not new_model and args.mass_bits is not None and parameters['mass_bits']!=args.mass_bits:
        parser.error('a saved model cannot be changed during resume or replay')
    if not new_model and args.prune_mixture and not parameters.get('prune_mixture',False):
        parser.error('a saved model cannot be pruned during resume or replay')
    if not new_model and args.variance_bins is not None and parameters['variance_bins']!=args.variance_bins:
        parser.error('a saved variance partition setting cannot be changed during resume or replay')
    if not new_model and args.union_cost_tilt is not None and parameters.get('union_cost_tilt')!=args.union_cost_tilt:
        parser.error('a saved comparison objective cannot be changed during resume or replay')
    model=build_model(args.precision,parameters)
    print('PAIRWISE DENSE domain',tuple(map(str,model.root)),flush=True)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    if args.replay:
        result=replay(model,saved,args.aggregate_bits)
        result.update(schema=SCHEMA+'-replay',ensemble=ENSEMBLE,
                      parameters=parameters,precision=args.precision,aggregate_bits=args.aggregate_bits)
        args.output.write_text(json.dumps(result,indent=2)+'\n')
        print('PAIRWISE DENSE fresh replay',result,flush=True)
        return
    def save(result):
        result.update(schema=SCHEMA,ensemble=ENSEMBLE,parameters=parameters,
                      precision=args.precision,target_bits=args.target_bits,
                      note='Partial dense cover unless no unresolved cells remain; fresh aggregate replay required.')
        if args.partition_from:
            result['partition_source_sha256']=hashlib.sha256(path.read_bytes()).hexdigest()
        args.output.write_text(json.dumps(result,indent=2)+'\n')
    result=sc.run(model,args.max_cells,args.max_depth,args.target_bits,save,saved,
                  coalesce=False,repropose=args.repropose_all or bool(args.partition_from))
    sc.partition(model,result['leaves'],result['unresolved'])
    save(result)
    print('PAIRWISE DENSE search finished, accepted',len(result['leaves']),
          'unresolved',len(result['unresolved']),flush=True)


if __name__=='__main__':main()
