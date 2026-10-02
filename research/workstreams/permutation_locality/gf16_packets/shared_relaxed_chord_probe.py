"""Compare isolated chord candidates with fixed-witness baseline cell bounds.

This regenerates the real shared-route model. Only the listed cells and
their endpoints are covered; this script never produces a whole-code claim.
"""
import argparse
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path

from flint import arb,ctx
import shared_relaxed_chord as chord
import shared_relaxed_strategy as proof


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source',type=Path)
    parser.add_argument('--index',type=int,default=0)
    parser.add_argument('--paths',nargs='+',required=True)
    parser.add_argument('--precision',type=int,default=256)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists() or args.output.resolve()==args.source.resolve():
        parser.error('new output path required; prior evidence is preserved')
    raw=args.source.read_bytes();record=json.loads(raw)
    row=proof.validate_record(record,args.index)
    model=proof.build_model(record,args.precision,args.index)
    cells=proof.sc.partition(model,row['cover']['leaves'],row['cover']['unresolved'])
    if any(path not in cells for path in args.paths):
        parser.error('requested path is not in the saved full partition')
    probes=[v for v in row.get('probes',[]) if 'witness' in v and not v.get('empty',False)]
    if not probes:parser.error('saved point witnesses are needed as fixed-parameter proposals')
    result=dict(schema='shared-gf16-chord-probes-1',source_sha256=hashlib.sha256(raw).hexdigest(),
        updates=2,precision=args.precision,threshold=row['threshold'],minimum_groups=record['minimum_groups'],
        note='Listed cells/endpoints only. Baseline uses the same fixed outer duals and output tilt, '
             'without importing point-specific regional or variance partitions.',rows=[])
    for path in args.paths:
        cell=cells[path];mid=sum(cell)/2
        source=min(probes,key=lambda v:abs(Q(v['mean'])-mid))
        old=source['witness']
        witness={key:old[key] for key in ('tilt','parameters','variance_dual')}
        for name,scope in (('cell',cell),('left',(cell[0],cell[0])),('right',(cell[1],cell[1]))):
            ctx.prec=args.precision
            baseline=model.outward(scope,witness)
            baseline_log=baseline.log()/arb(2).log()
            ctx.prec=args.precision
            upper,details=chord.outward(model,scope,witness)
            chord_log=chord.aq(upper).log()/arb(2).log()
            checked=dict(path=path,kind=name,cell=list(map(str,scope)),witness=witness,
                source_mean=source['mean'],baseline_upper=[int(v) for v in baseline.upper().man_exp()],
                chord_upper=proof.compact_endpoint(upper),baseline_log2=str(baseline_log),
                chord_log2=str(chord_log),details=details)
            result['rows'].append(checked)
            print('CHORD',path,name,'baseline',baseline_log,'new',chord_log,flush=True)
            args.output.parent.mkdir(parents=True,exist_ok=True)
            args.output.write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
