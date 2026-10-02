"""Complete local boundary-cell bounds, not a complete dense-domain cover."""
import argparse
from fractions import Fraction as Q
import json
from pathlib import Path
from time import monotonic

from flint import arb,ctx
from atlas import Model
from mixture import actual_components,G,N
import kernel

SCHEMA='four-bit-complete-boundary-cells-1'


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--distance',default='1/20')
    parser.add_argument('--minimum-groups',type=int,default=65)
    parser.add_argument('--excess',nargs='+',default=['1/64','1/32','1/16'])
    parser.add_argument('--precision',type=int,default=192)
    parser.add_argument('--output',type=Path)
    parser.add_argument('--replay',type=Path)
    args=parser.parse_args()
    if args.precision<128:parser.error('precision >=128 required')
    tilt=[Q(1),Q(1,32),Q(1,64),Q(1,128),Q(1,256)]
    record=None
    if args.replay:
        record=json.loads(args.replay.read_text())
        if record.get('schema')!=SCHEMA:parser.error('unsupported record')
        threshold=record['threshold'];q_min=record['minimum_groups'];tilt=list(map(Q,record['input_tilt']))
    else:
        distance=Q(args.distance)
        if not 0<distance<1:parser.error('distance in (0,1) required')
        threshold=int(distance*N);q_min=args.minimum_groups
    components=actual_components();data=kernel.actual();ctx.prec=args.precision
    model=Model(components,data,q_min,threshold,tilt);model.posterior_shuffle=True
    if record is not None:
        for row in record['cells']:
            if 'witness' not in row:
                print('UNRESOLVED CELL',row['excess'],flush=True);continue
            upper=model.outward(tuple(map(Q,row['cell'])),row['witness'])
            print('REPLAY CELL',row['excess'],'log2 upper',upper.log()/arb(2).log(),flush=True)
    else:
        record=dict(schema=SCHEMA,threshold=threshold,minimum_groups=q_min,
                    input_tilt=list(map(str,tilt)),precision=args.precision,cells=[])
        anchor=[r[0] for r in components].index('0003')
        first=Q(q_min,G)*model.features[anchor][1];epsilon=Q(1,2**32)
        for excess in map(Q,args.excess):
            if excess<=0:parser.error('positive excess required')
            start=monotonic();cell=(first-epsilon,first+excess/G,Q(0),epsilon,Q(0),epsilon,Q(0),epsilon)
            if cell[0]<0:parser.error('cell lower bound must be nonnegative')
            complete=model.boundary(cell)
            row=dict(excess=str(excess),cell=list(map(str,cell)),
                     composition_count=None if complete is None else len(complete))
            result=model.composition_cell_proposal(cell)
            if result is not None:
                _,row['witness']=result
                upper=model.outward(cell,row['witness'])
                row['outward_log2']=str(upper.log()/arb(2).log())
                row['log2_upper']=float(upper.log()/arb(2).log())
            row['seconds']=monotonic()-start;record['cells'].append(row)
            if args.output:
                args.output.parent.mkdir(parents=True,exist_ok=True)
                args.output.write_text(json.dumps(record,indent=2)+'\n')
            print('CELL',str(excess),'compositions',row['composition_count'],
                  'log2 upper',row.get('outward_log2','unresolved'),flush=True)
    print('Complete listed cells only; not a full-code certificate.',flush=True)


if __name__=='__main__':main()
