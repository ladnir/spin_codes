"""Diagnostic only: retune the proof's base activity tilt, not the encoder."""
import argparse
from fractions import Fraction as Q
import json
from pathlib import Path
import scalar_cover as sc
import birth_classes as inner


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--distance',default='9/100')
    parser.add_argument('--minimum-groups',type=int,default=97)
    parser.add_argument('--tilts',nargs='+',default=['1/16','3/32','1/8','3/16','1/4'])
    parser.add_argument('--activities',nargs='+',default=['.15','.2','.25','.3','.35'])
    parser.add_argument('--variance-bins',type=int,default=8)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    components=sc.actual_components(129,Q(21,50),Q(1001,1000),True);data=inner.actual();rows=[]
    for tilt in map(Q,args.tilts):
        model=sc.Model(components,data,int(Q(args.distance)*sc.N),args.minimum_groups,
            tilt=tilt,inner=inner,variance_shuffle=True,variance_bins=args.variance_bins)
        for p in map(Q,args.activities):
            if not 0<p<1:parser.error('activities must lie strictly between zero and one')
            x=tilt*p/(1-p+tilt*p)
            score,witness=model.proposal((x,x))
            row=dict(base_tilt=str(tilt),activity=str(p),x=str(x),proposal=score,witness=witness)
            rows.append(row)
            print(json.dumps({k:v for k,v in row.items() if k!='witness'}),flush=True)
            if args.output:args.output.write_text(json.dumps(dict(diagnostic_only=True,
                distance=args.distance,minimum_groups=args.minimum_groups,variance_bins=args.variance_bins,
                kernel='birth-classes',rows=rows),indent=2)+'\n')


if __name__=='__main__':main()
