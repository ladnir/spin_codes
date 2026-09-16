"""Discover and outward-check the IMT 11% cover using the refined BA bound."""
import argparse
from collections import Counter
from fractions import Fraction as F
import json
import importlib.util
import math
from pathlib import Path
import sys

from flint import ctx

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
import screen
import certify_imt_dense as dense

_spec=importlib.util.spec_from_file_location('imt_cover_discovery',HERE.parent/'cover.py')
cover=importlib.util.module_from_spec(_spec)
sys.modules[_spec.name]=cover
_spec.loader.exec_module(cover)


def authenticate(record):
    for name,digest in record['source_sha256'].items():
        assert screen.sha(screen.ROOT/name)==digest,name


def sources(extra):
    paths={Path(__file__),HERE.parent.parent/'NO_CONSTANT_MAP.json'}
    paths.update(p.resolve() for p in extra)
    paths.update(Path(mod.__file__).resolve() for mod in list(sys.modules.values())
                 if getattr(mod,'__file__',None) and Path(mod.__file__).suffix=='.py'
                 and Path(mod.__file__).resolve().is_relative_to(screen.ROOT))
    return {p.relative_to(screen.ROOT).as_posix():screen.sha(p) for p in paths}


def discover(segments):
    model=screen.Model()
    prior=json.loads((HERE.parent/'COVER_D11_BOUNDED.json').read_text())
    authenticate(prior)
    bank=cover.Witnesses(model,.11)
    for box in prior['leaves']:
        bank.add(box['witness'])
        if 'point_diagnostic' in box:bank.add(box['point_diagnostic'])
    pending=[]
    for box in prior['leaves']:
        x0,x1=map(F,box['row_density'])
        for segment in segments:
            lo,hi=max(x0,segment.lower),min(x1,segment.upper)
            if lo<hi:
                pending.append(dict(segment=segment.index,alpha=box['alpha'],
                                    row_density=[str(lo),str(hi)],depth=0))
    dense.check_geometry(pending,segments)
    accepted=[];splits=0
    while pending:
        box=pending.pop()
        # Exact uniform-input moment; binary64 only proposes its witnesses.
        al,ah=map(F,box['alpha']);xl,xh=map(F,box['row_density'])
        for a in (al,(al+ah)/2,ah):
            for x in (xl,(xl+xh)/2,xh):
                a,x=float(a),float(x)
                p=min(1-1e-12,.5+.5*a*(1-x)/(1-a*x))
                bank.add(dict(p=p,y=.5/p,lam=math.log(.89/.11),family='statefree'))
        value,w,vertices=bank.best(box,segments)
        if value < -1e-10:
            accepted.append(dict(**box,exponent=value,witness=w,vertices=vertices))
            continue
        assert box['depth']<24,('discovery unresolved',box,value)
        options=[]
        for axis,lo,hi in (('alpha',al,ah),('row_density',xl,xh)):
            mid=(lo+hi)/2
            children=[dict(box,**{axis:[str(lo),str(mid)]},depth=box['depth']+1),
                      dict(box,**{axis:[str(mid),str(hi)]},depth=box['depth']+1)]
            options.append((max(bank.best(b,segments)[0] for b in children),axis,children))
        pending.extend(min(options,key=lambda item:item[0])[2]);splits+=1
        assert splits<10000,'bounded discovery budget'
    dense.check_geometry(accepted,segments)
    return dict(status='BINARY64_REFINED_IMT_D11_COVER_NOT_CERTIFICATE',delta='11/100',
                leaves=accepted,boxes=len(accepted),unresolved=0,splits=splits,
                transfer_counts=dict(Counter(b['witness']['family'] for b in accepted)),
                maximum_exponent=max(b['exponent'] for b in accepted))


def main():
    if not __debug__:raise RuntimeError('Do not use -O.')
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode',choices=['cover','certify','replay'],required=True)
    parser.add_argument('--outer',type=Path,default=HERE/'OUTER_REFINED.json')
    parser.add_argument('--input',type=Path)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists():raise FileExistsError(args.output)
    outer=json.loads(args.outer.read_text());authenticate(outer)
    assert outer['status']=='proved'
    segments=screen.load_segments(args.outer)
    if args.mode=='cover':
        result=discover(segments)
        inputs=[args.outer,HERE.parent/'COVER_D11_BOUNDED.json']
    else:
        assert args.input is not None
        saved=json.loads(args.input.read_text());authenticate(saved)
        assert F(saved['delta'])==F(11,100)
        assert saved['outer_sha256']==screen.sha(args.outer)
        ctx.prec=512 if args.mode=='replay' else 256
        leaves=saved['leaves'];dense.check_geometry(leaves,segments)
        checker=dense.Checker();checked=[]
        if args.mode=='replay':assert saved['inner']==checker.engine.identity()['inner']
        else:assert saved['unresolved']==0
        for i,b in enumerate(leaves):
            box={k:b[k] for k in ('segment','alpha','row_density')}
            box['rational_witness']=b['rational_witness'] if args.mode=='replay' else checker.propose(b['witness'])
            val=checker.box(box,segments,F(11,100))
            assert val<0,(i,float(val))
            if args.mode=='replay':assert val<=F(b['exponent_upper'])
            checked.append(dict(**box,exponent_upper=str(val)))
            if (i+1)%100==0:print(ctx.prec,'checked',i+1,'/',len(leaves),flush=True)
        result=dict(status='OUTWARD_REFINED_IMT_D11_DENSE' if args.mode=='certify' else 'REPLAYED_REFINED_IMT_D11_DENSE',
                    delta='11/100',precision_bits=ctx.prec,alpha_range=['1/10000','1'],
                    row_density_range=['13/125','112/125'],exact_geometry_checked=True,
                    leaves=checked,boxes=len(checked),inner=checker.engine.identity()['inner'],
                    maximum_exponent_upper=str(max(F(b['exponent_upper']) for b in checked)),
                    full_asymptotic_theorem=False)
        inputs=[args.outer,args.input]
    result['outer_sha256']=screen.sha(args.outer)
    result['source_sha256']=sources(inputs)
    with args.output.open('x',encoding='utf-8') as stream:
        json.dump(result,stream,indent=2);stream.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('leaves','source_sha256','inner')}),flush=True)


if __name__=='__main__':main()
