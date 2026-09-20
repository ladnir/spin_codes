"""Exact full-dual audit of disjoint bases in the restricted expansion map."""
import argparse
from collections import Counter
import json
from pathlib import Path
import subprocess
import sys

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
import bch_model as model


def run_native(exe,a,b,s):
    text=f'{len(a)} {s}\n'+' '.join(map(str,a))+'\n'+' '.join(map(str,b))+'\n'
    process=subprocess.run([str(exe)],input=text,text=True,capture_output=True,check=True)
    return json.loads(process.stdout)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--exe',type=Path,default=HERE/'packing_audit.exe')
    p.add_argument('--output',type=Path,default=HERE/'PACKING_AUDIT.json')
    p.add_argument('--verify',action='store_true')
    args=p.parse_args()
    saved=model.base.read(args.output) if args.verify else None
    if saved:model.authenticate(saved)
    else:assert not args.output.exists()
    engine=model.Engine(20)
    data=run_native(args.exe.resolve(),engine.a_columns,engine.columns,19)
    assert data['states']==(1<<19)-1
    counts=Counter()
    for w,d,e,n in data['groups']:
        assert 0<=d<=w//19 and 0<=e<=(128-w)//19 and n>0
        counts[w]+=n
    assert dict(counts)==engine.b_spectrum
    if saved:
        assert saved['audit']==data
        model.base.write_new(args.output.with_name(args.output.stem+'_replay.json'),
                             dict(status='EXACT_PACKING_AUDIT_REPLAY_PASSED',producer_sha256=model.base.sha(args.output),
                                  executable_sha256=model.base.sha(args.exe)))
    else:
        paths=[Path(__file__),HERE/'packing_audit.cpp',args.exe.resolve(),
               model.ASYMMETRIC.parent/'NO_CONSTANT_MAP.json',model.ASYMMETRIC/'FEEDBACK_SCREEN.json']
        model.base.write_new(args.output,dict(status='EXACT_ALL_DUAL_SUPPORT_BASIS_PACKING',t=128,s=19,
                            policy='first-fit disjoint bases; max forward/reverse column order; low-pivot rank cross-check',
                            audit=data,source_sha256={path.relative_to(model.ROOT).as_posix():model.base.sha(path) for path in paths}))
    print('Audited',data['states'],'dual states in',len(data['groups']),'groups',flush=True)


if __name__=='__main__':main()
