"""Resume shell optimization in new folders with nearby verified warm bases.

Timeouts remain recorded and do not discard prior certificates or stop the
remaining independent objectives. Every successful solve is rationally audited.
"""
import argparse
import math
from pathlib import Path
import bridge as base
import exact_joint_shell_caps as exact


def run(weights,tag,seconds,warm_name):
    assert tag.isidentifier() and Path(warm_name).name==warm_name
    reference=base.BCH/'generated/shift_rank_oa29_joint'
    model,scales=exact.model_build()
    assert model==base.read(reference/'model.json') and scales==base.read(reference/'scales.json')
    original,_=exact.export(model,scales)
    warm=base.HERE/'generated'/warm_name
    receipt=base.read(warm/'cap.json');assert receipt['rational_primal_dual_checks_passed']
    for name,digest in receipt['local_sha256'].items():assert base.sha(base.HERE/name)==digest
    for w in weights:
        folder=base.HERE/'generated'/f'joint_shell_{tag}_w{w}'
        folder.mkdir(exist_ok=False)
        lp=original.replace(' obj: h_38\n',f' obj: q_{w} + 31 h_{w}\n')
        with (folder/'h_38.lp').open('x') as stream:stream.write(lp)
        with (folder/'warm.bas').open('x') as stream:stream.write((warm/'h_38.bas').read_text())
        print('Exact continuation shell',w,'warm',warm.name,flush=True)
        exact.solve(folder,seconds,True)
        if not (folder/'h_38.sol').exists():
            print('No solution for shell',w,'attempt retained; continuing next objective',flush=True)
            continue
        result=exact.audit(folder,w,model,scales)
        result['local_sha256']={str(p.relative_to(base.HERE)):base.sha(p) for p in
            (Path(__file__),Path(exact.__file__),folder/'h_38.lp',folder/'h_38.sol',folder/'warm.bas')}
        result['outer_sha256']={str(p.relative_to(base.BCH)):base.sha(p) for p in
            (reference/'model.json',reference/'scales.json',reference/'audit.json',
             base.BCH/'code/prepare_bch_shift_rank_oa29_probe.py',base.BCH/'code/verify_scaled_rational_solution.py')}
        base.write_new(folder/'cap.json',result)
        print('Exact continued shell',w,'cap bits',math.log2(result['cap']),flush=True)
        warm=folder


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--weights',nargs='+',type=int,required=True)
    p.add_argument('--tag',required=True);p.add_argument('--seconds',type=int,default=300)
    p.add_argument('--warm',default='joint_shell_exact_w62');a=p.parse_args();run(a.weights,a.tag,a.seconds,a.warm)
