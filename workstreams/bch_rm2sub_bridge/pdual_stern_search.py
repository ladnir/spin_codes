"""Bounded four-row information-set search; only explicit words are certified."""
import argparse
import json
from pathlib import Path
import subprocess
import bridge as base
import pdual_low_weight as low


def run(args):
    assert not args.output.exists()
    checks,rows,values=low.data()
    stdin=''.join(f'{hex(w)} {low.syndrome(w,values)}\n' for w in rows)
    exe=base.HERE/'generated/pdual_stern_search'
    wexe='/mnt/c/'+str(exe)[3:].replace('\\','/')
    completed=subprocess.run(['wsl.exe','-d','Ubuntu-24.04','--',wexe,str(args.seconds),str(args.seed)],
        input=stdin,text=True,capture_output=True,timeout=args.seconds+20,check=True)
    result=json.loads(completed.stdout)
    audit={name:low.verify_word(int(result[key],16)) for name,key in
           (('best','best_word'),('best_nonzero_F7','best_nonzero_F7_word'))}
    assert audit['best']['weight']==result['best_weight']
    assert audit['best_nonzero_F7']['weight']==result['best_nonzero_F7_weight'] and audit['best_nonzero_F7']['F7']!=0
    payload=dict(status='BOUNDED_FOUR_ROW_SEARCH_NOT_PROOF_OF_ABSENCE',seed=args.seed,seconds=args.seconds,
        result=result,audit=audit,weight30_found=audit['best']['weight']==30,
        source_sha256={p.relative_to(base.ROOT).as_posix():base.sha(p) for p in
            (Path(__file__),base.HERE/'pdual_stern_search.cpp',base.HERE/'pdual_low_weight.cpp',
             Path(low.__file__),base.HERE/'pdual_rank_refinement.py',base.BCH/'code/certify_bch_shift_rank.py',
             base.BCH/'code/affine_wambach.py',base.BCH/'code/bch_quotient.py')})
    base.write_new(args.output,payload);print(result,audit,flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--seconds',type=int,default=45)
    p.add_argument('--seed',type=int,default=20260908)
    p.add_argument('--output',type=Path,required=True)
    run(p.parse_args())
