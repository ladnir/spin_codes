"""Run one bounded information-set counterexample search and audit its words."""
import argparse
from pathlib import Path
import subprocess
import json
import bridge as base
import pdual_rank_refinement as rank
from affine_wambach import gf_pow


def data():
    code=rank.code_data(37);g=int(code['generator_hex'],16);gd=int(code['dual_generator_hex'],16)
    checks=[g<<i for i in range(code['dimension'])]
    checks=[w|((w.bit_count()%2)<<255) for w in checks]
    rows=[gd<<i for i in range(code['dual_dimension'])]+[(1<<256)-1]
    values=[gf_pow(2,7*i) for i in range(255)]+[0]
    return checks,rows,values


def syndrome(word,values):
    value=0
    for i,v in enumerate(values):
        if word>>i&1:value ^= v
    return value


def verify_word(word):
    checks,rows,values=data()
    assert 0<word<1<<256 and all((word&r).bit_count()%2==0 for r in checks)
    return dict(weight=word.bit_count(),F7=syndrome(word,values),
                exact_extended_Pdual_membership_verified=True,word_hex=hex(word))


def run(args):
    if args.verify:
        saved=base.read(args.output)
        for name,digest in saved['source_sha256'].items():assert base.sha(base.ROOT/name)==digest
        for name in ('best','best_nonzero_F7'):
            assert verify_word(int(saved['audit'][name]['word_hex'],16))==saved['audit'][name]
        print(saved['audit'],flush=True);return
    assert not args.output.exists()
    checks,rows,values=data()
    stdin=''.join(f'{hex(w)} {syndrome(w,values)}\n' for w in rows)
    exe=base.HERE/'generated/pdual_low_weight_search'
    wexe='/mnt/c/'+str(exe)[3:].replace('\\','/')
    command=['wsl.exe','-d','Ubuntu-24.04','--',wexe,str(args.seconds),str(args.seed)]
    completed=subprocess.run(command,input=stdin,text=True,capture_output=True,timeout=args.seconds+20,check=True)
    result=json.loads(completed.stdout)
    audit={name:verify_word(int(result[key],16)) for name,key in
           (('best','best_word'),('best_nonzero_F7','best_nonzero_F7_word'))}
    assert audit['best']['weight']==result['best_weight']
    assert audit['best_nonzero_F7']['weight']==result['best_nonzero_F7_weight'] and audit['best_nonzero_F7']['F7']!=0
    payload=dict(status='BOUNDED_CODEWORD_SEARCH_NOT_PROOF_OF_ABSENCE',seed=args.seed,seconds=args.seconds,
        result=result,audit=audit,weight30_found=audit['best']['weight']==30,
        source_sha256={p.relative_to(base.ROOT).as_posix():base.sha(p) for p in
            (Path(__file__),base.HERE/'pdual_low_weight.cpp',Path(rank.__file__),
             base.BCH/'code/certify_bch_shift_rank.py',base.BCH/'code/affine_wambach.py',base.BCH/'code/bch_quotient.py')})
    base.write_new(args.output,payload);print(payload['result'],payload['audit'],flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--seconds',type=int,default=45)
    p.add_argument('--seed',type=int,default=20260907)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--verify',action='store_true')
    run(p.parse_args())
