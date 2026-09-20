"""Fixed-budget weight-40/42 experiments with frozen inputs and full replay.

No adaptive retries, early acceptance, or use of the preflight as a sample
prefix. A nonzero hit count gives an inconclusive outcome for that shell.
"""
from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone
from fractions import Fraction
from pathlib import Path

sys.dont_write_bytecode = True
from analyze_endpoint_sampling import log_integer, neg_log_one_minus, ratio_record
from run_endpoint_certificate import ExperimentLock
from run_higher_endpoint_preflight import (EXE, PROTOCOL, ROOT, frozen_inputs,
                                          sha, verify as verify_preflight,
                                          verify_trace, write_new)
from verify_higher_shell_endpoint_reduction import sample_plans

PREFLIGHT = ROOT/'generated/higher_endpoint_preflight_20260904'


def inputs():
    result = frozen_inputs()
    paths = [Path(__file__), PREFLIGHT/'manifest.json', PREFLIGHT/'preflight.json',
             PREFLIGHT/'vectors.txt', PREFLIGHT/'vectors.json', PREFLIGHT/'kernel_check.txt']
    result.update({str(p.relative_to(ROOT)).replace('\\','/'):sha(p) for p in paths})
    return result


def fraction(record):
    return Fraction(int(record['numerator']), int(record['denominator']))


def verify_shell(directory, plan):
    weight, n = plan['weight'], plan['fixed_zero_hit_sample_budget']
    primary = json.loads((directory/'sample.json').read_text())
    replay = json.loads((directory/'replay.json').read_text())
    assert primary['mode']=='BCryptGenRandom' and replay['mode']=='independent_replay'
    keys = ('weight','samples','endpoint_hits','singular_queries','bytes_consumed',
            'bytes_loaded','zero_bytes_rejected','duplicate_bytes_rejected','agreement_histogram')
    assert all(primary[k]==replay[k] for k in keys), 'Replay disagrees with primary'
    assert primary['weight']==weight and primary['samples']==n
    for key in keys[:-1]:
        assert isinstance(primary[key],int) and primary[key]>=0
    histogram=primary['agreement_histogram']
    assert len(histogram)==weight and all(isinstance(v,int) and v>=0 for v in histogram)
    assert sum(histogram)+primary['singular_queries']==n
    assert histogram[weight-1]==primary['endpoint_hits'] and sum(histogram[:weight-20])==0
    assert primary['bytes_consumed']==(weight-19)*n+primary['zero_bytes_rejected']+primary['duplicate_bytes_rejected']
    size=(directory/'random.bin').stat().st_size
    assert size==primary['bytes_loaded'] and 0<=size-primary['bytes_consumed']<(1<<20)
    assert sha(directory/'sample.json.audit.jsonl')==sha(directory/'replay.json.audit.jsonl')
    stage=json.loads((directory/'sample.receipt.json').read_text())
    assert stage['weight']==weight and stage['samples']==n
    for name,expected in stage['sha256'].items():
        assert sha(directory/name)==expected, 'Primary artifacts changed during replay'
    trace=verify_trace(directory,weight,n)
    assert trace['endpoint_records']==primary['endpoint_hits'], 'Every endpoint must have an audit record'
    p_bad=fraction(plan['bad_shell_probability_lower'])
    decay_lo,decay_hi=neg_log_one_minus(p_bad)
    log2_lo,log2_hi=log_integer(2)
    assert n*decay_lo>=41*log2_hi
    assert (n-1)*decay_hi<41*log2_lo
    accepted=primary['endpoint_hits']==0
    return dict(status='conditional_statistical_acceptance' if accepted else 'inconclusive',
        weight=weight,target_cap=plan['target_cap'],samples=n,endpoint_hits=primary['endpoint_hits'],
        singular_queries=primary['singular_queries'],
        full_independent_algorithm_replay_matches=True,independent_python_trace=trace,
        sample_seconds=primary['seconds'],replay_seconds=replay['seconds'],
        random_tape_bytes=size,bytes_consumed=primary['bytes_consumed'],
        false_accept_bound_under_ideal_IID='2^-41' if accepted else None,
        no_hit_log2_probability_upper=ratio_record(-n*decay_lo/log2_hi),
        sha256={name:sha(directory/name) for name in ('random.bin','sample.json','replay.json',
            'sample.json.audit.jsonl','replay.json.audit.jsonl','sample.receipt.json')})


def verify(directory):
    manifest=json.loads((directory/'manifest.json').read_text())
    assert manifest['protocol']==PROTOCOL and manifest['diagnostic_only'] is False
    assert manifest['input_sha256']==inputs(), 'A frozen input has changed'
    assert manifest['plans']==sample_plans(), 'Plan changed'
    results={str(p['weight']):verify_shell(directory/f"w{p['weight']}",p) for p in manifest['plans']}
    accepted=all(r['status']=='conditional_statistical_acceptance' for r in results.values())
    return dict(status='conditional_statistical_acceptance' if accepted else 'inconclusive',
        protocol=PROTOCOL,manifest_sha256=sha(directory/'manifest.json'),results=results,
        familywise_false_accept_bound_for_these_two_tests_under_IID='2^-40',
        familywise_false_accept_bound_including_prior_A38_test_under_IID='2^-39',
        randomness_qualification='Actual bytes are Windows CNG cryptographic pseudorandomness. The statistical bounds assume ideal IID bytes. Under computational replacement, add the applicable distinguishing advantage. No unconditional information-theoretic randomness claim.',
        scope='Two fixed-code shell-count tests only. These are not posterior probabilities, deterministic spectrum bounds, or full-SPIN setup guarantees.',
        next_application_step='Aggregate the accepted shell caps with directed transfer coefficients; retain separate obligations for any omitted shell and occupations Q>=2.')


def progress(sub, mode, start):
    tape=sub/'random.bin'
    audit=sub/('sample.json.audit.jsonl' if mode=='run' else 'replay.json.audit.jsonl')
    last=None
    if audit.exists():
        with audit.open('rb') as stream:
            stream.seek(max(0,audit.stat().st_size-16384))
            tail=stream.read().splitlines()
        for line in reversed(tail):
            try:
                last=json.loads(line)['sample']+1
                break
            except (json.JSONDecodeError,KeyError,UnicodeDecodeError):
                pass
    print(json.dumps(dict(stage=mode,weight=sub.name,elapsed_seconds=round(time.monotonic()-start,1),
                          last_flushed_audit_sample=last,tape_bytes=tape.stat().st_size if tape.exists() else 0)),flush=True)


def run_stage(sub,mode,plan):
    name='sample.json' if mode=='run' else 'replay.json'
    assert not (sub/name).exists() and not (sub/(name+'.audit.jsonl')).exists()
    if mode=='run':
        assert not (sub/'random.bin').exists()
    start=time.monotonic()
    child=subprocess.Popen([str(EXE),mode,str(plan['weight']),str(plan['fixed_zero_hit_sample_budget']),
                            str(sub/'random.bin'),str(sub/name),PROTOCOL],
                           creationflags=subprocess.CREATE_NO_WINDOW)
    try:
        while True:
            try:
                code=child.wait(timeout=30)
                break
            except subprocess.TimeoutExpired:
                progress(sub,mode,start)
        if code:
            raise RuntimeError(f'{mode} failed, exit {code}; artifacts retained, no automatic retry')
    except BaseException:
        if child.poll() is None:
            child.terminate()
            child.wait()
        raise


def execute(directory):
    with ExperimentLock():
        assert shutil.disk_usage(ROOT).free>=20*(1<<30), 'Need at least 20 GiB free for retained tapes and reserve'
        assert verify_preflight(PREFLIGHT)['status']=='validated_preflight_only'
        directory.mkdir(parents=True,exist_ok=False)
        plans=sample_plans()
        manifest=dict(protocol=PROTOCOL,created_utc=datetime.now(timezone.utc).isoformat(),
            diagnostic_only=False,plans=plans,input_sha256=inputs(),
            execution_order=['w40 primary','w40 replay','w40 verification','w42 primary','w42 replay','w42 verification'],
            acceptance_rule='For each weight, exactly zero endpoint hits at the frozen sample count, plus successful replay and verification; otherwise inconclusive. No automatic retry.',
            randomness='Fresh Windows CNG bytes, full 1 MiB chunks retained; preflight tapes never used as prefixes',
            sigma_sampling='One fresh byte, masked with 31, giving the 32-element binary subspace S={0,...,31}',
            subset_sampling='Read bytes, rejecting zero and within-base duplicates until w-20 distinct nonzero elements; singular systems count as non-hits and are never resampled',
            error_model='Ideal IID bytes for the exact binomial bounds; actual CNG requires computational randomness qualification',
            all_endpoint_records_audited=True)
        write_new(directory/'manifest.json',manifest)
        print('FROZEN',directory/'manifest.json',sha(directory/'manifest.json'),flush=True)
        for plan in plans:
            sub=directory/f"w{plan['weight']}"
            sub.mkdir()
            run_stage(sub,'run',plan)
            write_new(sub/'sample.receipt.json',dict(weight=plan['weight'],samples=plan['fixed_zero_hit_sample_budget'],
                sha256={name:sha(sub/name) for name in ('random.bin','sample.json','sample.json.audit.jsonl')}))
            run_stage(sub,'replay',plan)
            assert inputs()==manifest['input_sha256']
            result=verify_shell(sub,plan)
            write_new(sub/'certificate.json',result)
            print('SHELL VERIFIED',json.dumps(result),flush=True)
        result=verify(directory)
        write_new(directory/'certificate.json',result)
        print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('operation',choices=('plan','execute','verify'))
    parser.add_argument('directory',type=Path,nargs='?')
    args=parser.parse_args()
    if args.operation=='plan':
        print(json.dumps(sample_plans(),indent=2))
    else:
        assert args.directory is not None
        if args.operation=='execute':
            execute(args.directory.resolve())
        else:
            print(json.dumps(verify(args.directory.resolve()),indent=2))
