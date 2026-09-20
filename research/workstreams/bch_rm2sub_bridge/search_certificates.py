"""Serial, resumable, wall-clock-bounded search; numerical failures are unresolved."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys
import time

import certificate_search_core as core
import verify_certificate_search as verifier

base=core.base


def task_key(task):
    return hashlib.sha256(json.dumps(task,sort_keys=True).encode()).hexdigest()[:20]


def worker(command,log,seconds):
    """One child at a time; terminate and reap it before returning on timeout."""
    with log.open('x',encoding='utf-8') as stream:
        child=subprocess.Popen(command,cwd=base.HERE,stdout=stream,stderr=subprocess.STDOUT)
        try:
            code=child.wait(timeout=max(.001,seconds))
            return 'COMPLETED' if code==0 else 'PROCESS_ERROR'
        except subprocess.TimeoutExpired:
            child.kill();child.wait();return 'TIMEOUT'
        except BaseException:
            child.kill();child.wait();raise


def prior_anchors(spec):
    """Hints only: every imported witness must be evaluated on the current instance."""
    records=[]
    for path in sorted((base.HERE/'generated').glob('inner_'+spec['configuration']+'_q*.json')):
        if path.stem.endswith('_replay'): continue
        data=base.read(path)
        if data.get('message_exponent')!=spec['message_exponent']: continue
        if data.get('rows')!=spec['rows'] or data.get('cutoff')!=spec['cutoff']: continue
        core.authenticate(data,base.ROOT)
        records.append(data)
    return records


def propose(spec,best,tried,history):
    def unseen(task):
        task=dict(instance=spec,**task)
        return task if task_key(task) not in tried else None
    if 1 not in best:
        task=unseen(dict(kind='q1',scaled_tilts=[264,295,332,376]))
        if task:return task
    for old in history:
        q=old['occupation']
        if q not in best and base.decode(old['upper'])<core.F(1,1<<40):
            task=unseen(dict(kind='range',interval=[q,q],tilt=old['tilt'],witness=old['witness']))
            if task:return task
    rows=spec['rows']
    low=next((q for q in range(2,rows+1) if q not in best),None)
    if low is None:return None
    prediction=max(-120,min(30,round((-76+6*math.log2(low)-10*(spec['message_exponent']-20)*math.log(2))/5)*5))
    hi=min(rows,max(low+7,int(low*1.22)),low+63)
    first=unseen(dict(kind='range',interval=[low,hi],tilt=prediction))
    if first and len(tried)%2==0:return first
    # Probe the dense endpoint before spending the whole budget on a sparse gap.
    if rows not in best:
        endpoint=[r for r in history if r['occupation']==rows]
        center=max(endpoint,key=lambda r:r['witness']['margin_bits_diagnostic'])['tilt'] if endpoint else 6
        used={r['tilt'] for r in endpoint}
        for delta in (0,-2,2,-4,4,-6,6):
            tilt=center+delta
            if tilt in used or not -120<=tilt<=30:continue
            task=unseen(dict(kind='range',interval=[rows,rows],tilt=tilt))
            if task:return task
    if first:return first
    # Refine a sparse witness, then narrow its interval if necessary.
    for end in dict.fromkeys((hi,(low+hi)//2,low)):
        for delta in (-5,5,-10,10):
            tilt=max(-120,min(30,prediction+delta))
            task=unseen(dict(kind='range',interval=[low,end],tilt=tilt))
            if task:return task
    return None


def collect(directory,spec,bits,reuse):
    best=core.legacy_baseline(spec) if reuse else {}
    core.require(best is not None,'Legacy baseline unavailable')
    accepted=[];tried=set();attempts=[]
    for folder in sorted(directory.glob('job_*')):
        task=base.read(folder/'task.json');core.require(task['instance']==spec,'Wrong cached instance')
        tried.add(task_key(task))
        outcome=base.read(folder/'outcome.json') if (folder/'outcome.json').exists() else dict(status='INTERRUPTED')
        attempts.append(dict(job=folder.name,task=task,outcome=outcome))
        certificate,replay=folder/'certificate.json',folder/'replay.json'
        if not (certificate.exists() and replay.exists()):continue
        values=verifier.load_verified(certificate,replay,spec)
        allowance=core.allocation(best,spec['rows'],bits)
        selected=[r for r in values if (base.decode(r['upper'])<=
            (core.F(1,1<<(bits+1)) if r['occupation']==1 else allowance)) or
            (r['occupation'] in best and base.decode(r['upper'])<best[r['occupation']])]
        if selected:
            core.merge_rows(best,selected,spec['rows'])
            accepted.append(dict(certificate=certificate.relative_to(directory).as_posix(),
                replay=replay.relative_to(directory).as_posix(),occupancies=[r['occupation'] for r in selected]))
    return best,accepted,tried,attempts


def run(directory,name,m,bits,seconds,job_seconds,max_jobs,reuse):
    core.require(math.isfinite(seconds) and math.isfinite(job_seconds) and seconds>0 and job_seconds>0
        and max_jobs>=0,'Invalid work budget')
    core.require(type(bits) is int and 1<=bits<=256,'Invalid target margin')
    directory=directory.resolve();directory.mkdir(parents=True,exist_ok=True)
    spec=core.instance(name,m)
    manifest=dict(schema='bch-search-policy-v1',instance=spec,target_bits=bits,reuse_legacy=reuse)
    path=directory/'manifest.json'
    if path.exists():core.require(base.read(path)==manifest,'Resume policy/instance mismatch')
    else:base.write_new(path,manifest)
    started=time.monotonic();history=prior_anchors(spec)
    for _ in range(max_jobs):
        best,accepted,tried,attempts=collect(directory,spec,bits,reuse)
        if core.summarize(best,spec['rows'],bits)['status']=='CERTIFIED':break
        remaining=seconds-(time.monotonic()-started)
        if remaining<=0:break
        hints=list(history)
        for attempt in attempts:
            path=directory/attempt['job']/'certificate.json'
            old=attempt['task']
            if old['kind']=='range' and old['interval'][0]==old['interval'][1] and path.exists():
                result=base.read(path)['result']
                hints.append(dict(occupation=old['interval'][0],tilt=old['tilt'],
                    witness=result['witness'],upper=result['bounds'][0]['upper']))
        task=propose(spec,best,tried,hints)
        if task is None:break
        folder=directory/f'job_{len(attempts):04d}_{task_key(task)}';folder.mkdir()
        base.write_new(folder/'task.json',task)
        print('trial',name,task['kind'],task.get('interval'),task.get('tilt'),flush=True)
        trial_started=time.monotonic();limit=min(job_seconds,remaining)
        status=worker([sys.executable,'-B','certificate_search_backend.py','--task',str(folder/'task.json'),
            '--output',str(folder/'certificate.json')],folder/'producer.log',limit)
        replay_status='NOT_RUN'
        if status=='COMPLETED':
            data=base.read(folder/'certificate.json');allowance=core.allocation(best,spec['rows'],bits)
            useful=any(base.decode(r['upper'])<=(core.F(1,1<<(bits+1)) if r['occupation']==1 else allowance)
                for r in data['result']['bounds'])
            remaining=min(seconds-(time.monotonic()-started),limit-(time.monotonic()-trial_started))
            if useful and remaining>0:
                replay_status=worker([sys.executable,'-B','verify_certificate_search.py','--certificate',
                    str(folder/'certificate.json'),'--output',str(folder/'replay.json')],folder/'replay.log',remaining)
            elif not useful:replay_status='SKIPPED_WEAK_BOUND'
        base.write_new(folder/'outcome.json',dict(status=status,replay_status=replay_status,
            seconds=time.monotonic()-trial_started))
        print('outcome',status,replay_status,flush=True)
    best,accepted,tried,attempts=collect(directory,spec,bits,reuse)
    report=dict(schema='bch-search-report-v1',instance=spec,target_bits=bits,reused_legacy_certificate=reuse,
        accepted_certificates=accepted,attempts=attempts,coverage=core.summarize(best,spec['rows'],bits),
        run_budget=dict(seconds=seconds,job_seconds=job_seconds,max_jobs=max_jobs),
        seconds=time.monotonic()-started,
        limitations=['Search exhaustion is not a code obstruction.','Legacy coverage is authenticated reuse, not a fresh producer run.',
            'Only Q1 and fixed-weight range backends are enabled; fixed-type refinement is not automated.'])
    index=len(list(directory.glob('report_*.json')))
    output=directory/f'report_{index:04d}.json';base.write_new(output,report)
    verifier.verify_report(output)
    coverage=report['coverage']
    print(json.dumps(dict(report=str(output),status=coverage['status'],
        covered_intervals=coverage['covered_intervals'],uncovered_intervals=coverage['uncovered_intervals'],
        full_margin_bits=coverage['covered_margin_bits'] if coverage['complete_coverage'] else None),indent=2),flush=True)
    return output


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--directory',type=Path,required=True);p.add_argument('--configuration',required=True)
    p.add_argument('--m',type=int,default=20);p.add_argument('--target-bits',type=int,default=40)
    p.add_argument('--seconds',type=float,default=180);p.add_argument('--job-seconds',type=float,default=120)
    p.add_argument('--max-jobs',type=int,default=4);p.add_argument('--reuse-legacy',action='store_true')
    a=p.parse_args();run(a.directory,a.configuration,a.m,a.target_bits,a.seconds,a.job_seconds,a.max_jobs,a.reuse_legacy)
