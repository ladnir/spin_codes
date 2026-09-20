"""Practical serial grid triage: shortlist first, optional bounded anchor checks."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import sys
import time

import certificate_search_core as core
import search_certificates as processes
import verify_certificate_search as verifier

base=core.base
CATALOG=tuple(sorted((*core.inputs.maps.NAMES,*base.CONFIGS)))
REFERENCE='t64_s20'


def key(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True).encode()).hexdigest()[:20]


def policy(names,exponents,target,headroom,tolerance,density,tilts,reuse):
    core.require(names and set(names)<=set(CATALOG),'Choose names from --list-configurations')
    core.require(exponents and all(type(m) is int and 13<=m<=20 for m in exponents),'Supported exponents: 13..20')
    core.require(type(target) is int and 1<=target<=256,'Invalid target')
    core.require(math.isfinite(headroom) and headroom>=0 and math.isfinite(tolerance) and tolerance>=0,'Invalid ranking thresholds')
    core.require(density>=3 and tilts>=3,'Screen grids require at least three points')
    return dict(schema='parameter-grid-policy-v1',configurations=sorted(set(names)|{REFERENCE}),
        message_exponents=sorted(set(exponents)),target_bits=target,q1_headroom_bits=headroom,
        max_reference_deficit_per_output_bit=tolerance,density_points=density,tilt_count=tilts,
        fractions=[.0625,.25,.5,1.],reuse_reference=reuse)


def screen_task(spec,settings):
    return dict(stage='screen',instance=spec,target_bits=settings['target_bits'],fractions=settings['fractions'],
        density_points=settings['density_points'],tilt_count=settings['tilt_count'],reuse_reference=settings['reuse_reference'])


def load_job(directory,task):
    folder=directory/'jobs'/key(task)
    if not (folder/'task.json').exists():return dict(state='PENDING')
    core.require(base.read(folder/'task.json')==task,'Job identity collision')
    if not (folder/'outcome.json').exists():return dict(state='INTERRUPTED')
    outcome=base.read(folder/'outcome.json')
    if outcome['status']!='COMPLETED':return dict(state=outcome['status'],outcome=outcome)
    result=base.read(folder/'result.json');core.authenticate(result,base.ROOT)
    core.require(result['task']==task and result['schema']=='parameter-grid-job-v1','Wrong cached job')
    if task['stage']=='anchor' and result['result']['replayed']:
        verifier.load_verified(folder/'certificate.json',folder/'replay.json',task['instance'])
    return dict(state='COMPLETED',result=result['result'],outcome=outcome)


def execute(directory,task,seconds):
    folder=directory/'jobs'/key(task);folder.mkdir(parents=True)
    base.write_new(folder/'task.json',task)
    print('grid job',task['instance']['configuration'],'K exponent',task['instance']['message_exponent'],
        task['stage'],task.get('occupation',''),flush=True)
    started=time.monotonic()
    status=processes.worker([sys.executable,'-B','parameter_grid_worker.py','--task',str(folder/'task.json'),
        '--output',str(folder/'result.json')],folder/'worker.log',seconds)
    base.write_new(folder/'outcome.json',dict(status=status,seconds=time.monotonic()-started,budget_seconds=seconds))


def classify(rows,settings):
    refs={r['message_exponent']:r for r in rows if r['configuration']==REFERENCE and r['screen_state']=='COMPLETED'}
    for row in rows:
        row.update(rank=None,reference_deficit_per_output_bit=None,priority='UNSCREENED',pareto=False)
        if row['screen_state']!='COMPLETED':continue
        if row['configuration']==REFERENCE:row['priority']='REFERENCE';continue
        reference=refs.get(row['message_exponent'])
        if reference is None:row['priority']='NO_REFERENCE';continue
        deficit=(reference['mixed_min_margin_bits']-row['mixed_min_margin_bits'])/row['output_bits']
        row['reference_deficit_per_output_bit']=deficit
        meets_q1=row['q1_margin_bits']>=settings['target_bits']+settings['q1_headroom_bits']
        row['priority']='SHORTLIST_HEURISTIC' if meets_q1 and deficit<=settings['max_reference_deficit_per_output_bit'] else 'WATCH_HEURISTIC'
    for m in settings['message_exponents']:
        eligible=[r for r in rows if r['message_exponent']==m and r['priority'] in ('SHORTLIST_HEURISTIC','WATCH_HEURISTIC')]
        eligible.sort(key=lambda r:(r['priority']!='SHORTLIST_HEURISTIC',r['epoch_updates'],r['state_bits'],
            -r['q1_margin_bits'],r['reference_deficit_per_output_bit'],r['configuration']))
        for i,row in enumerate(eligible):
            row['rank']=i+1
            # A descriptive Pareto flag, not a certificate or measured speed claim.
            row['pareto']=not any(other is not row and
                other['epoch_updates']<=row['epoch_updates'] and other['state_bits']<=row['state_bits'] and
                other['q1_margin_bits']>=row['q1_margin_bits'] and other['mixed_min_margin_bits']>=row['mixed_min_margin_bits'] and
                (other['epoch_updates']<row['epoch_updates'] or other['state_bits']<row['state_bits'] or
                 other['q1_margin_bits']>row['q1_margin_bits'] or other['mixed_min_margin_bits']>row['mixed_min_margin_bits'])
                for other in eligible)
    return rows


def anchor_tasks(spec,settings):
    rows=spec['rows'];m=spec['message_exponent']
    for q in sorted(set((2,max(2,rows//3)))):
        prediction=max(-120,min(30,round((-76+6*math.log2(q)-10*(m-20)*math.log(2))/5)*5))
        for delta in (0,5):
            yield dict(stage='anchor',instance=spec,target_bits=settings['target_bits'],occupation=q,
                tilt=max(-120,min(30,prediction+delta)))


def collect(directory,specs,settings):
    result=[]
    for cell in specs:
        if 'error' in cell:
            result.append(dict(configuration=cell['configuration'],message_exponent=cell['message_exponent'],
                screen_state='UNSUPPORTED_GEOMETRY',error=cell['error']));continue
        spec=cell['instance'];screen=load_job(directory,screen_task(spec,settings))
        row=dict(configuration=spec['configuration'],message_exponent=spec['message_exponent'],
            step_bits=spec['step_bits'],state_bits=spec['state_bits'],output_bits=256*spec['rows'],
            epoch_updates=256*spec['rows']//spec['step_bits'],screen_state=screen['state'],
            instance_fingerprint=spec['fingerprint'],full_margin_bits=None,anchor_evidence=[])
        if screen['state']=='COMPLETED':row.update(screen['result'])
        for task in anchor_tasks(spec,settings):
            record=load_job(directory,task)
            if record['state']!='PENDING':
                row['anchor_evidence'].append(dict(state=record['state'],occupation=task['occupation'],tilt=task['tilt'])
                    |record.get('result',{}))
        row['replayed_occupancies']=sorted(set(a['occupation'] for a in row['anchor_evidence'] if a.get('replayed')))
        row['evidence_level']='FULL_CERTIFICATE_REUSED' if row['full_margin_bits'] is not None else (
            'PARTIAL_REPLAYED_ANCHORS' if row['replayed_occupancies'] else (
                'SCREEN_ONLY' if screen['state']=='COMPLETED' else 'NO_COMPLETED_SCREEN'))
        attempted=set(a['occupation'] for a in row['anchor_evidence'])
        unresolved=attempted-set(row['replayed_occupancies'])
        row['refinement_signal']=('MIXED_EVIDENCE' if unresolved else 'PARTIAL_POSITIVE') if row['replayed_occupancies'] else (
            'UNRESOLVED_PROBES' if attempted else 'NOT_REFINED')
        result.append(row)
    return classify(result,settings)


def select_refinements(rows,top,uncertain):
    ordered=sorted(rows,key=lambda r:(-r['message_exponent'],r.get('rank') or 10**6,r['configuration']))
    return [r for r in ordered if r['priority']=='SHORTLIST_HEURISTIC'][:top]+[
        r for r in ordered if r['priority']=='WATCH_HEURISTIC'][:uncertain]


def write_report(directory,specs,settings,summary):
    rows=collect(directory,specs,settings)
    completed=sum(r['screen_state']=='COMPLETED' for r in rows)
    terminal=sum(r['screen_state']!='PENDING' for r in rows)
    report=dict(schema='parameter-grid-report-v1',settings=settings,rows=rows,
        screening_complete=terminal==len(rows),successful_screens=completed,total_cells=len(rows),
        success_means='A ranked, evidence-labelled grid; not certification of every cell.',**summary)
    existing=list(directory.glob('summary_*.json'));index=max([int(p.stem.split('_')[1]) for p in existing]+[-1])+1
    path=directory/f'summary_{index:04d}.json';base.write_new(path,report)
    lines=['# Parameter-grid shortlist','',
        'Ranks are heuristic exploration priorities, not security margins. Epoch updates are an operation-count proxy, not timings.','',
        '| K exponent | Map | Rank | Priority | Q1 bound bits (unreplayed) | Mixed-screen deficit / output bit | Epoch updates | Evidence |',
        '|---:|---|---:|---|---:|---:|---:|---|']
    for r in sorted(rows,key=lambda r:(r['message_exponent'],r.get('rank') or 0,r['configuration'])):
        q=f"{r['q1_margin_bits']:.3f}" if 'q1_margin_bits' in r else '-'
        d=f"{r['reference_deficit_per_output_bit']:.5f}" if r['reference_deficit_per_output_bit'] is not None else '-'
        lines.append(f"| {r['message_exponent']} | {r['configuration']} | {r['rank'] or '-'} | {r['priority']} | {q} | {d} | {r.get('epoch_updates','-')} | {r.get('evidence_level',r['screen_state'])} |")
    lines+=['','Q1 is only one occupancy. A weak sampled bound is not an obstruction. Full margins appear only in the JSON full_margin_bits field when backed by reused full coverage.',
        '',f"Completed screens: {completed}/{len(rows)}. Jobs executed this invocation: {summary['jobs_executed']}."]
    with path.with_suffix('.md').open('x',encoding='utf-8') as out:out.write('\n'.join(lines)+'\n')
    print('grid report',path,'successful screens',completed,'of',len(rows),flush=True)
    return path


def verify_report(path):
    report=base.read(path);directory=path.parent
    core.require(report['schema']=='parameter-grid-report-v1','Unknown report schema')
    settings=base.read(directory/'manifest.json');specs=base.read(directory/'instances.json')
    core.require(report['settings']==settings,'Report policy mismatch')
    for cell in specs:
        if 'instance' in cell:
            spec=cell['instance']
            core.require(spec==core.instance(spec['configuration'],spec['message_exponent']),'Changed instance')
    rows=collect(directory,specs,settings)
    core.require(report['rows']==rows,'Report ranking/evidence mismatch')
    core.require(report['successful_screens']==sum(r['screen_state']=='COMPLETED' for r in rows)
        and report['total_cells']==len(rows)
        and report['screening_complete']==all(r['screen_state']!='PENDING' for r in rows),'Report coverage mismatch')
    return 'GRID_REPORT_VERIFIED_NOT_A_SECURITY_CERTIFICATE'


def run(directory,settings,seconds,screen_seconds,refine_seconds,top,uncertain,max_jobs):
    core.require(all(math.isfinite(v) and v>0 for v in (seconds,screen_seconds,refine_seconds)),'Positive finite time budgets required')
    core.require(min(top,uncertain,max_jobs)>=0,'Nonnegative job counts required')
    directory=directory.resolve();directory.mkdir(parents=True,exist_ok=True)
    manifest=directory/'manifest.json'
    if manifest.exists():core.require(base.read(manifest)==settings,'Changed grid policy: use a new directory')
    else:base.write_new(manifest,settings)
    specs=[]
    for m in settings['message_exponents']:
        for name in settings['configurations']:
            try:specs.append(dict(instance=core.instance(name,m)))
            except ValueError as e:
                if str(e)!='Complete epochs per region required':raise
                specs.append(dict(configuration=name,message_exponent=m,error=str(e)))
    identity=directory/'instances.json'
    if identity.exists():core.require(base.read(identity)==specs,'Changed selected map/outer inputs')
    else:base.write_new(identity,specs)
    started=time.monotonic();count=0
    def dispatch(task,limit):
        nonlocal count
        if load_job(directory,task)['state']!='PENDING':return True
        remaining=seconds-(time.monotonic()-started)
        if count>=max_jobs or remaining<min(5.,limit):return False
        execute(directory,task,min(limit,remaining));count+=1;return True
    for cell in specs:
        if 'instance' in cell and not dispatch(screen_task(cell['instance'],settings),screen_seconds):break
    rows=collect(directory,specs,settings)
    selected=select_refinements(rows,top,uncertain)
    for row in selected:
        spec=next(c['instance'] for c in specs if 'instance' in c and c['instance']['fingerprint']==row['instance_fingerprint'])
        solved=set(row['replayed_occupancies'])
        for task in anchor_tasks(spec,settings):
            if task['occupation'] in solved:continue
            if not dispatch(task,refine_seconds):break
            record=load_job(directory,task)
            if record['state']=='COMPLETED' and record['result']['replayed']:solved.add(task['occupation'])
    return write_report(directory,specs,settings,dict(jobs_executed=count,seconds=time.monotonic()-started,
        refinement_selection=[dict(configuration=r['configuration'],message_exponent=r['message_exponent']) for r in selected],
        run_budget=dict(seconds=seconds,screen_seconds=screen_seconds,refine_job_seconds=refine_seconds,
            refine_top=top,refine_uncertain=uncertain,max_jobs=max_jobs)))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--list-configurations',action='store_true');p.add_argument('--directory',type=Path)
    p.add_argument('--verify-report',type=Path)
    p.add_argument('--configurations',nargs='+',default=list(CATALOG));p.add_argument('--m',nargs='+',type=int,default=[16,18,20])
    p.add_argument('--target-bits',type=int,default=40);p.add_argument('--q1-headroom',type=float,default=2)
    p.add_argument('--reference-tolerance',type=float,default=.05)
    p.add_argument('--density-points',type=int,default=17);p.add_argument('--tilts',type=int,default=13)
    p.add_argument('--seconds',type=float,default=300);p.add_argument('--screen-seconds',type=float,default=30)
    p.add_argument('--refine-seconds',type=float,default=45);p.add_argument('--refine-top',type=int,default=2)
    p.add_argument('--refine-uncertain',type=int,default=0);p.add_argument('--max-jobs',type=int,default=24)
    p.add_argument('--no-legacy',action='store_true');a=p.parse_args()
    if a.verify_report:print(verify_report(a.verify_report));sys.exit(0)
    if a.list_configurations:print('\n'.join(CATALOG));sys.exit(0)
    if a.directory is None:p.error('--directory required')
    settings=policy(a.configurations,a.m,a.target_bits,a.q1_headroom,a.reference_tolerance,a.density_points,a.tilts,not a.no_legacy)
    run(a.directory,settings,a.seconds,a.screen_seconds,a.refine_seconds,a.refine_top,a.refine_uncertain,a.max_jobs)
