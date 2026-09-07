"""Search-independent numerical replay and exact final coverage verification."""
import argparse
from pathlib import Path
import certificate_search_core as core
import certificate_search_backend as backend

base=core.base


def verify_certificate(path):
    data=base.read(path);core.authenticate(data,base.ROOT)
    core.require(data['schema']=='bch-certificate-search-v1','Unknown certificate schema')
    task=data['task'];spec=task['instance'];old=data['result']
    actual=backend.compute(task,512,old)
    expected=[1] if task['kind']=='q1' else list(range(task['interval'][0],task['interval'][1]+1))
    core.require([r['occupation'] for r in old['bounds']]==expected,'Certificate occupancy gap or duplicate')
    core.require([r['occupation'] for r in actual['bounds']]==expected,'Replay occupancy gap')
    for a,b in zip(actual['bounds'],old['bounds']):
        core.require(0<base.decode(a['upper'])<=base.decode(b['upper']),'Numerical upper-bound replay failed')
    return dict(status='NUMERICAL_512_BIT_REPLAY_PASSED',certificate_sha256=base.sha(path),
        instance_fingerprint=spec['fingerprint'],source_sha256={
            Path(__file__).relative_to(base.ROOT).as_posix():base.sha(Path(__file__))})


def load_verified(path,replay_path,spec):
    data=base.read(path);receipt=base.read(replay_path)
    core.authenticate(data,base.ROOT);core.authenticate(receipt,base.ROOT)
    core.require(data['task']['instance']==spec,'Cross-instance certificate reuse')
    core.require(receipt['status']=='NUMERICAL_512_BIT_REPLAY_PASSED' and
        receipt['certificate_sha256']==base.sha(path) and receipt['instance_fingerprint']==spec['fingerprint'],
        'Missing or mismatched numerical replay')
    expected=[1] if data['task']['kind']=='q1' else list(range(data['task']['interval'][0],data['task']['interval'][1]+1))
    core.require([r['occupation'] for r in data['result']['bounds']]==expected,'Invalid certificate coverage')
    return data['result']['bounds']


def verify_report(path,fresh=False):
    report=base.read(path);spec=report['instance']
    core.require(spec==core.instance(spec['configuration'],spec['message_exponent']),'Wrong report instance')
    best=core.legacy_baseline(spec) if report['reused_legacy_certificate'] else {}
    core.require(best is not None,'No baseline exists for this instance')
    for item in report['accepted_certificates']:
        certificate=path.parent/item['certificate'];replay_path=path.parent/item['replay']
        core.require(certificate.resolve().is_relative_to(path.parent.resolve()) and
            replay_path.resolve().is_relative_to(path.parent.resolve()),'Receipt outside search directory')
        if fresh: verify_certificate(certificate)
        values=load_verified(certificate,replay_path,spec)
        selected=item['occupancies']
        core.require(len(selected)==len(set(selected)) and set(selected)<=set(r['occupation'] for r in values),
            'Invalid selected occupancy subset')
        core.merge_rows(best,[r for r in values if r['occupation'] in selected],spec['rows'])
    result=core.summarize(best,spec['rows'],report['target_bits'])
    core.require(result==report['coverage'],'Report aggregation or status mismatch')
    return dict(status='SEARCH_REPORT_VERIFIED',numerical_replays_rerun=fresh,
        reused_legacy_producers_rerun=False,coverage=result)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--certificate',type=Path);p.add_argument('--report',type=Path)
    p.add_argument('--output',type=Path);p.add_argument('--fresh',action='store_true');a=p.parse_args()
    core.require((a.certificate is None)!=(a.report is None),'Choose certificate or report')
    result=verify_certificate(a.certificate) if a.certificate else verify_report(a.report,a.fresh)
    if a.output: base.write_new(a.output,result)
    print(result['status'],flush=True)
