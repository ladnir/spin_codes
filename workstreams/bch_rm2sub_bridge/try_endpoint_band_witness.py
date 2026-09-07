"""Bounded outward trial of a witness selected by pure-band diagnostics."""
import argparse
from pathlib import Path
import sys
import time

import certificate_search_core as core
import search_certificates as search

base=core.base


def run(discovery,directory,seconds):
    data=base.read(discovery);core.authenticate(data,base.ROOT)
    core.require(data['status']=='PURE_BAND_BINARY64_DIAGNOSTIC_ONLY','Wrong diagnostic input')
    spec=data['instance'];core.require(spec==core.instance(spec['configuration'],spec['message_exponent']),'Wrong instance')
    core.require(seconds>0 and not directory.exists(),'Use a fresh trial directory and positive budget')
    selected=max(data['results'],key=lambda r:r['pure_proxy_after_label_count_bits'])
    witness=dict(anchor=spec['rows'],tilt=selected['tilt'],
        p=[base.encode(core.F.from_float(r['p'])) for r in selected['bands']])
    task=dict(instance=spec,kind='range',interval=[spec['rows'],spec['rows']],tilt=selected['tilt'],witness=witness,
        discovery_sha256=base.sha(discovery))
    directory.mkdir(parents=True);base.write_new(directory/'task.json',task)
    started=time.monotonic()
    status=search.worker([sys.executable,'-B','certificate_search_backend.py','--task',str(directory/'task.json'),
        '--output',str(directory/'certificate.json')],directory/'producer.log',seconds)
    replay='NOT_RUN';margin=None
    if status=='COMPLETED':
        result=base.read(directory/'certificate.json')['result'];margin=result['diagnostic_margins'][0]
        upper=base.decode(result['bounds'][0]['upper']);remaining=seconds-(time.monotonic()-started)
        if upper<=core.F(1,1<<53) and remaining>0:
            replay=search.worker([sys.executable,'-B','verify_certificate_search.py','--certificate',
                str(directory/'certificate.json'),'--output',str(directory/'replay.json')],directory/'replay.log',remaining)
        else:replay='SKIPPED_WEAK_OR_OUT_OF_BUDGET'
    summary=dict(status=status,replay_status=replay,diagnostic_margin_bits=margin,tilt=selected['tilt'],
        seconds=time.monotonic()-started,full_coverage=False)
    base.write_new(directory/'outcome.json',summary);print(summary,flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--discovery',type=Path,required=True);p.add_argument('--directory',type=Path,required=True)
    p.add_argument('--seconds',type=float,default=240);a=p.parse_args()
    run(a.discovery.resolve(),a.directory.resolve(),a.seconds)
