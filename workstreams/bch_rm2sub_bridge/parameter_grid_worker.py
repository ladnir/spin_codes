"""One in-process grid job; parent controls its hard wall-clock limit."""
import argparse
import math
from pathlib import Path
import time

import certificate_search_core as core
import certificate_search_backend as backend
import verify_certificate_search as verifier

base=core.base


def sources():
    result=backend.sources()
    result[Path(__file__).relative_to(base.ROOT).as_posix()]=base.sha(Path(__file__))
    result[Path(verifier.__file__).relative_to(base.ROOT).as_posix()]=base.sha(Path(verifier.__file__))
    return result


def margin(upper):
    return math.log2(upper.denominator)-math.log2(upper.numerator)


def run(task,output):
    spec=task['instance'];core.require(spec==core.instance(spec['configuration'],spec['message_exponent']),'Wrong grid instance')
    started=time.monotonic()
    if task['stage']=='screen':
        q1=backend.compute(dict(instance=spec,kind='q1'))
        q1_upper=base.decode(q1['bounds'][0]['upper'])
        path=output.parent/'mixed_screen.json'
        qs=sorted(set(max(2,min(spec['rows'],round(spec['rows']*f))) for f in task['fractions']))
        core.inputs.run(spec['configuration'],spec['message_exponent'],qs,path,1e9,
            task['density_points'],task['tilt_count'])
        mixed=base.read(path)
        core.require(mixed['completed_requested_samples'],'Incomplete mixed screen')
        full=None
        if task['reuse_reference'] and spec['configuration']=='t64_s20' and spec['message_exponent']==20:
            full=core.summarize(core.legacy_baseline(spec),spec['rows'],task['target_bits'])
        result=dict(q1_margin_bits=margin(q1_upper),q1_upper=q1['bounds'][0]['upper'],
            q1_evidence='OUTWARD_PRODUCER_WITHOUT_REPLAY',
            mixed_min_margin_bits=min(r['margin_bits'] for r in mixed['results']),
            mixed_samples=mixed['results'],minimum_A_weight=mixed['minimum_A_weight'],
            minimum_kernel_weight=mixed['minimum_kernel_weight'],full_certificate=full,
            full_margin_bits=full['covered_margin_bits'] if full and full['status']=='CERTIFIED' else None)
    else:
        core.require(task['stage']=='anchor','Unknown grid stage')
        q=task['occupation'];certificate=output.parent/'certificate.json';replay=output.parent/'replay.json'
        backend.produce(dict(instance=spec,kind='range',interval=[q,q],tilt=task['tilt']),certificate)
        data=base.read(certificate);upper=base.decode(data['result']['bounds'][0]['upper'])
        useful=upper<=core.F(1,1<<(task['target_bits']+math.ceil(math.log2(spec['rows']))))
        if useful:base.write_new(replay,verifier.verify_certificate(certificate))
        result=dict(occupation=q,tilt=task['tilt'],diagnostic_margin_bits=data['result']['diagnostic_margins'][0],
            upper=data['result']['bounds'][0]['upper'],replayed=useful,
            evidence='REPLAYED_SINGLE_OCCUPANCY' if useful else 'WEAK_UPPER_BOUND_NOT_AN_OBSTRUCTION',
            certificate_sha256=base.sha(certificate),replay_sha256=base.sha(replay) if useful else None)
    base.write_new(output,dict(schema='parameter-grid-job-v1',task=task,result=result,
        seconds=time.monotonic()-started,source_sha256=sources()))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--task',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();run(base.read(a.task),a.output)
