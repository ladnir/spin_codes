"""Bounded sparse-only scheduling around the frozen certificate backend.

Dense endpoints are handled by the explicit split controller, never pooled
into this scheduler's task queue. Producer/replay jobs run strictly serially.
"""
import argparse
import math
from pathlib import Path
import sys
import time

import search_certificates as old

core, base, verifier = old.core, old.base, old.verifier


def propose(spec, best, tried, last):
    def unseen(task):
        task = dict(instance=spec, **task)
        return None if old.task_key(task) in tried else task
    if 1 not in best:
        return unseen(dict(kind='q1', scaled_tilts=[264, 295, 332, 376]))
    low = next((q for q in range(2, last+1) if q not in best), None)
    if low is None:
        return None
    prediction = max(-120, min(30, round((-76+6*math.log2(low)-10*(spec['message_exponent']-20)*math.log(2))/5)*5))
    hi = min(last, max(low+7, int(low*1.22)), low+63)
    for end in dict.fromkeys((hi, (low+hi)//2, low)):
        for delta in (0, -5, 5, -10, 10):
            task = unseen(dict(kind='range', interval=[low, end], tilt=max(-120, min(30, prediction+delta))))
            if task:
                return task
    return None


def interrupted_retry(attempts, tried):
    """One fresh producer/replay retry per interrupted logical task."""
    for attempt in attempts:
        task, outcome = attempt['task'], attempt['outcome']
        if 'retry_of' in task:
            continue
        if outcome.get('status') not in ('TIMEOUT', 'INTERRUPTED') and outcome.get('replay_status') != 'TIMEOUT':
            continue
        retry = dict(task, retry_of=attempt['job'])
        if old.task_key(retry) not in tried:
            return retry
    return None


def run(directory, m=18, last=479, seconds=240, jobs=40, job_seconds=40):
    core.require(seconds > 0 and jobs >= 0 and job_seconds > 0, 'Invalid budget')
    spec = core.instance('t128_s19', m)
    core.require(1 <= last < spec['rows'], 'Sparse endpoint must exclude dense endpoint')
    directory = directory.resolve()
    directory.mkdir(parents=True, exist_ok=True)
    manifest = dict(schema='sparse-only-search-v1', instance=spec, last=last, target_bits=40)
    path = directory/'manifest.json'
    if path.exists():
        core.require(base.read(path) == manifest, 'Changed sparse search instance/policy')
    else:
        base.write_new(path, manifest)
    best, accepted, tried, attempts = old.collect(directory, spec, 40, False)
    started = time.monotonic()
    for _ in range(jobs):
        remaining = seconds-(time.monotonic()-started)
        task = interrupted_retry(attempts, tried) or propose(spec, best, tried, last)
        if remaining <= 0 or task is None:
            break
        folder = directory/f'job_{len(attempts):04d}_{old.task_key(task)}'
        folder.mkdir()
        base.write_new(folder/'task.json', task)
        print('sparse trial', task['kind'], task.get('interval'), task.get('tilt'), flush=True)
        trial_start = time.monotonic()
        limit = min(job_seconds, remaining)
        status = old.worker([sys.executable, '-B', 'certificate_search_backend.py', '--task', str(folder/'task.json'),
            '--output', str(folder/'certificate.json')], folder/'producer.log', limit)
        replay_status = 'NOT_RUN'
        if status == 'COMPLETED':
            data = base.read(folder/'certificate.json')
            allowance = core.allocation(best, spec['rows'], 40)
            selected = [v for v in data['result']['bounds'] if base.decode(v['upper']) <=
                (core.F(2)**-41 if v['occupation'] == 1 else allowance)]
            remaining = min(seconds-(time.monotonic()-started), limit-(time.monotonic()-trial_start))
            if selected and remaining > 0:
                replay_status = old.worker([sys.executable, '-B', 'verify_certificate_search.py', '--certificate',
                    str(folder/'certificate.json'), '--output', str(folder/'replay.json')], folder/'replay.log', remaining)
            elif not selected:
                replay_status = 'SKIPPED_WEAK_BOUND'
            if replay_status == 'COMPLETED':
                verifier.load_verified(folder/'certificate.json', folder/'replay.json', spec)
                core.merge_rows(best, selected, spec['rows'])
                accepted.append(dict(certificate=(folder/'certificate.json').relative_to(directory).as_posix(),
                    replay=(folder/'replay.json').relative_to(directory).as_posix(), occupancies=[v['occupation'] for v in selected]))
        outcome = dict(status=status, replay_status=replay_status, seconds=time.monotonic()-trial_start)
        base.write_new(folder/'outcome.json', outcome)
        attempts.append(dict(job=folder.name, task=task, outcome=outcome))
        tried.add(old.task_key(task))
        print(status, replay_status, 'coverage', core.intervals(best), flush=True)
    report = dict(schema='bch-search-report-v1', instance=spec, target_bits=40, reused_legacy_certificate=False,
        accepted_certificates=accepted, attempts=attempts, coverage=core.summarize(best, spec['rows'], 40),
        sparse_last=last, seconds=time.monotonic()-started)
    output = directory/f'report_{len(list(directory.glob("report_*.json"))):04d}.json'
    base.write_new(output, report)
    verifier.verify_report(output)
    print('Sparse report', output, core.intervals(best), flush=True)
    return output


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--directory', type=Path, required=True)
    p.add_argument('--m', type=int, default=18)
    p.add_argument('--last', type=int, default=479)
    p.add_argument('--seconds', type=float, default=240)
    p.add_argument('--max-jobs', type=int, default=40)
    p.add_argument('--job-seconds', type=float, default=40)
    a = p.parse_args()
    run(a.directory, a.m, a.last, a.seconds, a.max_jobs, a.job_seconds)
