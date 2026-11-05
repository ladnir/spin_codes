"""One pinned positive-remainder calculation; no regional bank or fitting."""
from hashlib import sha256
import json
from pathlib import Path
from time import perf_counter

from hd_post_refresh_nonquiet import post_refresh_remainder


DIRECTORY = Path(__file__).resolve().parent
GATE = 'hd_post_refresh_moment_gate_float.json'
GATE_SHA = '81f8d1457d1900dacbc311ef1c55e12f446beaa8a9408e20c57f0173b2fa8b04'
HELPER_SHA = 'e5325a36c8a8e33dbc467b0847dca5dd43e0890803604004b8bc1bb413805a79'


def run():
    raw = (DIRECTORY/GATE).read_bytes()
    if (sha256(raw).hexdigest() != GATE_SHA or
            sha256((DIRECTORY/'hd_post_refresh_nonquiet.py').read_bytes()).hexdigest() != HELPER_SHA):
        raise ValueError('Reviewed source or frozen gate changed')
    data = json.loads(raw)
    if (not data['complete_requested_control'] or data['whole_code_certificate'] or
            data['construction'] != 'proposed-g16-post-refresh-coupled-inner'):
        raise ValueError('Completed scoped post-refresh FLOAT gate required')
    for name, pin in data['reviewed_source_sha256'].items():
        if sha256((DIRECTORY/name).read_bytes()).hexdigest() != pin:
            raise ValueError('Gate dependency changed: '+name)
    started = perf_counter()
    def check_time():
        if perf_counter()-started > 30:
            raise TimeoutError('Positive remainder exceeded thirty seconds')
    answer = dict(gate_artifact=GATE, gate_sha256=GATE_SHA, helper_sha256=HELPER_SHA,
                  complete_requested_component=False, whole_code_certificate=False)
    try:
        result = post_refresh_remainder(data['state_measure'], data['prefix_log2'], check_time=check_time)
        if abs(result['full_holder_log2']-data['log2_screen']) > 1e-8:
            raise ArithmeticError('Positive flags do not reproduce the frozen complete gate')
        answer.update(result, complete_requested_component=True)
    except TimeoutError as error:
        answer['incomplete_reason'] = str(error)
    answer.update(elapsed_seconds=perf_counter()-started,
                  source_sha256=sha256(Path(__file__).read_bytes()).hexdigest())
    return answer


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Do not overwrite a frozen result')
    answer = run()
    with args.output.open('x', encoding='utf-8') as handle:
        json.dump(answer, handle, indent=2, allow_nan=False)
        handle.write('\n')
    print(json.dumps(answer, allow_nan=False))
