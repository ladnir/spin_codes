"""Bounded exact comparison of all 65 V57 local return caps at V52 witness0.

Exact rational inequalities are checked before logarithms are printed.
This is a local diagnostic, not a regional or whole-code certificate.
"""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from math import log2
from pathlib import Path
import resource
import shutil
from time import monotonic

import hd_field128_return_tail_moments_v57 as candidate

HERE = Path(__file__).resolve().parent
SOURCE = 'hd_field128_reversed_dual_marker_v52_20261010.json'
SOURCE_SHA = '1ec14fce6aede5f7f45d301fdeae5cb7057670a09f8e244bad52e9bd1b69b3e5'
PROVIDER_SHA = '8c1ed22d4ccd3a7033bbd384f23e46291785c0c1bafc6f332a277b08aec61bef'
LIMIT, RESERVE = 48 << 10, 20 << 30


def log_fraction(value):
    return None if not value else log2(value.numerator)-log2(value.denominator)


def rational_digest(value):
    digest = sha256()
    for integer in (value.numerator, value.denominator):
        raw = integer.to_bytes(max(1, (integer.bit_length()+7)//8), 'big')
        digest.update(len(raw).to_bytes(8, 'big'))
        digest.update(raw)
    return digest.hexdigest()


def run(output_name):
    if HERE != Path('/tmp/spin-field128-proof-52GBmlZD/bundle/theory'):
        raise RuntimeError('This bounded numerical job runs only in remote scratch')
    if Path(output_name).name != output_name or not output_name.endswith('.json'):
        raise ValueError('Expected a result basename')
    output = HERE/output_name
    if output.exists() or shutil.disk_usage(HERE).free-LIMIT < RESERVE:
        raise ValueError('Exclusive output and planned reserve required')
    resource.setrlimit(resource.RLIMIT_AS, (1 << 30, 1 << 30))
    resource.setrlimit(resource.RLIMIT_CPU, (150, 155))
    resource.setrlimit(resource.RLIMIT_FSIZE, (LIMIT, LIMIT))
    source = HERE/SOURCE
    if source.stat().st_size > 96 << 10:
        raise ValueError('Oversized input')
    raw = source.read_bytes()
    if sha256(raw).hexdigest() != SOURCE_SHA:
        raise ValueError('V52 result changed')
    pins = candidate.source_pins()
    if pins[candidate.NAME] != PROVIDER_SHA:
        raise ValueError('V57 provider changed')
    own = Path(__file__).read_bytes()
    spec = json.loads(raw)['witness_specs'][0]
    z = F(spec['z'])
    start, rows = monotonic(), []
    for j in range(65):
        if monotonic()-start > 145:
            raise TimeoutError('Local-cap wall budget exceeded')
        values = candidate.cap_components(j, z)
        old = candidate.provider.exact_moments(j, z)
        new = candidate.exact_moments(j, z)
        mass = candidate.factorial_bounds(j)[0]
        if not z**64*mass <= new[2] <= old[2]:
            raise ArithmeticError('Exact return enclosure failed')
        if (old[0], old[1], old[3]) != (new[0], new[1], new[3]):
            raise ArithmeticError('A non-return packet entry changed')
        if values['selected'] != new[2]:
            raise ArithmeticError('Local interfaces disagree')
        rows.append(dict(occupancy=j, maximum_zeros=values['maximum_zeros'],
            improved=new[2] < old[2],
            log2_gain=log_fraction(old[2]/new[2]) if new[2] else 0.,
            cap_logs={key:log_fraction(value) for key, value in values.items()
                      if key != 'maximum_zeros'},
            selected_exact_sha256=rational_digest(new[2])))
    candidate.check_source_pins(pins)
    if source.read_bytes() != raw or Path(__file__).read_bytes() != own:
        raise RuntimeError('Inputs/source changed during comparison')
    pins[Path(__file__).name] = sha256(own).hexdigest()
    result = dict(schema='field128-return-tail-local-v58', source_sha256=pins,
        input=dict(path=SOURCE, sha256=SOURCE_SHA), witness_index=0,
        z=spec['z'], alpha=spec['alpha'], occupancies=rows,
        improved_occupancies=sum(row['improved'] for row in rows),
        maximum_log2_gain=max(row['log2_gain'] for row in rows),
        exact_rational_comparisons=True, log_values_are_float=True,
        digest_format='SHA256 of numerator then denominator, each u64be length plus unsigned big-endian bytes',
        construction_changed=False, principal_setup_unchanged=True,
        original_input_normalization=True, whole_code_certificate=False,
        regional_replay=False, elapsed_seconds=monotonic()-start)
    data = (json.dumps(result, sort_keys=True, separators=(',', ':'),
                       allow_nan=False)+'\n').encode()
    if len(data) > LIMIT or shutil.disk_usage(HERE).free-len(data) < RESERVE:
        raise RuntimeError('Output budget exceeded')
    with output.open('xb') as stream:
        stream.write(data)
    return dict(path=output_name, bytes=len(data), sha256=sha256(data).hexdigest(),
        improved_occupancies=result['improved_occupancies'],
        maximum_log2_gain=result['maximum_log2_gain'],
        gains=[(r['occupancy'], r['log2_gain']) for r in rows if r['improved']],
        elapsed_seconds=result['elapsed_seconds'])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True)
    print(json.dumps(run(parser.parse_args().output), allow_nan=False), flush=True)
