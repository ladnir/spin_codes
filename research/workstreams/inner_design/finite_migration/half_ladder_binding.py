"""Extend the retained half-rate timing binding to the K16 IMT certificate."""
import argparse
import json
from pathlib import Path

import half_binding as prior
import verify_budget


def verify(directory):
    result = prior.verify(directory)
    here,model = prior.HERE,prior.model
    proof = verify_budget.verify(16,here/'Q1_LOWER.json',here/'SPARSE_M16_v2.json',
                                 here/'DENSE_M16_coupled_v1.json')
    implementation = model.base.read(prior.IMPLEMENTATION/'IMPLEMENTATION.json')
    assert proof['instance']['inner'] == implementation['instance']['inner']
    assert proof['instance']['outer_manifest_sha256'] == implementation['instance']['outer_manifest_sha256']
    assert proof['instance']['message_bits'] == 1 << 16
    assert result['full_certificate_exponents'] == [18,20]
    first = result['cells'][0]
    assert first['message_exponent'] == 16 and not first['distance_certificate_available']
    first.update(distance_certificate_available=True,certified_margin_bits=proof['margin_bits'])
    result.update(status='VERIFIED_IMT_HALF_RATE_LADDER_TIMINGS',full_certificate_exponents=[16,18,20],
        new_serial_measurements=False,reuses_authenticated_series=True,measurement_series='half_20260916')
    result['source_sha256'].update(proof['source_sha256'])
    result['source_sha256'][Path(__file__).resolve().relative_to(model.ROOT).as_posix()] = model.base.sha(Path(__file__).resolve())
    return result


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--directory',type=Path,default=prior.HERE/'measurements/half_20260916')
    p.add_argument('--output',type=Path,required=True)
    a = p.parse_args()
    result = verify(a.directory.resolve())
    prior.model.base.write_new(a.output.resolve(),result)
    print(json.dumps({k:v for k,v in result.items() if k != 'source_sha256'},indent=2))
