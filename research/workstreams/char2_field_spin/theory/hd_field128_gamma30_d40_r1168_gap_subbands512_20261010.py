"""Four explicit smaller bands after the preserved broad partial, using the live bank."""

from hashlib import sha256
import json
from pathlib import Path

import hd_field128_gamma30_d40_r1168_gap_directed_20261010 as original

HERE = Path(__file__).resolve().parent
SOURCE = Path(__file__).resolve()
SOURCE_HASH = sha256(SOURCE.read_bytes()).hexdigest()
ORIGINAL_HASH = 'ccb8ed1f5361e66f0144cfb4296d18a34df0306cbfe961eb41f897d00b98d703'
FAILURE_HASH = 'c105375d56ff3ff09f17a96acbbbd46c3e219561e5329d7dbf87bf015966dd61'


def run(bank, original_check):
    previous = HERE/'hd_field128_gamma30_d40_r1168_gap_directed256_20261010'
    failure = previous/'band_1145_1191/summary.json'
    bank_path = previous/'bank.json'
    bank_raw = bank_path.read_bytes()
    metadata = json.loads(bank_raw)
    def check():
        original_check()
        if (sha256(SOURCE.read_bytes()).hexdigest() != SOURCE_HASH
                or sha256(original.RUNNER.read_bytes()).hexdigest() != ORIGINAL_HASH
                or sha256(failure.read_bytes()).hexdigest() != FAILURE_HASH
                or bank_path.read_bytes() != bank_raw):
            raise RuntimeError('Retry source or previous receipt changed')
    check()
    if type(bank) is not original.p.AliasBank or bank.metadata()['witness'] != metadata['witness']:
        raise ValueError('The exact prior bank must be retained unchanged')
    if metadata['calibration_active_groups'] != 1168:
        raise ValueError('Unexpected bank calibration')
    output = HERE/'hd_field128_gamma30_d40_r1168_gap_subbands512_directed256_20261010'
    output.mkdir(exist_ok=False)
    metadata['retry_provenance'] = dict(source_sha256={SOURCE.name: SOURCE_HASH},
        previous_bank=dict(path=str(bank_path), sha256=sha256(bank_raw).hexdigest()),
        excluded_prior_partial=dict(path=str(failure), sha256=FAILURE_HASH),
        bank_reused_without_reconstruction=True)
    original.write_new(output/'bank.json', metadata)
    result = original.run_bands(bank, check, output, metadata['witness_source'], ((1176, 1191), (1168, 1175), (1160, 1167), (1152, 1159)))
    check()
    pin = original.write_new(output/'retry_execution.json', dict(
        source_sha256={SOURCE.name: SOURCE_HASH, original.RUNNER.name: ORIGINAL_HASH},
        provenance=metadata['retry_provenance'], requested_ranges=result['ranges'],
        links=dict(path=str(output/'band_links.json'), sha256=sha256((output/'band_links.json').read_bytes()).hexdigest()),
        whole_code_certificate=False, union_probability_claim=False))
    print(json.dumps(dict(retry_execution=pin)), flush=True)
    return result
