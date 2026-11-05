"""Retry the three exact ranges after the recorded storage stop, using the live bank."""

from hashlib import sha256
import json
from pathlib import Path

import hd_field128_gamma30_d40_gap2336_2944_directed_20261010 as original

HERE = Path(__file__).resolve().parent
SOURCE = Path(__file__).resolve()
SOURCE_HASH = sha256(SOURCE.read_bytes()).hexdigest()
ORIGINAL_HASH = '8cd9ed5def4d3b8803fb1482a36f95073df23e7d5784ba8396003678da746228'
FAILURE_HASH = '5d7ba28272946fdc8664c96b6714c709348bde8b2d8b1767dfa672de1e238dee'


def run(bank, original_check):
    previous = HERE/'hd_field128_gamma30_d40_r2336_gap_directed256_20261010'
    failure = previous/'band_2304_2367/failure.json'
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
    if metadata['calibration_active_groups'] != 2336:
        raise ValueError('Unexpected bank calibration')
    output = HERE/'hd_field128_gamma30_d40_r2336_gap_retry02_directed256_20261010'
    output.mkdir(exist_ok=False)
    metadata['retry_provenance'] = dict(source_sha256={SOURCE.name: SOURCE_HASH},
        previous_bank=dict(path=str(bank_path), sha256=sha256(bank_raw).hexdigest()),
        excluded_resource_stop=dict(path=str(failure), sha256=FAILURE_HASH),
        bank_reused_without_reconstruction=True)
    original.write_new(output/'bank.json', metadata)
    result = original.run_bands(bank, check, output, metadata['witness_source'], original.CONFIG[2336][2])
    check()
    pin = original.write_new(output/'retry_execution.json', dict(
        source_sha256={SOURCE.name: SOURCE_HASH, original.RUNNER.name: ORIGINAL_HASH},
        provenance=metadata['retry_provenance'], requested_ranges=result['ranges'],
        links=dict(path=str(output/'band_links.json'), sha256=sha256((output/'band_links.json').read_bytes()).hexdigest()),
        whole_code_certificate=False, union_probability_claim=False))
    print(json.dumps(dict(retry_execution=pin)), flush=True)
    return result
