"""Eight explicit FLOAT weights from exact r960 literals; no interval claim."""

import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'

from hashlib import sha256
import json
from pathlib import Path
from shutil import disk_usage

import hd_field128_gamma_midpoint_neighbors as neighbors


def main():
    here, runner = Path(__file__).resolve().parent, Path(__file__).resolve()
    source = here/'hd_field128_gamma30_d40_r960_gap_refit_20261010/proposal.json'
    expected = '80e0f727716f6c0044f81db404c7e80b121340d1f37c50d5ae3a4edfbab424b4'
    runner_hash = sha256(runner.read_bytes()).hexdigest()
    if sha256(source.read_bytes()).hexdigest() != expected:
        raise RuntimeError('Exact positive source changed')
    result = neighbors.run(source, [913, 928, 944, 952, 960, 975, 988, 999], max_seconds=90.)
    result['source_sha256'][runner.name] = runner_hash
    if (sha256(source.read_bytes()).hexdigest() != expected
            or sha256(runner.read_bytes()).hexdigest() != runner_hash):
        raise RuntimeError('Source changed during diagnostic')
    if disk_usage(here).free < 750 * 1024**2:
        raise OSError('Artifact guard: less than750MiB disk space remains')
    pin = neighbors.provider.search.write_new(here/'hd_field128_gamma30_d40_r960_gap_neighbors_float_20261010.json', result)
    print(json.dumps(dict(artifact=pin, status=result['status'])), flush=True)


if __name__ == '__main__':
    main()
