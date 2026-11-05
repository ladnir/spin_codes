"""Seven explicit FLOAT weights from exact r1168 literals; no interval claim."""

import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'

from hashlib import sha256
import json
from pathlib import Path
from shutil import disk_usage

import hd_field128_gamma_midpoint_neighbors as neighbors


def main():
    here, runner = Path(__file__).resolve().parent, Path(__file__).resolve()
    source = here/'hd_field128_gamma30_d40_r1168_gap_calibration_20261010/proposal.json'
    expected = 'fd0ac0024464d8d2f852dde4bb3bdc5566eac2bcbd1ab5f0c5577876da9bff81'
    runner_hash = sha256(runner.read_bytes()).hexdigest()
    if sha256(source.read_bytes()).hexdigest() != expected:
        raise RuntimeError('Exact positive source changed')
    result = neighbors.run(source, [1145, 1152, 1160, 1168, 1176, 1184, 1191], max_seconds=90.)
    result['source_sha256'][runner.name] = runner_hash
    if (sha256(source.read_bytes()).hexdigest() != expected
            or sha256(runner.read_bytes()).hexdigest() != runner_hash):
        raise RuntimeError('Source changed during diagnostic')
    if disk_usage(here).free < 750 * 1024**2:
        raise OSError('Artifact guard: less than750MiB disk space remains')
    pin = neighbors.provider.search.write_new(here/'hd_field128_gamma30_d40_r1168_gap_neighbors_float_20261010.json', result)
    print(json.dumps(dict(artifact=pin, status=result['status'])), flush=True)


if __name__ == '__main__':
    main()
