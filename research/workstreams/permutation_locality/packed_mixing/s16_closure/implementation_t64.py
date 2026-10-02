"""Bind the checked research encoder's explicit maps to the t64 proof maps.

This is map/source/log provenance, not a proof of the C++ program or of its
deterministic setup generator. The independent C++ oracles check the encoder.
The distance theorem quantifies over independent ideal uniform setup.
"""
import argparse
import contextlib
import hashlib
import io
import json
import math
from pathlib import Path
import re
import statistics
import sys

import kernel_t64

PERMUTATION = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PERMUTATION))
import gfni_t64_codegen as codegen


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def verify(directory):
    directory = Path(directory).resolve()
    wrapper, maps = kernel_t64.prepare()
    rows, columns, _, _, _, cpp_hash, source_hash = codegen.prepare()
    if (list(map(hex, rows)) != maps['expansion_rows_hex']
            or columns != maps['feedback_columns']
            or source_hash != maps['source']['sha256']):
        raise ArithmeticError('C++ generator and freshly enumerated proof maps differ')
    generated = io.StringIO()
    with contextlib.redirect_stdout(generated), contextlib.redirect_stderr(io.StringIO()):
        codegen.main()
    if (directory/'GfniT64.h').read_text() != generated.getvalue():
        raise ArithmeticError('measured generated header is not the freshly checked circuit')
    files = ('GfniT64.h','FusedR4Gfni.h','gfni_t64_probe.cpp','gfni_t64_codegen.py','gfni_t64_run.sh')
    for name in files:
        if digest(directory/name) != digest(PERMUTATION/name):
            raise ArithmeticError('local/measured source mismatch: '+name)
    build_logs = list(directory.glob('measurements/t64-build-*/hashes.txt'))
    if len(build_logs) != 1:
        raise ValueError('one immutable build provenance log required')
    build = {}
    for line in build_logs[0].read_text().splitlines():
        value, path = line.split(maxsplit=1)
        name = Path(path.strip()).name
        if name in build or not re.fullmatch('[0-9a-f]{64}', value):
            raise ValueError('duplicate build filename or malformed digest')
        build[name] = value
    for name in ('GfniT64.h','FusedR4Gfni.h','gfni_t64_probe.cpp','gfni-t64'):
        if build.get(name) != digest(directory/name):
            raise ArithmeticError('build/log hash mismatch: '+name)
    checks = list(directory.glob('measurements/t64-check-*/full-*.txt'))
    expected = {'full-14-1.txt':512,'full-14-17.txt':512,'full-20-1.txt':32768}
    if len(checks) != 3 or {p.name for p in checks} != set(expected):
        raise ValueError('both small seeds and full-size independent checks required')
    for path in checks:
        text = path.read_text()
        for token in ('checks PASS:', 'inner+FULL adjoint', 'complete output/suffix',
                      'step=64; state=16; refresh=uniform-gl16', 'map_sha256='+cpp_hash,
                      'physical_epochs='+str(expected[path.name])):
            if token not in text:
                raise ValueError('missing C++ check scope: '+token)
    confirm = list(directory.glob('measurements/t64-confirm-*'))
    if len(confirm) != 1:
        raise ValueError('one matched serial confirmation batch required')
    measurements = {'control':[], 't64':[]}
    checksums = {}
    expected_files = {f'{seed}-{rep}-{kind}.csv' for seed in (1,17)
                      for rep in (0,1) for kind in measurements}
    if {p.name for p in confirm[0].iterdir()} != expected_files:
        raise ValueError('exactly eight matched confirmation files required')
    for seed in (1,17):
        for rep in (0,1):
            for kind in measurements:
                text = (confirm[0]/f'{seed}-{rep}-{kind}.csv').read_text()
                if 'checks PASS:' not in text:
                    raise ValueError('timing must also pass correctness checks')
                if kind == 't64' and any(token not in text for token in (
                        'inner+FULL adjoint', 'complete output/suffix',
                        'step=64; state=16; refresh=uniform-gl16',
                        'map_sha256='+cpp_hash, 'physical_epochs=32768',
                        'hot_update_bytes=1048576')):
                    raise ValueError('timed t64 check has different map/geometry scope')
                data_rows = [line for line in text.splitlines() if re.match(r'^\d+,',line)]
                prefix = (['20','1',str(seed),'101'] if kind == 'control'
                          else ['20','64','16','uniform-gl16',str(seed),'101'])
                suffix = ['4','1'] if kind == 'control' else ['1']
                if len(data_rows) != 1:
                    raise ValueError('missing matched101-call process measurement')
                fields = data_rows[0].split(',')
                if (len(fields) != len(prefix)+4+len(suffix)
                        or fields[:len(prefix)] != prefix or fields[-len(suffix):] != suffix):
                    raise ValueError('timing row has different configuration')
                middle = fields[len(prefix):len(prefix)+4]
                latency, low, high = map(float,middle[:3])
                checksum = middle[3]
                if (not all(math.isfinite(value) and value > 0 for value in (latency,low,high))
                        or not low <= latency <= high or not re.fullmatch('[0-9a-f]{1,16}',checksum)):
                    raise ValueError('invalid timing quantiles or checksum')
                previous = checksums.setdefault((kind,seed),checksum)
                if previous != checksum:
                    raise ValueError('matched repetitions disagree on encoded checksum')
                measurements[kind].append(latency)
    medians = {name:statistics.median(values) for name,values in measurements.items()}
    return dict(schema='paired-t64-implementation-binding-1',
        K=1 << 20, N=1 << 21, physical_step_bits=64, state_bits=16, element_bits=128,
        map_proof_sha256=wrapper['map_sha256'], map_cpp_sha256=cpp_hash,
        map_hashes_bind_identical_explicit_matrices=True,
        source_map_sha256=source_hash, binary_sha256=digest(directory/'gfni-t64'),
        source_sha256={name:digest(directory/name) for name in files},
        linked_build_provenance=build,
        checks={p.name:digest(p) for p in checks},
        measured_process_medians_ms=measurements, median_of_process_medians_ms=medians,
        relative_latency=medians['t64']/medians['control'],
        timing_scope='Precomputed, in-place transpose; K=2^20; 128-bit elements; single CPU15. '
                     'Four serial101-call processes per candidate, alternating order; setup excluded.',
        limitation='Exact maps and recurrence are checked by independent C++ oracles. '
                   'The benchmark uses deterministic setup; the mathematical certificate uses ideal independent uniform setup.',
        whole_code_certificate=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--directory',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('fresh output required')
    result = verify(args.directory)
    with args.output.open('x') as stream:
        json.dump(result,stream,indent=2); stream.write('\n')
    print('T64 implementation binding and serial measurements checked:',
          result['median_of_process_medians_ms'],flush=True)


if __name__ == '__main__':
    main()
