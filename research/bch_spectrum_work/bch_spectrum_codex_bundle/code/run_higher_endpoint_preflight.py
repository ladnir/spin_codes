"""Validate and preflight the weight-40/42 kernels; never accept a shell cap.

The two fixed-size preflights run sequentially under the existing experiment
mutex. Sources, executable, and Python positive controls are hashed before CNG
sampling. Independent replay and full Gaussian trace checks follow each run.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import random
import subprocess
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

sys.dont_write_bytecode = True
from run_endpoint_certificate import ExperimentLock
from verify_higher_shell_endpoint_reduction import Field, endpoint_roots, witness_words

ROOT = Path(__file__).resolve().parents[1]
EXE = ROOT / 'generated/test_bch_higher_endpoints.exe'
PROTOCOL = 'fixed-bch256-higher-endpoint-v1'


def sha(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def write_new(path, payload):
    with Path(path).open('x', encoding='utf-8') as stream:
        json.dump(payload, stream, indent=2)
        stream.write('\n')


def full_solve(field, nodes, sigma, weight):
    r, m = (weight-38)//2, 18+(weight-38)//2
    matrix = [[field.power(x, j) for j in range(1, m+1)] +
              [field.power(x, 146+j) for j in range(1, r+1)] for x in nodes]
    rhs = [1 ^ field.mul[sigma][field.power(x, 146)] for x in nodes]
    solution = field.solve(matrix, rhs)
    if solution is None:
        return dict(valid=False, agreements=0, g=[0]*(m+1), h=[0]*(r+1))
    g, h = [1]+solution[:m], [sigma]+solution[m:]
    return dict(valid=True, agreements=len(endpoint_roots(field, g, h, 146)), g=g, h=h)


def prepare_vectors(directory):
    field = Field(256, 0x14d)
    powers = [field.power(2, j) for j in range(255)]
    rng = random.Random(2026090601)
    records = []
    sources = set()
    coverage = Counter()

    def add(weight, nodes, sigma, category, expected=None):
        result = full_solve(field, nodes, sigma, weight)
        if expected is not None:
            g, h = expected
            assert result == dict(valid=True, agreements=weight-1, g=g, h=h)
        records.append(dict(weight=weight, sigma=sigma, nodes=nodes, **result))
        coverage[f'w{weight}_{category}'] += 1
        coverage[f'w{weight}_singular' if not result['valid'] else
                 f'w{weight}_endpoints' if result['agreements']==weight-1 else
                 f'w{weight}_nonendpoints'] += 1

    for weight in (40, 42):
        n = weight-20
        for family in ('p', 'q'):
            path, words = witness_words(weight, family)
            sources.add(path)
            for word in words:
                support = [powers[j] for j in range(255) if (word>>j)&1]
                if (word>>255)&1:
                    support.append(0)
                assert len(support)==weight
                anchor = support[0]
                punctured = [x^anchor for x in support if x!=anchor]
                locator = field.product(punctured, locator=True)
                sigma = locator[37]
                if sigma:
                    scale = field.power(field.inv[sigma], pow(37, -1, 255))
                    punctured = [field.mul[scale][x] for x in punctured]
                    locator = field.product(punctured, locator=True)
                sigma = locator[37]
                assert sigma in (0, 1)
                g, h = locator[::2], locator[37::2]
                roots = [field.power(x, 253) for x in punctured]
                for _ in range(64):
                    add(weight, rng.sample(roots, n), sigma, f'positive_sigma{sigma}', (g,h))
        # Fixed-seed diagnostic vectors, with no role in statistical acceptance.
        for _ in range(1024):
            add(weight, rng.sample(range(1, 256), n), rng.randrange(32), 'arbitrary')
        assert coverage[f'w{weight}_singular'] > 0
        assert coverage[f'w{weight}_nonendpoints'] > 0
        assert coverage[f'w{weight}_positive_sigma0'] > 0
        assert coverage[f'w{weight}_positive_sigma1'] > 0
    with (directory/'vectors.txt').open('x', encoding='ascii') as stream:
        for row in records:
            fields = [row['weight'],row['sigma'],int(row['valid']),row['agreements'],
                      *row['nodes'],*row['g'],*row['h']]
            stream.write(' '.join(map(str, fields))+'\n')
    summary = dict(classification='deterministic structural fixtures, not statistical evidence',
                   vectors=len(records), coverage=dict(coverage),
                   vectors_sha256=sha(directory/'vectors.txt'),
                   witness_sha256={str(p.relative_to(ROOT)):sha(p) for p in sorted(sources)})
    write_new(directory/'vectors.json', summary)
    print(json.dumps(summary), flush=True)


def verify_trace(directory, weight, samples):
    field = Field(256, 0x14d)
    checked, singular, endpoints = 0, 0, 0
    previous = -1
    seen = set()
    with (directory/'random.bin').open('rb') as tape, \
         (directory/'sample.json.audit.jsonl').open() as stream:
        for line in stream:
            row = json.loads(line)
            index = row['sample']
            assert previous < index < samples
            previous = index
            seen.add(index)
            begin, end = row['byte_offset'], row['next_byte_offset']
            assert 0 <= begin < end
            tape.seek(begin)
            raw = tape.read(end-begin)
            assert len(raw) == end-begin and (raw[0]&31)==row['sigma']
            nodes = []
            for offset, x in enumerate(raw[1:], start=1):
                if x and x not in nodes:
                    nodes.append(x)
                if len(nodes)==weight-20:
                    assert offset == len(raw)-1
                    break
            assert len(nodes)==weight-20 and nodes==row['nodes']
            result = full_solve(field, nodes, row['sigma'], weight)
            assert all(row[key]==value for key,value in result.items())
            checked += 1
            singular += not result['valid']
            endpoints += result['valid'] and result['agreements']==weight-1
    assert set(range(min(256,samples))) <= seen and samples-1 in seen
    return dict(records=checked, singular_records=singular, endpoint_records=endpoints,
                tape_decode_and_full_gaussian_checks_passed=True)


def frozen_inputs():
    names = ['code/test_bch_higher_endpoints.cpp', 'generated/test_bch_higher_endpoints.exe',
             'code/run_higher_endpoint_preflight.py', 'code/verify_higher_shell_endpoint_reduction.py',
             'code/run_endpoint_certificate.py', 'code/analyze_endpoint_sampling.py',
             'code/bch_quotient.py', 'code/affine_wambach.py']
    return {name:sha(ROOT/name) for name in names}


def verify(directory):
    manifest = json.loads((directory/'manifest.json').read_text())
    assert manifest['protocol']==PROTOCOL and manifest['diagnostic_only'] is True
    assert manifest['input_sha256']==frozen_inputs()
    vector_report = json.loads((directory/'vectors.json').read_text())
    assert sha(directory/'vectors.txt')==vector_report['vectors_sha256']==manifest['vectors_sha256']
    assert sha(directory/'vectors.json')==manifest['vectors_report_sha256']
    assert sha(directory/'kernel_check.txt')==manifest['kernel_check_sha256']
    for name, expected in vector_report['witness_sha256'].items():
        assert sha(ROOT/name)==expected
    results = {}
    for weight in (40,42):
        sub = directory/f'w{weight}'
        primary = json.loads((sub/'sample.json').read_text())
        replay = json.loads((sub/'replay.json').read_text())
        assert primary['mode']=='BCryptGenRandom' and replay['mode']=='independent_replay'
        keys = ('weight','samples','endpoint_hits','singular_queries','bytes_consumed','bytes_loaded',
                'zero_bytes_rejected','duplicate_bytes_rejected','agreement_histogram')
        assert all(primary[k]==replay[k] for k in keys)
        assert primary['weight']==weight and primary['samples']==manifest['samples_per_weight']
        n, histogram = primary['samples'],primary['agreement_histogram']
        assert len(histogram)==weight and all(isinstance(v,int) and v>=0 for v in histogram)
        assert sum(histogram)+primary['singular_queries']==n
        assert histogram[weight-1]==primary['endpoint_hits'] and sum(histogram[:weight-20])==0
        assert primary['bytes_consumed']==(weight-19)*n+primary['zero_bytes_rejected']+primary['duplicate_bytes_rejected']
        size = (sub/'random.bin').stat().st_size
        assert size==primary['bytes_loaded'] and 0<=size-primary['bytes_consumed']<1<<20
        assert sha(sub/'sample.json.audit.jsonl')==sha(sub/'replay.json.audit.jsonl')
        results[str(weight)] = dict(samples=n, endpoint_hits=primary['endpoint_hits'],
            singular_queries=primary['singular_queries'],
            full_independent_algorithm_replay_matches=True,
            python_trace=verify_trace(sub,weight,n), primary_seconds=primary['seconds'],
            replay_seconds=replay['seconds'], random_tape_bytes=size,
            sha256={name:sha(sub/name) for name in ('random.bin','sample.json','replay.json',
                                                   'sample.json.audit.jsonl','replay.json.audit.jsonl')})
    return dict(status='validated_preflight_only', shell_caps_established=False,
                protocol=PROTOCOL, manifest_sha256=sha(directory/'manifest.json'),
                positive_and_arbitrary_vectors=vector_report['vectors'],
                results=results,
                scope='Implementation validation only; no statistical shell acceptance, no deterministic spectrum bound, and no full SPIN claim.')


def execute(directory, samples):
    assert 1<=samples<=1_000_000, 'This driver is preflight-only, capped at one million trials per weight'
    with ExperimentLock():
        directory.mkdir(parents=True, exist_ok=False)
        prepare_vectors(directory)
        checked = subprocess.run([str(EXE),'check',str(directory/'vectors.txt')],
                                 check=True,capture_output=True,text=True)
        with (directory/'kernel_check.txt').open('x') as stream:
            stream.write(checked.stdout)
        print(checked.stdout,flush=True)
        manifest = dict(protocol=PROTOCOL, created_utc=datetime.now(timezone.utc).isoformat(),
            diagnostic_only=True, samples_per_weight=samples, weights=[40,42],
            stopping_rule='Exactly the predeclared sample count; singular queries count as non-hits, no resampling or retry',
            randomness='Windows CNG; entire 1 MiB chunks retained; conditional IID model only',
            sigma_sampling='Read one byte and take its low five bits, giving S={0,...,31}',
            subset_sampling='Reject zero and within-subset duplicates until w-20 nonzero distinct bytes; do not reject singular systems',
            primary='Newton divided differences, shifted polynomial remainders, closed 1x1/2x2 solve, AVX2 term-table root count',
            replay='Separate incremental interpolation for each remainder, Gaussian 1x1/2x2 solve, eight-lane scalar Horner root count',
            trace_checker='Full 20x20/22x22 Python Gaussian solve, all-field root count, and independent tape decoding',
            input_sha256=frozen_inputs(), vectors_sha256=sha(directory/'vectors.txt'),
            vectors_report_sha256=sha(directory/'vectors.json'),
            kernel_check_sha256=sha(directory/'kernel_check.txt'))
        write_new(directory/'manifest.json',manifest)
        print('FROZEN',directory/'manifest.json',sha(directory/'manifest.json'),flush=True)
        for weight in (40,42):
            sub = directory/f'w{weight}'
            sub.mkdir()
            for mode,name in (('run','sample.json'),('replay','replay.json')):
                subprocess.run([str(EXE),mode,str(weight),str(samples),str(sub/'random.bin'),
                                str(sub/name),PROTOCOL],check=True)
        report = verify(directory)
        write_new(directory/'preflight.json',report)
        print(json.dumps(report,indent=2),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('operation',choices=('execute','verify'))
    parser.add_argument('directory',type=Path)
    parser.add_argument('--samples',type=int,default=65536)
    args=parser.parse_args()
    if args.operation=='execute':
        execute(args.directory.resolve(),args.samples)
    else:
        print(json.dumps(verify(args.directory.resolve()),indent=2))
