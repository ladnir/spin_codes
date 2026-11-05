"""Exact GF256 structural probe; not a distance certificate or benchmark.

Compare A128 A128^T and A128 A64^T at state32. Enumerate every disjoint
2-by-2 minor and record bounded examples. Root runs this only on Linux.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import shutil
import sys
from time import monotonic

HERE = Path(__file__).resolve().parent
LIMIT, RESERVE = 16 << 10, 20 << 30


def mul(a, b):
    out = 0
    for _ in range(8):
        if b & 1:
            out ^= a
        b >>= 1
        a <<= 1
        if a & 256:
            a ^= 0x11b
    return out


def power(a, n):
    out = 1
    while n:
        if n & 1:
            out = mul(out, a)
        a = mul(a, a)
        n >>= 1
    return out


def lagrange(nodes):
    import numpy as np
    result = np.empty((64, 32), dtype=np.uint8)
    for x in range(64):
        for j, b in enumerate(nodes):
            numerator, denominator = 1, 1
            for a in nodes:
                if a != b:
                    numerator = mul(numerator, x ^ a)
                    denominator = mul(denominator, b ^ a)
            result[x, j] = mul(numerator, power(denominator, 254))
    return result


def product(a, b, table):
    import numpy as np
    out = np.zeros((a.shape[0], b.shape[1]), dtype=np.uint8)
    for k in range(a.shape[1]):
        out ^= table[a[:, k, None], b[None, k, :]]
    return out


def rank(a, table):
    import numpy as np
    a = np.array(a, copy=True)
    pivot = 0
    for col in range(a.shape[1]):
        candidates = np.flatnonzero(a[pivot:, col])
        if len(candidates) == 0:
            continue
        chosen = pivot+int(candidates[0])
        a[[pivot, chosen]] = a[[chosen, pivot]]
        a[pivot] = table[a[pivot], power(int(a[pivot, col]), 254)]
        factors = a[:, col].copy()
        factors[pivot] = 0
        a ^= table[factors[:, None], a[pivot][None, :]]
        pivot += 1
        if pivot == a.shape[0]:
            break
    return pivot


def minor2_report(p, table):
    import numpy as np
    columns = np.array([(a, b) for a in range(64) for b in range(a+1, 64)])
    totals, examples = {}, []
    for a in range(64):
        for b in range(a+1, 64):
            eligible = np.all((columns != a) & (columns != b), axis=1)
            c = columns[eligible]
            det = table[p[a, c[:, 0]], p[b, c[:, 1]]] ^ table[p[a, c[:, 1]], p[b, c[:, 0]]]
            for halves in (0, 1, 2):
                mask = ((c[:, 0] >= 32).astype(int)+(c[:, 1] >= 32)) == halves
                key = str(int(a >= 32)+int(b >= 32))+','+str(halves)
                item = totals.setdefault(key, dict(checked=0, singular=0))
                item['checked'] += int(np.count_nonzero(mask))
                item['singular'] += int(np.count_nonzero(mask & (det == 0)))
            if len(examples) < 8:
                for pair in c[det == 0][:8-len(examples)]:
                    cc, dd = map(int, pair)
                    examples.append(dict(rows=[a, b], columns=[cc, dd],
                        matrix=p[np.ix_([a, b], [cc, dd])].tolist()))
    return dict(counts_by_number_in_high_half=totals, examples=examples,
                total=sum(x['checked'] for x in totals.values()),
                singular=sum(x['singular'] for x in totals.values()))


def run(output):
    if sys.platform != 'linux':
        raise RuntimeError('Remote Linux only')
    import resource
    import signal
    resource.setrlimit(resource.RLIMIT_AS, (1 << 30, 1 << 30))
    resource.setrlimit(resource.RLIMIT_CPU, (40, 45))
    resource.setrlimit(resource.RLIMIT_FSIZE, (LIMIT, LIMIT))
    signal.alarm(50)
    target = (HERE/output).resolve()
    if target.parent != HERE or target.exists() or shutil.disk_usage(HERE).free-LIMIT < RESERVE:
        raise ValueError('Exclusive bounded output and reserve required')
    started = monotonic()
    source = Path(__file__).read_bytes()
    import numpy as np
    table = np.array([[mul(a, b) for b in range(256)] for a in range(256)], dtype=np.uint8)
    if not all(mul(a, power(a, 254)) == 1 for a in range(1, 256)):
        raise ArithmeticError('GF inverse check')
    emission = lagrange(tuple(range(128, 160)))
    feedback = lagrange(tuple(range(64, 96)))
    old = product(emission, emission.T, table)
    new = product(emission, feedback.T, table)
    if np.any(product(feedback.T, emission, table)) or np.any(product(new, new, table)):
        raise ArithmeticError('CA/P2 identity')
    derivative = 1
    for u in range(1, 32):
        derivative = mul(derivative, u)
    inverse_derivative = power(derivative, 254)
    for x in range(64):
        for y in range(64):
            w, value = 192 ^ x ^ y, 1
            for u in range(32):
                value = mul(value, w ^ u)
            expected = mul(mul(value, inverse_derivative), power(w, 254))
            if new[x, y] != expected:
                raise ArithmeticError('Explicit kernel identity')
    result = dict(schema='field128-feedback-rank-probe-v41', source_sha256=sha256(source).hexdigest(),
        field_polynomial='0x11b', state=32, emission_nodes=[128, 159], feedback_nodes=[64, 95],
        ranks=[rank(old, table), rank(new, table)],
        old_zero_entries=int(np.count_nonzero(old == 0)), new_zero_entries=int(np.count_nonzero(new == 0)),
        CA_zero=True, P_squared_zero=True, all4096_explicit_entries_equal=True,
        old_disjoint_minor2=minor2_report(old, table), new_disjoint_minor2=minor2_report(new, table),
        all_support_rank_claim=False, directed=False, whole_code_certificate=False,
        elapsed_seconds=monotonic()-started,
        limits=dict(memory_bytes=1 << 30, wall_seconds=50, output_bytes=LIMIT))
    data = (json.dumps(result, sort_keys=True, indent=2, allow_nan=False)+'\n').encode()
    if len(data) > LIMIT or Path(__file__).read_bytes() != source:
        raise RuntimeError('Source/output guard')
    with target.open('xb') as stream:
        stream.write(data)
    signal.alarm(0)
    return dict(path=str(target), sha256=sha256(data).hexdigest(), bytes=len(data),
                ranks=result['ranks'], old_singular=result['old_disjoint_minor2']['singular'],
                new_singular=result['new_disjoint_minor2']['singular'],
                new_examples=result['new_disjoint_minor2']['examples'], elapsed_seconds=result['elapsed_seconds'])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    print(json.dumps(run(args.output), allow_nan=False), flush=True)
