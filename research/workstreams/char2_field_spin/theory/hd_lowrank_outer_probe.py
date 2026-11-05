"""Exact local checks for stronger low-rank affine-coset outer mixers.

H_d=I+V_d V_d^T, V[i,j]=u_i^j, j=0..d-1, on the existing16 labels.
This probes algebraic prerequisites only: no distance or timing claim.
The candidate changes the code and is not the frozen G16 construction.
"""
from collections import Counter
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
from time import perf_counter

from gf8_packet_maps import MUL, rank, matmul, xor_sum
from hd_two_sum_gf16_sparse import LABELS, matrix as reference_matrix


def vandermonde(degree):
    if type(degree) is not int or not 2 <= degree <= 4:
        raise ValueError('Only the bounded rank2/3/4 designs are in scope')
    rows = []
    for u in LABELS:
        row = [1]
        for _ in range(1, degree):
            row.append(MUL[row[-1]][u])
        rows.append(tuple(row))
    return tuple(rows)


def matrix(degree):
    v = vandermonde(degree)
    return tuple(tuple(int(i == j)^xor_sum(MUL[a][b] for a, b in zip(v[i], v[j]))
                       for j in range(16)) for i in range(16))


def principal_histogram(degree, check_time):
    """All principal nullities through a degree<=4 determinant lemma.

    rank(I_a+V_A V_A^T)=a-d+rank(I_d+V_A^T V_A).
    Gray-code updates change one outer product, without reusing setup.
    """
    v = vandermonde(degree)
    products = [tuple(MUL[a][b] for a in row for b in row) for row in v]
    small = [int(i == j) for i in range(degree) for j in range(degree)]
    histogram, examples = Counter(), {}
    previous, active = 0, 0
    for index in range(1 << 16):
        if index % 1024 == 0:
            check_time()
        gray = index^(index >> 1)
        if index:
            bit = (gray^previous).bit_length()-1
            active += 1 if gray >> bit & 1 else -1
            small = [a^b for a, b in zip(small, products[bit])]
        nullity = degree-rank([small[j:j+degree] for j in range(0, degree*degree, degree)])
        histogram[active, nullity] += 1
        if nullity:
            examples.setdefault(nullity, gray)
        previous = gray
    if sum(histogram.values()) != 1 << 16:
        raise ArithmeticError('Missing principal subsets')
    return [dict(arity=a, nullity=n, count=c) for (a, n), c in sorted(histogram.items())], examples


def probe(degree, check_time):
    h, v = matrix(degree), vandermonde(degree)
    vt = tuple(zip(*v))
    gram = matmul(vt, v)
    identity = tuple(tuple(int(i == j) for j in range(16)) for i in range(16))
    if any(x for row in gram for x in row) or matmul(h, h) != identity or rank(h) != 16:
        raise ArithmeticError('Moment cancellation or invertibility failed')
    pairs = list(combinations(range(16), 2))
    zero_minors = sum((MUL[h[i][a]][h[j][b]]^MUL[h[i][b]][h[j][a]]) == 0
                     for i, j in pairs for a, b in pairs)
    pair_histograms = Counter()
    for a, b in pairs:
        check_time()
        distribution = Counter(sum(bool(h[i][a]^MUL[h[i][b]][ratio]) for i in range(16))
                               for ratio in range(1, 256))
        pair_histograms[tuple(sorted(distribution.items()))] += 1
    principal, examples = principal_histogram(degree, check_time)
    return dict(degree=degree, matrix=h, gram_zero=True, involution=True, full_rank=16,
        column_weights=[sum(bool(h[i][j]) for i in range(16)) for j in range(16)],
        two_by_two_zero_minors=zero_minors, total_two_by_two_minors=len(pairs)**2,
        pair_weight_laws=[dict(pair_count=count, counts_by_weight=dict(hist))
                          for hist, count in sorted(pair_histograms.items())],
        principal_nullities=principal, singular_principal_examples=examples,
        all_principal_subsets_invertible=not examples,
        same_as_frozen_rank2=(h == reference_matrix()),
        whole_code_certificate=False, native_implementation=False, measured_performance=False)


def run(max_seconds=30, progress=None):
    if type(max_seconds) not in (int, float) or not 0 < max_seconds <= 30:
        raise ValueError('An explicit at-most30s algebra budget is required')
    started = perf_counter()
    def check_time():
        if perf_counter()-started >= max_seconds:
            raise TimeoutError('Low-rank local algebra budget exhausted')
    result = dict(labels=list(LABELS), field_order=256, field_polynomial='0x11b',
        family='H_d=I+V_d V_d^T; V[i,j]=u_i^j', requested_degrees=[2, 3, 4],
        complete=False, compute_budget_seconds=max_seconds, entries=[])
    try:
        for degree in (2, 3, 4):
            row = probe(degree, check_time)
            result['entries'].append(row)
            if progress:
                progress(dict(degree=degree, singular_principal_examples=row['singular_principal_examples'],
                    minimum_column_weight=min(row['column_weights']),
                    zero_2x2_minors=row['two_by_two_zero_minors'], elapsed_seconds=perf_counter()-started))
        result.update(status='completed', complete=True)
    except TimeoutError as error:
        result.update(status='compute-budget-exceeded', reason=str(error))
    result['elapsed_seconds'] = perf_counter()-started
    return result


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--max-seconds', type=float, default=30)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Do not overwrite a frozen artifact')
    result = run(args.max_seconds, progress=lambda row: print(json.dumps(row), flush=True))
    result['source_sha256'] = sha256(Path(__file__).read_bytes()).hexdigest()
    with args.output.open('x', encoding='utf-8') as target:
        json.dump(result, target, indent=2)
        target.write('\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'entries'}))
