"""Check scope and aggregation of a completed shared-route replay receipt.

This integrity check does not reprove the numerical leaf inequalities.
shared_relaxed_strategy.py regenerates those inequalities and the sparse
complement. Here the exact endpoints, coverage, source hash, and claimed
distance/margin are checked independently of the receipt's status labels.
"""
import argparse
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path

import shared_relaxed_strategy as proof


def audit(report, dense_raw, distance, bits=20, index=0):
    distance = Q(distance)
    if not 0 < distance < Q(1, 2) or type(bits) is not int or bits < 1:
        raise ValueError('exact relative distance and positive integer margin required')
    dense = json.loads(dense_raw)
    claim = proof.validate_record(dense, index, complete=True)
    cutoff = int(distance*(1 << 21))
    if (claim['threshold'] != cutoff or report.get('schema') != proof.SCHEMA
            or report.get('ensemble') != proof.ENSEMBLE or report.get('complete') is not True
            or type(report.get('precision')) is not int or report['precision'] < 256
            or report.get('dense_sha256') != hashlib.sha256(dense_raw).hexdigest()):
        raise ValueError('matching freshly replayed construction, cutoff, and source hash required')
    through = dense['minimum_groups']-1
    integers = dict(message_length=1 << 20, output_length=1 << 21,
                    updates=2, threshold=cutoff, minimum_distance=cutoff+1)
    if (any(type(report.get(k)) is not int or report[k] != v for k, v in integers.items())
            or report.get('sparse_occupancies') != [1, through]
            or report.get('dense_occupancies') != [through+1, 2048]):
        raise ValueError('disjoint exhaustive occupancy scope and correct strict distance required')
    leaves = claim['cover']['leaves']
    checked = report.get('dense_checked')
    if (not isinstance(checked, list) or len(checked) != len(leaves)
            or type(report.get('dense_leaves')) is not int or report['dense_leaves'] != len(leaves)
            or any(not isinstance(r, dict) or not isinstance(r.get('path'), str) for r in checked)
            or len({r['path'] for r in checked}) != len(checked)
            or {r['path'] for r in checked} != set(leaves)):
        raise ValueError('exactly one checked result for each exhaustive dense leaf required')
    dense_sum = Q(0)
    for row in checked:
        if Q(row['output_tilt']) != Q(leaves[row['path']]['witness']['parameters'][0]):
            raise ValueError('checked output tilt differs from the source witness')
        bound = proof.dyadic(row['upper'])
        if bound <= 0 or Q(row['output_tilt']) <= 0:
            raise ValueError('strictly positive dense endpoint and output tilt required')
        dense_sum += bound
    sparse_sum = proof.sum_sparse(report.get('sparse_checked'), through)
    dense_upper = proof.dyadic(report.get('dense_upper'))
    sparse_upper = proof.dyadic(report.get('sparse_upper'))
    total = proof.dyadic(report.get('total_upper'))
    if (dense_upper < dense_sum or sparse_upper < sparse_sum
            or total < dense_upper+sparse_upper or not total < Q(2)**-bits):
        raise ValueError('outward component aggregation and strict requested margin required')
    return dict(verified_scope_and_sum=True, numerical_replay_required=True,
        ensemble=proof.ENSEMBLE, minimum_distance=cutoff+1, output_length=1 << 21,
        margin_bits_at_least=bits, dense_leaves=len(leaves), sparse_through=through)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('report', type=Path)
    parser.add_argument('dense', type=Path)
    parser.add_argument('--distance', required=True)
    parser.add_argument('--bits', type=int, default=20)
    parser.add_argument('--index', type=int, default=0)
    args = parser.parse_args()
    result = audit(json.loads(args.report.read_bytes()), args.dense.read_bytes(),
                   args.distance, args.bits, args.index)
    print(json.dumps(result, indent=2))


if __name__ == '__main__': main()
