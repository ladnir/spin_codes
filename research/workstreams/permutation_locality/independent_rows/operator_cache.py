"""Local memoization of outward region operators, not an independent proof.

Keys bind the generator sources, fixed maps, precision, and all parameters.
Entries preserve exact dyadic upper endpoints, never binary64 approximations.
Omit --operator-cache when independently reproducing the construction.
"""
import hashlib
import json
from pathlib import Path

from flint import arb, arb_mat


def encoded(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':')).encode()


def source_digest():
    here = Path(__file__).resolve().parent
    parent = here.parent
    paths = sorted(parent.glob('*.py')) + sorted(here.glob('*.py')) + [
        parent.parents[2]/'spin/src/kernels/generated/SelectedMaps.h',
        parent.parents[1]/'workstreams/inner_design/NO_CONSTANT_MAP.json']
    digest = hashlib.sha256()
    for path in paths:
        digest.update(str(path.relative_to(parent.parents[2])).encode())
        digest.update(b'\0')
        digest.update(path.read_bytes())
        digest.update(b'\0')
    return digest.hexdigest()


def parameters(args, tilt, penalty, sources):
    return dict(schema=2, sources=sources, precision=args.precision,
                groups=args.groups, tilt=tilt, penalty=penalty,
                full_feedback=args.full_feedback,
                window_histogram=args.window_histogram,
                joint_cancellation=args.joint_cancellation,
                column_density=getattr(args,'column_density',False),
                feedback_density=getattr(args,'feedback_density',0),
                weight_tilt=getattr(args,'weight_tilt','1'), rounds=2)


def cache_path(directory, key):
    return Path(directory)/('region-'+hashlib.sha256(encoded(key)).hexdigest()+'.json')


def save(directory, key, matrices):
    values = []
    for matrix in matrices:
        assert matrix.nrows() == matrix.ncols()
        rows = []
        for i in range(matrix.nrows()):
            row = []
            for j in range(matrix.ncols()):
                value = matrix[i,j]
                assert value.is_finite() and value >= 0
                row.append([int(x) for x in value.upper().man_exp()])
            rows.append(row)
        values.append(rows)
    body = dict(key=key, matrices=values)
    record = dict(body=body, checksum=hashlib.sha256(encoded(body)).hexdigest())
    path = cache_path(directory, key)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(encoded(record))


def load(directory, key):
    path = cache_path(directory, key)
    if not path.exists():
        return None
    record = json.loads(path.read_bytes())
    body = record['body']
    if body['key'] != key or hashlib.sha256(encoded(body)).hexdigest() != record['checksum']:
        raise ValueError('operator cache key/checksum mismatch')
    values = body['matrices']
    if len(values) != key['groups']+1:
        raise ValueError('operator cache degree mismatch')
    result = []
    for rows in values:
        if len(rows) != 11 or any(len(row) != 11 for row in rows):
            raise ValueError('operator cache dimension mismatch')
        exact = []
        for row in rows:
            converted = []
            for mantissa, exponent in row:
                if not isinstance(mantissa, int) or not isinstance(exponent, int) or mantissa < 0:
                    raise ValueError('invalid cached dyadic')
                value = arb(mantissa)*arb(2)**exponent
                if not value.is_exact() or not value.is_finite():
                    raise ValueError('cache precision does not preserve exact upper endpoint')
                converted.append(value)
            exact.append(converted)
        result.append(arb_mat(exact))
    return result
