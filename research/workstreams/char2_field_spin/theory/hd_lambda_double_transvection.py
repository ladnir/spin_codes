"""Specialized64-XOR circuit for Hlambda+v0123*v0123^T+v0167*v0167^T.

The two quiet supports intersect in {0,1}. In the original circuit write
p_i=x_(2i)+x_(2i+1). The parity q0=p0+p1 already exists, as does r13=p1+p3
inside the a1 syndrome tree. Compute only t=p0+p3. Add r13 to outputs0,1,
q0 to outputs2,3, and t to outputs6,7. The shared outputs receive q0+t=r13.
This is one extra parity XOR and six output XORs over the57-XOR Hlambda.

The construction is the same labelled matrix screened by the frozen
multi-transvection helper. This file does not transfer any old matrix's
moments, low-arity bounds, certificate, or native performance result.
Scalar circuit counts exclude the random diagonals, setup, and routing.
"""
from collections import Counter
from hashlib import sha256
from pathlib import Path

import hd_lambda_multi_transvection as core


DIRECTORY = Path(__file__).resolve().parent
SUPPORTS = ((0, 1, 2, 3), (0, 1, 6, 7))
DEFAULT_LAMBDA = 0xbc
CONSTRUCTION = core.construction(SUPPORTS)
CORE_SHA256 = '42491167749463d4570e5685e7b69ff606746589dce03f4a5030312b215c7ab9'


def checked_dependencies():
    if sha256((DIRECTORY/'hd_lambda_multi_transvection.py').read_bytes()).hexdigest() != CORE_SHA256:
        raise ValueError('Frozen matrix helper changed')
    return {**core.checked_dependencies(), Path(__file__).name: sha256(Path(__file__).read_bytes()).hexdigest()}


def matrix(lambda_value=DEFAULT_LAMBDA, q=256):
    return core.matrix(SUPPORTS, lambda_value, q)


def _values(values, lambda_value, q):
    core.original._parameters(lambda_value, q)
    values = tuple(values)
    if len(values) != 16 or any(type(x) is not int or not 0 <= x < q for x in values):
        raise ValueError('Sixteen values in the declared field are required')
    return values


def _circuit(x, add, multiply, lambda_value):
    p = tuple(add(x[2*i], x[2*i+1]) for i in range(8))
    a0 = add(add(add(x[1], x[3]), add(x[5], x[7])), add(add(x[9], x[11]), add(x[13], x[15])))
    r13 = add(p[1], p[3])
    r57 = add(p[5], p[7])
    a1 = add(r13, r57)
    q = tuple(add(p[2*i], p[2*i+1]) for i in range(4))
    a2, a3 = add(q[1], q[3]), add(q[2], q[3])
    total = multiply(lambda_value, add(add(q[0], q[1]), a3))
    syndromes = (a0, a1, a2, a3)
    correction = [total]*16
    for i in range(1, 16):
        low = i & -i
        correction[i] = add(correction[i ^ low], syndromes[low.bit_length()-1])
    output = [add(x[i], correction[i]) for i in range(16)]
    t = add(p[0], p[3])
    output[0] = add(output[0], r13)
    output[1] = add(output[1], r13)
    output[2] = add(output[2], q[0])
    output[3] = add(output[3], q[0])
    output[6] = add(output[6], t)
    output[7] = add(output[7], t)
    return tuple(output)


def circuit(values, lambda_value=DEFAULT_LAMBDA, q=256):
    values = _values(values, lambda_value, q)
    return _circuit(values, lambda a, b: a ^ b,
                    lambda a, b: core.original._multiply(a, b, q), lambda_value)


def direct(values, lambda_value=DEFAULT_LAMBDA, q=256):
    values = _values(values, lambda_value, q)
    result = []
    for row in matrix(lambda_value, q):
        value = 0
        for coefficient, x in zip(row, values):
            value ^= core.original._multiply(coefficient, x, q)
        result.append(value)
    return tuple(result)


def primitive_counts():
    counts = Counter()
    def add(a, b):
        counts['binary_xors'] += 1
        return a ^ b
    def multiply(a, b):
        counts['fixed_products'] += 1
        return core.original._multiply(a, b, 256)
    _circuit((0,)*16, add, multiply, DEFAULT_LAMBDA)
    if (counts['binary_xors'], counts['fixed_products']) != (64, 1):
        raise ArithmeticError('Specialized circuit count changed')
    return dict(construction=CONSTRUCTION, supports=SUPPORTS,
        binary_xors=64, fixed_products=1, variable_products=0, extra_xors_over_Hlambda=7,
        generic_unfused_extra_xors=9, reused_existing_parities=('q0=p0+p1', 'r13=p1+p3'),
        new_parity='t=p0+p3', corrected_outputs=(0, 1, 2, 3, 6, 7),
        includes_diagonals=False, includes_setup=False, includes_route=False,
        native_implementation_exists=False, native_performance_claim=False,
        old_receipts_transferred=False, whole_code_certificate=False)
