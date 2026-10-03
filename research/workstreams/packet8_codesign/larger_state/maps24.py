"""Literal GF(256)^3 maps and bounded-memory, exact integer profile census."""
from collections import Counter
import hashlib
from itertools import combinations, combinations_with_replacement
import json
from pathlib import Path

import numpy as np


def multiply(a, b, *, bits=8, modulus=0x11b):
    value = 0
    while b:
        if b & 1:
            value ^= a
        b >>= 1
        a <<= 1
        if a & (1 << bits):
            a ^= modulus
    return value


def apply(columns, word):
    value = 0
    while word:
        low = word & -word
        value ^= int(columns[low.bit_length()-1])
        word ^= low
    return value


def rank(columns):
    pivots = {}
    for value in map(int, columns):
        while value:
            top = value.bit_length()-1
            if top not in pivots:
                pivots[top] = value
                break
            value ^= pivots[top]
    return len(pivots)


def make_maps(*, packet_bits=8, modulus=0x11b, windows=8):
    b, Q = packet_bits, 1 << packet_bits
    if windows > Q or packet_bits > 8 or windows*packet_bits > 64:
        raise ValueError('small field and at most64 output bits required')
    mul = lambda a, c: multiply(a, c, bits=b, modulus=modulus)
    coefficients = tuple((1, h, mul(h, h)) for h in range(windows))
    tables = np.array([[[mul(coefficient, x) for x in range(Q)]
        for coefficient in row] for row in coefficients], dtype=np.uint8)
    # Transposes use the literal polynomial-basis binary dot product.
    transpose = np.zeros_like(tables)
    for h in range(windows):
        for k in range(3):
            for x in range(Q):
                transpose[h, k, x] = sum(((int(tables[h, k, 1 << i]) & x).bit_count() & 1) << i
                                           for i in range(b))
    rows = tuple(sum(int(tables[h, i//b, 1 << (i % b)]) << (b*h)
                     for h in range(windows)) for i in range(3*b))
    columns = tuple(sum(int(tables[h, k, 1 << i]) << (b*k) for k in range(3))
                    for h in range(windows) for i in range(b))
    restrictions = {j: dict(sorted(Counter(rank(columns[b*h+i] for h in support for i in range(b))
        for support in combinations(range(windows), j)).items()))
        for j in range(1, min(windows, 4)+1)}
    identity = dict(packet_bits=b, windows=windows, state_bits=3*b, modulus=hex(modulus),
        expansion_rows_hex=list(map(hex, rows)), feedback_columns=list(columns))
    record = dict(identity, map_sha256=hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest(),
        expansion_rank=rank(rows), feedback_rank=rank(columns),
        CA_columns=[apply(columns, row) for row in rows],
        packet_restriction_rank_counts=restrictions,
        order='y=X+A*a; a_next=M*a+C*X',
        update=f'fresh independent uniform GL(3,GF({Q}))',
        zero_initial_state=True, final_flush=False, whole_code_certificate=False)
    if rank(rows) != 3*b or rank(columns) != 3*b:
        raise ArithmeticError('full state rank required')
    if b == 8 and (any(record['CA_columns']) or restrictions[3] != {24: 56}):
        raise ArithmeticError('literal candidate must have CA=0 and full-rank triples')
    # Four bits suffice to encode one byte Hamming weight, including8.
    choices = list(combinations_with_replacement(range(b+1), windows))
    keys = np.array([sum(int(v) << (4*h) for h, v in enumerate(choice))
                     for choice in choices], dtype=np.uint32)
    order = np.argsort(keys)
    choices = np.array(choices, dtype=np.uint8)[order]
    profiles = np.array([np.bincount(choice, minlength=b+1) for choice in choices], dtype=np.uint8)
    return dict(bits=3*b, packet_bits=b, windows=windows, rows=rows, columns=columns,
        tables=tables, transpose_tables=transpose, profile_keys=keys[order],
        profiles=profiles, profile_weights=choices.sum(axis=1), record=record)


def profile_indices(data, states, *, character=False):
    b, W = data['packet_bits'], data['windows']
    states = np.asarray(states, dtype=np.uint32)
    mask = (1 << b)-1
    coordinates = [(states >> (b*k)) & mask for k in range(3)]
    tables = data['transpose_tables'] if character else data['tables']
    weights = np.empty((len(states), W), dtype=np.uint8)
    for h in range(W):
        word = tables[h, 0, coordinates[0]] ^ tables[h, 1, coordinates[1]] ^ tables[h, 2, coordinates[2]]
        weights[:, h] = np.bitwise_count(word)
    weights.sort(axis=1)
    keys = np.zeros(len(states), dtype=np.uint32)
    for h in range(W):
        keys |= weights[:, h].astype(np.uint32) << (4*h)
    indices = np.searchsorted(data['profile_keys'], keys)
    if np.any(data['profile_keys'][indices] != keys):
        raise ArithmeticError('incomplete profile catalog')
    return indices.astype(np.uint16)


def census(data, *, chunk_bits=16, include_character=True):
    S, count = 1 << data['bits'], len(data['profiles'])
    expansion = np.empty(S, dtype=np.uint16)
    character = np.empty(S, dtype=np.uint16) if include_character else None
    histogram = np.zeros(count, dtype=np.int64)
    chunk = 1 << chunk_bits
    for start in range(0, S, chunk):
        states = np.arange(start, min(S, start+chunk), dtype=np.uint32)
        index = profile_indices(data, states)
        expansion[start:start+len(states)] = index
        histogram += np.bincount(index, minlength=count)
        if include_character:
            character[start:start+len(states)] = profile_indices(data, states, character=True)
    spectrum = np.bincount(data['profile_weights'].astype(int), weights=histogram,
                           minlength=data['packet_bits']*data['windows']+1).astype(np.int64)
    if histogram.sum() != S or spectrum[0] != 1:
        raise ArithmeticError('census must include exactly one zero expansion')
    return dict(expansion_indices=expansion, character_indices=character,
                histogram=histogram, spectrum=spectrum)


def sources():
    path = Path(__file__).resolve()
    return {str(path): hashlib.sha256(path.read_bytes()).hexdigest()}
