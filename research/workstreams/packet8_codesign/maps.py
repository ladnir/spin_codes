"""Literal maps and fresh binary censuses for the byte-native candidate.

State a+256*b represents (a,b) in the AES polynomial basis. The forward
maps are A(a,b)[h]=a+h*b and C(x)=(sum x[h],sum h*x[h]), h=0,...,7.
C is not assumed to be the binary transpose of A. No distance claim is made.
"""
from __future__ import annotations

from collections import Counter
import hashlib
from itertools import combinations
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
NAMES = ('byte_native', 'byte_native_A_scaled', 'byte_native_t32', 'old_rankfix')
SCALES = (1, 3, 5, 15, 17, 51, 85, 255)
OLD_MAP = HERE.parent/'rate_quarter_bch/inner_calibration/maps/t64_s16_selected.json'


def multiply(a, b):
    result = 0
    while b:
        if b & 1:
            result ^= a
        b >>= 1
        a <<= 1
        if a & 256:
            a ^= 0x11b
    return result


def rank(columns):
    pivots = {}
    for value in map(int, columns):
        while value:
            bit = value.bit_length()-1
            if bit not in pivots:
                pivots[bit] = value
                break
            value ^= pivots[bit]
    return len(pivots)


def apply(columns, word):
    result = 0
    while word:
        low = word & -word
        result ^= int(columns[low.bit_length()-1])
        word ^= low
    return result


def images_from_rows(rows):
    result = np.zeros(1 << len(rows), dtype=np.uint64)
    for i, row in enumerate(rows):
        count = 1 << i
        result[count:2*count] = result[:count] ^ np.uint64(row)
    return result


def profiles(weights, packet_bits):
    result = np.zeros((len(weights), packet_bits+1), dtype=np.int16)
    for h in range(weights.shape[1]):
        result[np.arange(len(weights)), weights[:, h]] += 1
    return np.unique(result, axis=0, return_inverse=True)


def prepare_maps(rows, columns, *, packet_bits=8):
    """Independent A-image and C-character profiles, including every state.

    This generic interface also supports small exhaustive test fixtures.
    It does not substitute A profiles for C profiles when the maps differ.
    """
    rows, columns = tuple(map(int, rows)), tuple(map(int, columns))
    bits, width = len(rows), len(columns)
    if (not 1 <= bits <= 20 or not 1 <= packet_bits <= 8 or not width
            or width > 64 or width % packet_bits
            or any(x < 0 or x >= 1 << width for x in rows)
            or any(x < 0 or x >= 1 << bits for x in columns)
            or rank(rows) != bits or rank(columns) != bits):
        raise ValueError('full-rank binary maps and a complete packet geometry required')
    S, W, mask = 1 << bits, width//packet_bits, (1 << packet_bits)-1
    images = images_from_rows(rows)
    image_weights = np.stack([np.bitwise_count((images >> (packet_bits*h)) & mask)
                              for h in range(W)], axis=1)
    expansion_profiles, expansion_indices = profiles(image_weights, packet_bits)
    chars = np.arange(S, dtype=np.uint32)
    character_weights = np.zeros((S, W), dtype=np.uint8)
    for h in range(W):
        for c in columns[packet_bits*h:packet_bits*(h+1)]:
            character_weights[:, h] += np.bitwise_count(chars & c) & 1
    character_profiles, character_indices = profiles(character_weights, packet_bits)
    restriction_counts = {}
    for j in range(1, min(3, W)+1):
        restriction_counts[j] = dict(sorted(Counter(
            rank(columns[packet_bits*h+b] for h in support for b in range(packet_bits))
            for support in combinations(range(W), j)).items()))
    identity = dict(packet_bits=packet_bits, state_bits=bits,
        expansion_rows_hex=list(map(hex, rows)), feedback_columns=list(columns))
    record = dict(identity, physical_t=width, windows=W,
        map_sha256=hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest(),
        CA_columns=[apply(columns, row) for row in rows],
        expansion_spectrum=dict(sorted(Counter(map(int, image_weights.sum(axis=1))).items())),
        packet_restriction_rank_counts=restriction_counts,
        state_count=S, expansion_profile_count=len(expansion_profiles),
        feedback_character_profile_count=len(character_profiles), whole_code_certificate=False)
    data = dict(rows=rows, columns=columns, bits=bits, width=width, windows=W,
        packet_bits=packet_bits, images=images, expansion_profiles=expansion_profiles,
        expansion_indices=expansion_indices, character_profiles=character_profiles,
        character_indices=character_indices,
        zero_feedback_forbidden=tuple(j for j, counts in restriction_counts.items()
                                     if set(counts) == {j*packet_bits}))
    return data, record


def prepare(name='byte_native'):
    if name in ('byte_native', 'byte_native_A_scaled', 'byte_native_t32'):
        windows = 4 if name == 'byte_native_t32' else 8
        scales = SCALES if name == 'byte_native_A_scaled' else (1,)*windows
        rows = tuple(sum(multiply(scales[h], (1 << bit) if bit < 8 else multiply(h, 1 << (bit-8))) << (8*h)
                         for h in range(windows)) for bit in range(16))
        columns = tuple((1 << bit) | (multiply(h, 1 << bit) << 8)
                        for h in range(windows) for bit in range(8))
        description = ('A(a,b)[h]=lambda[h]*(a+h*b); C(x)=(sum x[h],sum h*x[h]); AES GF256'
                       if name == 'byte_native_A_scaled' else
                       'A(a,b)[h]=a+h*b; C(x)=(sum x[h],sum h*x[h]); AES GF256')
    elif name == 'old_rankfix':
        source = json.loads(OLD_MAP.read_text(encoding='utf-8'))
        old = tuple(int(row, 16) for row in source['generator_rows_hex'])
        permutation = list(range(64))
        for i in (7, 23, 39, 55):
            permutation[i], permutation[i+8] = permutation[i+8], permutation[i]
        rows = tuple(sum(((row >> old_bit) & 1) << new_bit
                         for new_bit, old_bit in enumerate(permutation)) for row in old)
        columns = tuple(sum(((row >> x) & 1) << i for i, row in enumerate(rows)) for x in range(64))
        description = 'historical selected RM2 subcode, four coordinate swaps, C=A^T'
    else:
        raise ValueError(f'unknown map {name!r}')
    data, record = prepare_maps(rows, columns)
    if ((name != 'byte_native_A_scaled' and any(record['CA_columns']))
            or record['packet_restriction_rank_counts'][1] != {8: data['windows']}):
        raise ArithmeticError('declared CA product and injective packets required')
    if (name != 'old_rankfix' and record['packet_restriction_rank_counts'][2]
            != {16: data['windows']*(data['windows']-1)//2}):
        raise ArithmeticError('every byte pair must have full feedback rank')
    record.update(name=name, description=description, original_order=True,
        field_modulus='0x11b' if name != 'old_rankfix' else None,
        alpha=list(range(data['windows'])) if name != 'old_rankfix' else None,
        coordinate_scales=list(SCALES) if name == 'byte_native_A_scaled' else 'all one',
        scales_apply_to='A only; C unchanged' if name == 'byte_native_A_scaled' else 'none',
        CA_zero=not any(record['CA_columns']),
        update='fresh independent uniform GL(2,GF256)' if name != 'old_rankfix'
               else 'fresh independent uniform GL(16,2)',
        fixed_nonzero_state_image='uniform on all 65535 nonzero binary states')
    return data, record


def sources(name='byte_native'):
    paths = [Path(__file__).resolve()]
    if name == 'old_rankfix':
        paths.append(OLD_MAP)
    return {str(path): hashlib.sha256(path.read_bytes()).hexdigest() for path in paths}
