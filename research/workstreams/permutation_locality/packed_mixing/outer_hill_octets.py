"""Exact small canonical-octet count; no inner or whole-code claim.

For the even production BCH code, minimum distance38 implies that a word
supported on five octets has weight38 or40. It is therefore the union of
those octets with zero or two positions deleted. Syndrome lookups enumerate
these possibilities exactly, without a rank computation for every block set.
"""
import argparse
from collections import Counter
import hashlib
from itertools import combinations
import json
from math import comb
from pathlib import Path
import re


P = 0x1C3B42FE115AA90747020D79AC738ADB
MODULUS = 0x14D


def gf_mul(a, b):
    result = 0
    while b:
        if b & 1:
            result ^= a
        b >>= 1
        a <<= 1
        if a & 256:
            a ^= MODULUS
    return result


def remainder(word, polynomial):
    if polynomial <= 0:
        raise ValueError('nonzero binary polynomial required')
    degree = polynomial.bit_length()
    while word.bit_length() >= degree:
        word ^= polynomial << (word.bit_length()-degree)
    return word


def syndrome_columns(rows, length):
    """Return a quotient-map column table whose kernel is exactly span(rows)."""
    if type(length) is not int or length < 1:
        raise ValueError('positive binary code length required')
    basis = {}
    for word in rows:
        if type(word) is not int or not 0 < word < 1 << length:
            raise ValueError('nonzero binary generator rows required')
        while word:
            bit = word.bit_length()-1
            if bit not in basis:
                basis[bit] = word
                break
            word ^= basis[bit]
    ordered = sorted(basis.items(), reverse=True)
    columns = []
    for i in range(length):
        word = 1 << i
        for bit, row in ordered:
            if word >> bit & 1:
                word ^= row
        columns.append(word)
    return columns, len(basis)


def verify_distance38(rows):
    """Check the even extension of a primitive BCH code with36 roots.

    The BCH bound gives punctured distance at least37. Every generator row
    is even, so extension gives distance at least38. No spectrum receipt or
    previous numerical shortening bound is used for this premise.
    """
    powers, value = [], 1
    for _ in range(255):
        powers.append(value)
        value = gf_mul(value, 2)
    if value != 1 or len(set(powers)) != 255 or 0 in powers:
        raise ArithmeticError('primitive degree8 field representation required')
    for root in powers[1:37]:
        value = 0
        for bit in range(P.bit_length()-1, -1, -1):
            value = gf_mul(value, root) ^ (P >> bit & 1)
        if value:
            raise ArithmeticError('BCH polynomial lacks a required consecutive root')
    columns, rank = syndrome_columns(rows, 256)
    if len(rows) != 128 or rank != 128:
        raise ArithmeticError('production binary dimension128 required')
    if any(word.bit_count() % 2 or remainder(word & ((1 << 255)-1), P) for word in rows):
        raise ArithmeticError('production code is not an even BCH extension')
    return columns


def near_full_words(columns, width, blocks):
    """Enumerate codewords obtained by deleting zero/two bits from full blocks.

    Completeness for *all* words in those blocks additionally requires an
    even code with minimum distance at least width*blocks-2. The production
    entry point checks that premise. Width>=3 prevents a deleted pair from
    emptying a selected block, so every result has exactly blocks active.
    """
    length = len(columns)
    if (type(width) is not int or width < 3 or length % width
            or type(blocks) is not int or not 1 <= blocks <= length//width
            or width*blocks % 2):
        raise ValueError('integral blocks of width>=3 and even selected length required')
    pairs = {}
    for a, b in combinations(range(length), 2):
        syndrome = columns[a] ^ columns[b]
        if not syndrome or syndrome in pairs:
            raise ValueError('no weight2/4 codewords required for unique pair syndromes')
        pairs[syndrome] = (1 << a) | (1 << b)
    masks = [((1 << width)-1) << (width*b) for b in range(length//width)]
    syndromes = []
    for start in range(0, length, width):
        value = 0
        for column in columns[start:start+width]:
            value ^= column
        syndromes.append(value)
    words = []
    for selected in combinations(range(len(masks)), blocks):
        syndrome, support = 0, 0
        for b in selected:
            syndrome ^= syndromes[b]
            support |= masks[b]
        if not syndrome:
            words.append(support)
        pair = pairs.get(syndrome)
        if pair is not None and pair & support == pair:
            words.append(support ^ pair)
    if len(set(words)) != len(words):
        raise ArithmeticError('a canonical support was counted twice')
    return words


def read_rows(header):
    raw = Path(header).read_bytes()
    source = raw.decode('utf-8')
    body = source.split('BchRows[128][4] = {', 1)[1].split('};', 1)[0]
    limbs = [int(value, 16) for value in re.findall(r'0x([0-9a-fA-F]+)ULL', body)]
    if len(limbs) != 512 or any(word >= 1 << 64 for word in limbs):
        raise ValueError('complete production128x256 generator table required')
    return [sum(limbs[4*i+j] << (64*j) for j in range(4)) for i in range(128)], hashlib.sha256(raw).hexdigest()


def check_production(header=None):
    if header is None:
        header = Path(__file__).resolve().parents[4]/'spin/src/kernels/generated/BchCircuit.h'
    rows, digest = read_rows(header)
    columns = verify_distance38(rows)
    words = near_full_words(columns, 8, 5)
    # No word can use fewer than5 octets because4*8<38. Two independent
    # words cannot share a40-position union: their three nonzero words
    # would have total weight>=3*38, but each union position contributes2.
    return dict(schema='canonical-octet-h5-exact-1',
        generator=dict(path=str(Path(header).resolve()), sha256=digest),
        dimension=128, length=256, block_width=8, minimum_distance=38,
        checked_block_sets=comb(32, 5), checked_pair_syndromes=comb(256, 2),
        nonzero_codewords_at_most5blocks=len(words),
        word_weight_counts=dict(sorted(Counter(word.bit_count() for word in words).items())),
        rank_cumulative_caps_at5=[15*len(words), 0, 0, 0],
        canonical_cdf_at5=15*len(words), words=[hex(word) for word in words],
        note='Exact local count for the production generator; no inner or whole-code claim.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--header', type=Path)
    args = parser.parse_args()
    print(json.dumps(check_production(args.header), indent=2))


if __name__ == '__main__':
    main()
