"""Generate compact, UNTRUSTED exact-count candidates; no Walsh transform.

Only forward Map22 images and 240/27000 one/two-packet inputs are visited.
No per-state arrays are serialized. These counts still require a Lean
relation to the actual map before they can authenticate physical rows.
"""
from __future__ import annotations
import argparse
from collections import Counter
import hashlib
from itertools import combinations, product
import json
from math import comb
from pathlib import Path
import re
import subprocess
import sys
import time
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
MAP = 'e2467b9d5de8403ccd8263f0f55c73df128a8a2dd4c4cc405bc3a8e057fc3581'


def run(args):
    if args.output.exists():
        raise FileExistsError('fresh candidate path required')
    start = time.monotonic()
    physical = json.loads(args.physical_receipt.read_text(encoding='utf-8'))
    record = physical['map_record']
    rows = [int(x, 16) for x in record['expansion_rows_hex']]
    lean = ROOT / 'research/lean/SpinCodes/UniformMixer/Map22.lean'
    text = lean.read_text(encoding='utf-8')
    literal = text.split('def rowMask : Fin 22 → ℕ := ![', 1)[1].split(']', 1)[0]
    lean_rows = [int(x, 16) for x in re.findall(r'0x[0-9a-fA-F]+', literal)]
    columns = [sum(((row >> p) & 1) << i for i, row in enumerate(rows)) for p in range(64)]
    identity = dict(bits=22, width=64, expansion_rows=list(map(hex, rows)), feedback_columns=columns)
    digest = hashlib.sha256(json.dumps(identity, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
    if rows != lean_rows or len(rows) != 22 or columns != record['feedback_columns'] or digest != MAP:
        raise ValueError('actual Lean Map22 and retained producer map must agree exactly')
    size = 1 << 22
    images = np.zeros(size, dtype=np.uint64)
    width = 1
    for row in rows:
        images[width:2*width] = images[:width] ^ np.uint64(row)
        width *= 2
    codes = np.zeros(size, dtype=np.uint32)
    lookup = np.array([17 ** i.bit_count() for i in range(16)], dtype=np.uint32)
    for offset in range(0, 64, 4):
        codes += lookup[((images >> offset) & 15).astype(np.uint8)]
    unique, counts = np.unique(codes, return_counts=True)
    profiles = [[int(code) // 17**r % 17 for r in range(5)] for code in unique]
    multiplicities = list(map(int, counts))
    if len(profiles) != 25 or sum(multiplicities) != size or any(sum(h) != 16 for h in profiles):
        raise ArithmeticError('complete 25-profile partition expected')
    levels = [0, 16, 24, 28, 32, 36, 40, 48, 64]
    if sorted({sum(r*n for r,n in enumerate(h)) for h in profiles}) != levels:
        raise ArithmeticError('production weight classes mismatch')
    feedback = []
    for slot in range(16):
        values = []
        for label in range(1, 16):
            state = 0
            for b in range(4):
                if label >> b & 1:
                    state ^= columns[4*slot+b]
            values.append((state, label.bit_count()))
        feedback.append(values)
    image_cache = {0: 0}
    def expansion_weight(state):
        if state not in image_cache:
            mask, image = state, 0
            while mask:
                bit = (mask & -mask).bit_length() - 1
                image ^= rows[bit]
                mask &= mask - 1
            image_cache[state] = image.bit_count()
        return image_cache[state]
    births = []
    for j in (1, 2):
        table = Counter()
        for slots in combinations(range(16), j):
            for labels in product(*(feedback[s] for s in slots)):
                state, weight = 0, 0
                for syndrome, w in labels:
                    state ^= syndrome
                    weight += w
                table[expansion_weight(state), weight] += 1
        denominator = comb(16, j) * 15**j
        if sum(table.values()) != denominator or any(w not in levels for w, _ in table):
            raise ArithmeticError('small exact-nonzero birth census mismatch')
        births.append(dict(occupancy=j, denominator=denominator,
            expansion_weight_levels=levels,
            coefficients=[[table[level, w] for w in range(4*j+1)] for level in levels]))
    result = dict(schema='turbo-coarse-packet-profile-candidate-1',
        scope='Untrusted compact exact integer candidates; not Lean-authenticated counts or a probability certificate.',
        map_sha256=MAP, physical_receipt_sha256=hashlib.sha256(args.physical_receipt.read_bytes()).hexdigest(),
        lean_map_sha256=hashlib.sha256(lean.read_bytes()).hexdigest(),
        expansion_rows_hex=list(map(hex, rows)), profiles=profiles, multiplicities=multiplicities,
        state_count=size, nonzero_state_count=size-1, packet_slots=16, small_births=births,
        walsh_transforms_performed=0, raw_state_arrays_persisted=False,
        elapsed_seconds=time.monotonic()-start)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x', encoding='utf-8') as out:
        json.dump(result, out, indent=2)
        out.write('\n')
    print(json.dumps(dict(profiles=len(profiles), birth_inputs=[b['denominator'] for b in births],
        output_bytes=args.output.stat().st_size, elapsed_seconds=result['elapsed_seconds'])))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--physical-receipt', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--max-seconds', type=float, default=30)
    parser.add_argument('--_worker', action='store_true', help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args._worker:
        run(args)
    else:
        result = subprocess.run([sys.executable, str(Path(__file__).resolve()),
            *sys.argv[1:], '--_worker'], timeout=args.max_seconds)
        raise SystemExit(result.returncode)
