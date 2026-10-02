"""Research-only exact state-basis search; does not run encoder benchmarks."""
from __future__ import annotations
import itertools
import json
from pathlib import Path
import random
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import packet_inner_quadratic_extension as maps

BASE = (1,2,4,8,16,32,64,0x3500,0x1180,0x5980,0x2c80,0xcf80,0x9400,0xe400,0xe00,0xd980)
MASKS = tuple(i for i in range(64) if i.bit_count() <= 2)

def multiply(left, right):
    return tuple(xor(right[j] for j in range(len(right)) if row >> j & 1) for row in left)

def xor(values):
    result = 0
    for value in values:
        result ^= value
    return result

def inverse(rows):
    n = len(rows)
    aug = [row | (1 << (n+i)) for i, row in enumerate(rows)]
    for i in range(n):
        p = next(j for j in range(i,n) if aug[j] >> i & 1)
        aug[i], aug[p] = aug[p], aug[i]
        for j in range(n):
            if j != i and aug[j] >> i & 1:
                aug[j] ^= aug[i]
    return tuple(row >> n for row in aug)

def transpose(rows, width):
    return tuple(sum(((row >> j) & 1) << i for i,row in enumerate(rows)) for j in range(width))

def popcost(rows):
    return sum(max(0,row.bit_count()-1) for row in rows)

def circuit(rows, inputs):
    remaining = list(rows)
    gates = []
    next_index = inputs
    while True:
        pairs = {}
        for row in remaining:
            terms = [j for j in range(next_index) if row >> j & 1]
            for pair in itertools.combinations(terms,2):
                pairs[pair] = pairs.get(pair,0)+1
        if not pairs:
            break
        pair, score = max(pairs.items(), key=lambda x: (x[1],-x[0][0],-x[0][1]))
        if score < 2:
            break
        mask = sum(1 << j for j in pair)
        for i, row in enumerate(remaining):
            if row & mask == mask:
                remaining[i] ^= mask | (1 << next_index)
        gates.append(pair)
        next_index += 1
    return len(gates)+popcost(remaining), gates, remaining

def original(bits):
    rows, _, _ = maps.construction(bits)
    coefficients = []
    for row in rows:
        values = [(row >> i)&1 for i in range(64)]
        for step in (1,2,4,8,16,32):
            for i in range(64):
                if i & step:
                    values[i] ^= values[i^step]
        assert all(not values[i] for i in range(64) if i not in MASKS)
        coefficients.append(sum(values[mask] << i for i,mask in enumerate(MASKS)))
    return tuple(coefficients)

def objective(expansion, feedback):
    return popcost(expansion)+popcost(feedback)

def change(expansion, feedback, basis, i, j):
    e = tuple(row ^ (((row>>i)&1)<<j) for row in expansion)
    f = list(feedback); f[i] ^= f[j]
    p = list(basis); p[i] ^= p[j]
    return e, tuple(f), tuple(p)

def search(bits=20, restarts=32):
    g = original(bits)
    initial = BASE + tuple(1<<i for i in range(16,bits))
    e0 = multiply(transpose(g,22), inverse(initial))
    f0 = multiply(initial, g)
    rng = random.Random(731)
    candidates = []
    for restart in range(restarts):
        e,f,p = e0,f0,initial
        for _ in range(restart % 12):
            i,j = rng.sample(range(bits),2)
            e,f,p = change(e,f,p,i,j)
        for _ in range(200):
            before = objective(e,f)
            best = before
            choices = []
            for i in range(bits):
                for j in range(bits):
                    if i==j: continue
                    candidate = change(e,f,p,i,j)
                    value = objective(candidate[0],candidate[1])
                    if value < best:
                        choices = [candidate];best = value
                    elif value == best and value < before:
                        choices.append(candidate)
            if not choices: break
            e,f,p = rng.choice(choices)
        ec,eg,er = circuit(e,bits)
        fc,fg,fr = circuit(f,22)
        candidates.append((ec+fc, objective(e,f),p,e,f,(ec,eg,er),(fc,fg,fr)))
    result=min(candidates,key=lambda item:item[:3])
    _,_,p,e,f,ec,fc=result
    assert multiply(p,inverse(p)) == tuple(1<<i for i in range(bits))
    assert e == multiply(transpose(g,22),inverse(p))
    assert f == multiply(p,g)
    return dict(bits=bits, basis=list(p),inverse_basis=list(inverse(p)),
                initial_popcost=objective(e0,f0),initial_greedy_cse=circuit(e0,bits)[0]+circuit(f0,22)[0],
                popcost=objective(e,f),greedy_cse=result[0],expansion_cost=ec[0],feedback_cost=fc[0],
                expansion=list(e),feedback=list(f),anf_masks=list(MASKS),
                expansion_gates=ec[1],expansion_outputs=ec[2],feedback_gates=fc[1],feedback_outputs=fc[2])

if __name__ == '__main__':
    print(json.dumps(search(int(sys.argv[1]) if len(sys.argv)>1 else 20),indent=2))
