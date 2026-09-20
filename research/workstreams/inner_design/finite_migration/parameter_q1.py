"""Fresh IMT Q1 diagnostics for every geometry in the paper's three slices.

These are nearest-binary64 bounds on Q1 only, not full distance certificates.
The RM-derived expansion chain is shared with the historical study; feedback
is a separately generated nested chain and mixing is one transvection.
"""
import argparse
from functools import lru_cache
import hashlib
import json
import math
from pathlib import Path
import random
import sys

import numpy as np
import ladder

ROOT = ladder.model.ROOT
sys.path.insert(0, str(ROOT/'workstreams/finite_asymptotic_theory/landscape_db'))
import run_activation_pilot as maps
wm = ladder.model.independent.search.wm
SEED_A, SEED_B = 3390173185, 2026091601


def geometries():
    keys = {(b, 64, 20, m) for b in (64, 128) for m in range(12, 27)}
    keys |= {(b, t, s, 20) for b in (64, 128) for t in (64, 128, 256)
             for s in range(t.bit_length(), 21)}
    keys |= {(b, 64, s, m) for b in (64, 128) for s in (10, 12, 16)
             for m in (16, 18, 20, 22, 24)}
    assert len(keys) == 130
    return sorted(keys, key=lambda x: (x[1], x[2], x[0], x[3]))


@lru_cache(maxsize=3)
def chains(t):
    assert t in (64, 128, 256)
    m, minimum = t.bit_length()-1, t.bit_length()
    monomials = [(i,j) for i in range(m) for j in range(i+1,m)]
    quadratics = maps.sample_independent_masks(random.Random(SEED_A+t), len(monomials), 20-m-1)
    a = maps.coordinate_columns(m, quadratics, monomials)
    rng = random.Random(SEED_B+t)
    b = rng.sample(range(1, 1 << minimum), t)
    assert maps.rank(maps.generator_words(b, minimum)) == minimum
    for s in range(minimum, 20):
        while True:
            proposed = [column | (rng.getrandbits(1) << s) for column in b]
            if maps.rank(maps.generator_words(proposed, s+1)) == s+1:
                b = proposed
                break
    return tuple(a), tuple(b)


@lru_cache(maxsize=40)
def inner(t, s):
    assert t.bit_length() <= s <= 20
    aa, bb = chains(t)
    mask = (1 << s)-1
    a, b = [c & mask for c in aa], [c & mask for c in bb]
    generators = maps.generator_words(a, s)
    assert maps.rank(generators) == maps.rank(maps.generator_words(b,s)) == s
    assert len(set(a)) == len(set(b)) == t and all(a) and all(b)
    counts = maps.enumerate_spectrum(generators, s, t)
    spectrum = {w:n for w,n in enumerate(counts) if w and n}
    assert sum(spectrum.values()) == (1 << s)-1
    images = [sum(((c&q).bit_count() & 1) << j for j,c in enumerate(a)) for q in b]
    cancellation = [[word.bit_count(), (word ^ (1 << j)).bit_count()] for j,word in enumerate(images)]
    return dict(t=t, s=s, spectrum=spectrum, levels=sorted(spectrum), cancellation=cancellation,
                expansion_columns=a, feedback_columns=b, transvection_rounds=1)


def coefficients(rz, ra, length):
    """Scaled positive batched products; one scale per tilt and coefficient."""
    if not np.all(np.isfinite(rz) | np.isneginf(rz)) or not np.all(np.isfinite(ra) | np.isneginf(ra)):
        raise ArithmeticError('Invalid region logarithms')
    scales = [x.max(axis=(1,2)) for x in (rz,ra)]
    matrices = [np.exp(x-scale[:,None,None]) for x,scale in zip((rz,ra), scales)]
    batch, n = len(rz), rz.shape[1]
    current = np.zeros((batch, length+1, n))
    current[:,0,0] = 1
    logs = np.full((batch, length+1), -np.inf)
    logs[:,0] = 0
    for step in range(length):
        count = step+1
        left = current[:,:count] @ matrices[0]
        right = current[:,:count] @ matrices[1]
        llogs, rlogs = logs[:,:count]+scales[0][:,None], logs[:,:count]+scales[1][:,None]
        merged = np.full((batch,count+1), -np.inf)
        merged[:,:count] = llogs
        merged[:,1:] = np.maximum(merged[:,1:], rlogs)
        updated = np.zeros((batch,count+1,n))
        updated[:,:count] = left*np.exp(llogs-merged[:,:count])[:,:,None]
        updated[:,1:] += right*np.exp(rlogs-merged[:,1:])[:,:,None]
        normalizers = updated.max(axis=2)
        if np.any(normalizers <= 0) or not np.all(np.isfinite(normalizers)):
            # Never interpret a vanished coefficient as a useful tiny bound.
            raise ArithmeticError('Vanished or invalid diagnostic coefficient')
        current[:,:count+1] = updated/normalizers[:,:,None]
        logs[:,:count+1] = merged+np.log(normalizers)
    return (logs+np.log(current.sum(axis=2))
            -np.array([math.log(math.comb(length,w)) for w in range(length+1)]))


@lru_cache(maxsize=2)
def outer(b):
    key = {64:'xbch64', 128:'ebch128'}[b]
    c = maps.CONSTITUENTS[key]
    path = c.spectrum_path
    if path.exists():
        counts = maps.load_spectrum(c)
    else:
        path = maps.SMALL/'spectra'/f'{key}_weight_counts.json'
        record = json.loads(path.read_text())
        assert (record['length'], record['dimension']) == (b,b//2)
        counts = {int(w):int(n) for w,n in record['weight_counts'].items()}
    assert counts[0] == 1 and sum(counts.values()) == 1 << (b//2)
    assert min(w for w,n in counts.items() if w and n) == c.minimum_distance
    return {w:n for w,n in counts.items() if w and n}, path


def screen(record, b, exponent, tilts):
    rows = (1 << exponent)//(b//2)
    t = record['t']
    assert rows % t == 0
    lam = np.exp(tilts)
    zero, one = wm.transfers(record,lam,1)
    rz, ra = wm.regions(zero,one,rows//t)
    moments = coefficients(rz,ra-math.log(rows//t),b)
    cutoff = b*rows//10
    values = np.minimum(0., moments+cutoff*lam[:,None])
    witnesses = np.argmin(values,axis=0)
    best = values[witnesses,np.arange(b+1)]
    counts,_ = outer(b)
    terms = [math.log(rows*n)+best[w] for w,n in counts.items()]
    dominant = list(counts)[int(np.argmax(terms))]
    return dict(b=b,t=t,s=record['s'],exponent=exponent,
        q1_margin_bits=-float(np.logaddexp.reduce(terms))/math.log(2),
        full_margin_bits=None, dominant_weight=dominant,
        dominant_tilt=float(tilts[witnesses[dominant]]),
        dominant_grid_edge=bool(witnesses[dominant] in (0,len(tilts)-1)),
        bad_weight=cutoff,outer_rows=rows, full_distance_proved=False)


def run(output, limit=None):
    assert not output.exists()
    tilts = np.arange(-180,1,dtype=float)/10
    selected = geometries() if limit is None else geometries()[:limit]
    records, identities = [], {}
    for b,t,s,m in selected:
        record = inner(t,s)
        identities[f't{t}_s{s}'] = record
        result = screen(record,b,m,tilts)
        records.append(result)
        print('IMT Q1',len(records),'/',len(selected),b,t,s,m,result['q1_margin_bits'],flush=True)
    sources = ladder.candidate.sources()
    for b in (64,128):
        _, path = outer(b)
        sources[path.relative_to(ROOT).as_posix()] = maps.sha(path)
    ladder.model.base.write_new(output,dict(status='BINARY64_IMT_PARAMETER_Q1_ONLY',
        rows=records, maps=identities,geometry_complete=selected == geometries(),
        full_distance_proved=False,tilts=list(tilts),source_sha256=sources,
        expansion_chain_seed=SEED_A,feedback_chain_seed=SEED_B))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--limit',type=int)
    a = p.parse_args()
    if a.limit is not None and not 1 <= a.limit <= 130:
        p.error('limit must be in 1..130')
    run(a.output.resolve(),a.limit)
