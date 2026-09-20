"""Nested IMT expansion chains excluding the all-one word."""
import argparse
from contextlib import contextmanager
from functools import lru_cache
import math
from pathlib import Path
import random

import numpy as np
import parameter_three_type_dense as dense

full = dense.prior.prior.prior.base
grid = full.grid


@lru_cache(maxsize=3)
def chains(t):
    assert t in (64,128,256)
    m,minimum = t.bit_length()-1,t.bit_length()
    monomials = [(i,j) for i in range(m) for j in range(i+1,m)]
    rng = random.Random(grid.SEED_A+t)
    quadratics = grid.maps.sample_independent_masks(rng,len(monomials),20-m)
    # Drop the constant generator, not a coordinate. Adding a fixed column
    # shift gives each generator a fixed constant coefficient; independence
    # of their nonconstant parts still excludes the all-one word.
    columns = [column >> 1 for column in grid.maps.coordinate_columns(m,quadratics,monomials)]
    mask = (1 << minimum)-1
    excluded = {column & mask for column in columns}
    choices = sorted(set(range(1 << minimum))-excluded)
    shift = rng.choice(choices) | (rng.getrandbits(20-minimum) << minimum)
    a = tuple(column ^ shift for column in columns)
    _,b = grid.chains(t)
    return a,b


@lru_cache(maxsize=40)
def inner(t,s):
    assert t.bit_length() <= s <= 20
    aa,bb = chains(t)
    mask = (1 << s)-1
    a,b = [v & mask for v in aa],[v & mask for v in bb]
    generators = grid.maps.generator_words(a,s)
    assert grid.maps.rank(generators) == grid.maps.rank(grid.maps.generator_words(b,s)) == s
    assert len(set(a)) == len(set(b)) == t and all(a) and all(b)
    counts = grid.maps.enumerate_spectrum(generators,s,t)
    assert counts[0] == 1 and counts[t] == 0
    spectrum = {w:n for w,n in enumerate(counts) if w and n}
    assert sum(spectrum.values()) == (1 << s)-1
    images = [sum(((c&q).bit_count() & 1) << j for j,c in enumerate(a)) for q in b]
    cancellation = [[word.bit_count(),(word ^ (1 << j)).bit_count()] for j,word in enumerate(images)]
    return dict(t=t,s=s,spectrum=spectrum,levels=sorted(spectrum),cancellation=cancellation,
        expansion_columns=a,feedback_columns=b,transvection_rounds=1)


@contextmanager
def use():
    original = grid.inner
    full.prepare.cache_clear()
    grid.inner = inner
    try:
        yield
    finally:
        full.prepare.cache_clear()
        grid.inner = original


def probe(output,seed):
    assert not output.exists()
    model = grid.ladder.model
    saved = model.base.read(seed)
    model.authenticate(saved)
    assert saved['status'] == 'BINARY64_IMT_THREE_TYPE_DENSE_COVER'
    b,t,s,exponent = saved['geometry']
    boxes = sorted(saved['dense']['selected_boxes'],key=lambda x:x['own_log_bound'],reverse=True)[:4]
    with use():
        checker = dense.Dense(t,s,b,exponent,sorted(set(np.arange(-7.,1.51,.5))|{.75}))
        rows = []
        for box in boxes:
            point = full.typed.vertices(box['lower'],box['upper'],checker.length)[0]
            result = checker.evaluate(point,point)
            rows.append(dict(point=point.tolist(),margin_bits=-result['own_log_bound']/math.log(2),
                             witness=result['witness']))
            print('no-constant point',rows[-1],flush=True)
        q1 = grid.screen(checker.record,b,exponent,np.arange(-180,1,dtype=float)/10)
        sources = grid.ladder.candidate.sources()
        _,path = grid.outer(b)
        sources[path.relative_to(grid.ROOT).as_posix()] = grid.maps.sha(path)
        sources[seed.relative_to(grid.ROOT).as_posix()] = model.base.sha(seed)
        model.base.write_new(output,dict(status='BINARY64_IMT_NO_CONSTANT_POINT_PROBE',
            geometry=saved['geometry'],inner=checker.record,points=rows,q1=q1,
            full_distance_proved=False,source_sha256=sources))
        print('no-constant Q1',q1['q1_margin_bits'],flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--seed',type=Path,required=True)
    a = p.parse_args()
    probe(a.output.resolve(),a.seed.resolve())
