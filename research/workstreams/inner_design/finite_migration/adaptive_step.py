"""One-round IMT length study. Binary64 diagnostics, never full certificates.

Preserves the old t>=64 diagnostic maps; builds smaller RM-derived maps with
the same deterministic recipe, capped by the available nonconstant dimension.
The optional injection-shell transfer isolates proof slack from step size.
"""
import argparse
from fractions import Fraction
from functools import lru_cache
import json
import math
from pathlib import Path
import random

import numpy as np
import parameter_no_constant as nc
import mixing_rounds as mixing

grid = nc.grid
HERE = Path(__file__).resolve().parent


def capacity(t):
    if t < 8 or t & (t-1):
        raise ValueError('Expected a power of two, at least eight')
    m = t.bit_length()-1
    return min(20, m*(m+1)//2)


@lru_cache(maxsize=None)
def inner(t, s):
    if not t.bit_length() <= s <= capacity(t):
        raise ValueError('State outside this map family')
    if t >= 64:
        return nc.inner(t,s)
    m, minimum, maximum = t.bit_length()-1, t.bit_length(), capacity(t)
    monomials = [(i,j) for i in range(m) for j in range(i+1,m)]
    rng = random.Random(grid.SEED_A+t)
    quadratics = grid.maps.sample_independent_masks(rng,len(monomials),maximum-m)
    a = [v >> 1 for v in grid.maps.coordinate_columns(m,quadratics,monomials)]
    low_mask = (1 << minimum)-1
    choices = sorted(set(range(1 << minimum))-{v & low_mask for v in a})
    shift = rng.choice(choices) | (rng.getrandbits(maximum-minimum) << minimum)
    mask = (1 << s)-1
    a = [(v ^ shift) & mask for v in a]
    rng = random.Random(grid.SEED_B+t)
    b = rng.sample(range(1,1 << minimum),t)
    for dim in range(minimum,s):
        while True:
            proposed = [v | (rng.getrandbits(1) << dim) for v in b]
            if grid.maps.rank(grid.maps.generator_words(proposed,dim+1)) == dim+1:
                b = proposed
                break
    generators = grid.maps.generator_words(a,s)
    assert grid.maps.rank(generators) == grid.maps.rank(grid.maps.generator_words(b,s)) == s
    assert len(set(a)) == len(set(b)) == t and all(a) and all(b)
    counts = grid.maps.enumerate_spectrum(generators,s,t)
    assert counts[0] == 1 and counts[t] == 0
    spectrum = {w:n for w,n in enumerate(counts) if w and n}
    images = [sum(((c&q).bit_count() & 1) << j for j,c in enumerate(a)) for q in b]
    return dict(t=t,s=s,spectrum=spectrum,levels=sorted(spectrum),
                cancellation=[[v.bit_count(),(v ^ (1 << j)).bit_count()] for j,v in enumerate(images)],
                expansion_columns=a,feedback_columns=b,transvection_rounds=1)


def adjacent_retention(epochs, gap=0):
    """Exact mean retained-state mass before the second singleton feedback.

The two active regions have `gap` empty regions between them. Singleton
epochs are independent uniform in 0..epochs-1. One transvection per update.
"""
    if epochs < 1 or gap < 0:
        raise ValueError('Invalid region geometry')
    return 2*Fraction(1,2**(gap*epochs))*(1-Fraction(1,2**epochs))**2/epochs**2


def cancel_probability(t,s,epochs,gap=0):
    rho = adjacent_retention(epochs,gap)
    return rho/t+(1-rho)/((1 << s)-1)


def epoch_logs(record,tilts,refresh=False,sharp=False):
    if not sharp:
        return grid.wm.transfers(record,np.exp(tilts),None if refresh else 1)
    matrices = []
    spectrum = {int(w):v for w,v in record['spectrum'].items()}
    for tilt in tilts:
        zero,one,n = mixing.transfers(spectrum,[v[0] for v in record['cancellation']],
                                     [v[1] for v in record['cancellation']],
                                     math.exp(-math.exp(float(tilt))),0. if refresh else .5)
        matrices.append(np.array([zero,one]).reshape(2,n,n))
    with np.errstate(divide='ignore'):
        logs = np.log(np.array(matrices))
    return logs[:,0],logs[:,1]


def evaluate(record,b,exponent,refresh=False,sharp=False,log_replay=False):
    rows = 2**(exponent+1)//b
    t,s = record['t'],record['s']
    if rows < t or rows % t:
        raise ValueError('Nonintegral epoch geometry')
    # Broad, length-centered bank; includes the old grid's .1 spacing.
    tilts = np.arange(-180,1,dtype=float)/10
    rz,ra = grid.wm.regions(*epoch_logs(record,tilts,refresh,sharp),rows//t)
    coefficients = grid.wm.coefficients if log_replay else grid.coefficients
    moments = coefficients(rz,ra-math.log(rows//t),b)
    cutoff = b*rows//10
    values = np.minimum(0.,moments+cutoff*np.exp(tilts)[:,None])
    ix = np.argmin(values,axis=0)
    best = values[ix,np.arange(b+1)]
    if b == 256:
        # Apply the retained deterministic BCH inequalities, not a fitted spectrum.
        model = mixing.model
        weighted = {w: Fraction.from_float(math.exp(max(-700.,math.log(rows)+best[w])))
                    for w in model.base.WEIGHTS}
        bound,_,rest = model.base.bch_bound(weighted)
        margin = math.log2(bound.denominator)-math.log2(bound.numerator)
        dominant = 38
    else:
        counts,_ = grid.outer(b)
        terms = [math.log(rows*n)+best[w] for w,n in counts.items()]
        margin = -float(np.logaddexp.reduce(terms))/math.log(2)
        dominant = list(counts)[int(np.argmax(terms))]
    return dict(b=b,t=t,s=s,exponent=exponent,refresh=refresh,sharp=sharp,
                q1_margin_bits=margin,outer_rows=rows,epochs_per_region=rows//t,
                dominant_weight=dominant,dominant_tilt=float(tilts[ix[dominant]]),
                dominant_grid_edge=bool(ix[dominant] in (0,len(tilts)-1)),
                full_distance_proved=False)


def run(args):
    if args.output.exists():
        raise FileExistsError('Use a fresh output')
    cells, maps = [], {}
    for t in args.t:
        s = min(args.s,capacity(t))
        record = inner(t,s)
        maps[f't{t}_s{s}'] = record
        for b in args.b:
            for exponent in args.m:
                if 2**(exponent+1)//b < t:
                    continue
                for refresh in ([False,True] if args.refresh else [False]):
                    row = evaluate(record,b,exponent,refresh,args.sharp)
                    cells.append(row)
                    print(row,flush=True)
    sources = grid.ladder.candidate.sources()
    for b in args.b:
        if b != 256:
            _,path = grid.outer(b)
            sources[path.relative_to(grid.ROOT).as_posix()] = grid.maps.sha(path)
    sources[Path(__file__).resolve().relative_to(grid.ROOT).as_posix()] = grid.maps.sha(Path(__file__))
    grid.ladder.model.base.write_new(args.output,dict(status='BINARY64_IMT_ADAPTIVE_STEP_Q1',
        cells=cells,maps=maps,full_distance_proved=False,source_sha256=sources))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--b',type=int,nargs='+',default=[64,128])
    p.add_argument('--t',type=int,nargs='+',default=[16,32,64])
    p.add_argument('--s',type=int,default=20)
    p.add_argument('--m',type=int,nargs='+',default=[12,14,16,18,20])
    p.add_argument('--refresh',action='store_true')
    p.add_argument('--sharp',action='store_true')
    run(p.parse_args())
