"""One-active-group output-weight diagnostic, including returns to zero.

Binary64 exploration only until replayed with outward arithmetic. It does not
cover messages spanning multiple groups. No files or production code changed.
"""

from collections import Counter, defaultdict
from itertools import combinations
from math import comb, log
import argparse
import json
import random
import re

import numpy as np
from scipy.special import logsumexp

import bch_joint_support as joint


def maps():
    source = (joint.REPO / 'spin/src/kernels/generated/SelectedMaps.h').read_text()
    source = source.split('struct Map128S19 {',1)[1]
    def columns(name):
        raw = re.search(r'\b' + name + r'\{([^}]+)\}',source).group(1)
        return [int(v.strip(),0) for v in raw.split(',')]
    a, b = columns('columns'), columns('feedbackColumns')
    record = json.loads((joint.RESEARCH / 'workstreams/inner_design/NO_CONSTANT_MAP.json').read_text())
    assert a == record['columns']
    pool = [sum(1 << j for j in supp) for supp in combinations(range(19),5)]
    assert b == random.Random(0).sample(pool,128)
    rows = [sum(((c >> j) & 1) << p for p,c in enumerate(a)) for j in range(19)]
    images = [0] * (1 << 19)
    for state in range(1,1 << 19):
        images[state] = images[state & (state-1)] ^ rows[(state & -state).bit_length()-1]
    spectrum = Counter(image.bit_count() for image in images)
    assert spectrum == {int(w):n for w,n in record['spectrum'].items()}
    assert spectrum[0] == 1
    del spectrum[0]
    return images,b,spectrum


def census():
    images,b,spectrum = maps()
    # Number of (nonzero state, window) pairs by total and window weight.
    windows = defaultdict(Counter)
    for image in images[1:]:
        weight = image.bit_count()
        for start in range(0,128,16):
            windows[weight][((image >> start) & 65535).bit_count()] += 1
    assert all(sum(windows[w].values()) == 8*n for w,n in spectrum.items())
    zero = [0] * 17
    atoms = [Counter() for _ in range(17)]
    cancellation = [defaultdict(Counter) for _ in range(17)]
    for start in range(0,128,16):
        syndrome = [0] * 65536
        for mask in range(1,65536):
            a = mask.bit_count()
            syndrome[mask] = syndrome[mask & (mask-1)] ^ b[start+(mask & -mask).bit_length()-1]
            value = syndrome[mask]
            atoms[a][value] += 1
            if value:
                image = images[value]
                cancellation[a][image.bit_count()][(image ^ (mask << start)).bit_count()] += 1
            else:
                zero[a] += 1
    assert zero == [0]*8 + [2] + [0]*8
    maximum = [0]+[max(atoms[a].values()) for a in range(1,17)]
    for a in range(1,17):
        assert sum(atoms[a].values()) == 8*comb(16,a)
        assert sum(sum(v.values()) for v in cancellation[a].values()) + zero[a] == 8*comb(16,a)
    print(f'Exact map/cancellation census: expansion levels {sorted(spectrum)}, maximum atoms {maximum[1:]}',flush=True)
    return spectrum,windows,zero,maximum,cancellation


def transfers(data, lam, active_weight=None):
    spectrum,windows,zero,maximum,cancellation = data
    levels = sorted(spectrum)
    n,m = len(levels)+2,(1 << 19)-1
    z = np.exp(-lam)
    powers = z**np.arange(145)
    inactive = np.zeros((n,n))
    inactive[0,0] = 1
    masses = np.array([spectrum[w]/m for w in levels])
    inactive[1,1] = powers[min(levels)]/2
    inactive[1,2:] = powers[min(levels)]*masses/2
    for i,v in enumerate(levels):
        inactive[i+2,i+2] += powers[v]/2
        inactive[i+2,2:] += powers[v]*masses/2
    matrices = []
    for a in ([active_weight] if active_weight else range(1,17)):
        choices = 8*comb(16,a)
        matrix = np.zeros((n,n))
        matrix[0,0] = zero[a]/choices*powers[a]
        matrix[0,1] = (1-zero[a]/choices)*powers[a]
        # Unknown nonzero state: deterministic triangle-inequality bound.
        d = powers[max(0,min(levels)-a)]
        min_cancel = min(w for row in cancellation[a].values() for w in row)
        cd = min(d,maximum[a]/choices*powers[min_cancel])
        matrix[1,0] = cd/2+d/(2*m)
        matrix[1,1] = d/2
        matrix[1,2:] = d*masses/2
        for i,v in enumerate(levels):
            # Level rows represent a pointwise uniform-density envelope, not
            # an assumption that the true conditioned state is uniform.
            moment = 0.
            for local,count in windows[v].items():
                for overlap in range(max(0,a+local-16),min(local,a)+1):
                    moment += count*comb(local,overlap)*comb(16-local,a-overlap)*powers[v+a-2*overlap]
            moment /= choices*spectrum[v]
            cancel = sum(count*powers[w] for w,count in cancellation[a][v].items())/(choices*spectrum[v])
            matrix[i+2,0] = cancel/2+moment/(2*m)
            matrix[i+2,1] = moment/2
            matrix[i+2,2:] = moment*masses/2
        matrices.append(matrix)
    # This permits adversarial active-lane counts across regions. It only
    # enlarges the nonnegative envelope for any fixed message's column types.
    return inactive,np.maximum.reduce(matrices)


def regions(zero,one,epochs=64):
    rz = np.eye(len(zero))
    ra = np.zeros_like(zero)
    for _ in range(epochs):
        ra = ra @ zero + rz @ one
        rz = rz @ zero
    return rz,ra/epochs


def log_moments(zero,one,length=256):
    # Log-domain polynomial recurrence: no underflow for high support ranks.
    with np.errstate(divide='ignore'):
        z,a = np.log(zero),np.log(one)
    current = np.full((1,len(zero)),-np.inf)
    current[0,0] = 0
    for _ in range(length):
        empty = logsumexp(current[:,:,None]+z[None,:,:],axis=1)
        active = logsumexp(current[:,:,None]+a[None,:,:],axis=1)
        updated = np.full((len(current)+1,len(zero)),-np.inf)
        updated[:-1] = empty
        updated[1:] = np.logaddexp(updated[1:],active)
        current = updated
    return logsumexp(current,axis=1)-np.array([log(comb(length,u)) for u in range(length+1)])


def log_weighted_caps(caps, log_weights):
    # Taking the suffix maximum yields a decreasing majorant, so summation
    # by parts with upper CDFs remains valid even if the raw curve increases.
    weights = np.maximum.accumulate(log_weights[::-1])[::-1]
    logs = []
    for u in range(len(caps)-1):
        if caps[u] and weights[u] > weights[u+1]:
            difference = weights[u]+np.log(-np.expm1(weights[u+1]-weights[u]))
            logs.append(log(caps[u])+difference)
    logs.append(log(caps[-1])+weights[-1])
    return logsumexp(logs)


def first_support_moments(zero,one,length=256):
    """Conditional log moments by remaining regions l and total support u.

    The first occupied region has transfer one. The preceding regions have
    zero input and zero state. The remaining u-1 occupied regions are uniform
    among the final l regions. Entry [l,u] stores this conditional moment.
    """
    with np.errstate(divide='ignore'):
        z,a = np.log(zero),np.log(one)
    current = np.zeros((1,len(zero)))  # terminal column, all states
    result = np.full((length,length+1),np.inf)
    for remaining in range(length):
        values = logsumexp(current+a[0,:],axis=1)
        result[remaining,1:remaining+2] = values-np.array([log(comb(remaining,w)) for w in range(remaining+1)])
        empty = logsumexp(z[None,:,:]+current[:,None,:],axis=2)
        active = logsumexp(a[None,:,:]+current[:,None,:],axis=2)
        updated = np.full((len(current)+1,len(zero)),-np.inf)
        updated[:-1] = empty
        updated[1:] = np.logaddexp(updated[1:],active)
        current = updated
    return result


def conditioned_rank_one(data,spectrum):
    subset_logs = []
    for a in range(1,17):
        best = np.full((256,257),np.inf)
        for lam in (.00016,.00025,.0004,.00064,.001,.0016,.0025):
            moment = first_support_moments(*regions(*transfers(data,lam,a)))
            best = np.minimum(best,moment+lam*209715)
        subset_logs.append(log(comb(16,a))+best)
        print(f'Conditioned rank-one diagnostic: processed active-row count {a}',flush=True)
    union = np.minimum(0,logsumexp(subset_logs,axis=0))
    terms = []
    for w in range(1,257):
        if not spectrum[w]:
            continue
        # P(first occupied region leaves l following regions).
        conditional = [log(comb(l,w-1))-log(comb(256,w))+union[l,w]
                       for l in range(w-1,256)]
        terms.append((w,log(512*spectrum[w])+logsumexp(conditional)))
    print(f'BINARY64 ONLY shared-placement rank1 log2 union={logsumexp([v for w,v in terms])/log(2):.3f}',flush=True)
    print('Dominant shells:',[(w,round(v/log(2),3)) for w,v in sorted(terms,key=lambda p:-p[1])[:5]],flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--rank-one',action='store_true')
    parser.add_argument('--condition-first',action='store_true')
    args = parser.parse_args()
    spectrum = joint.authenticated_caps()
    data = census()
    if args.condition_first:
        conditioned_rank_one(data,spectrum)
        return
    if args.rank_one:
        total_logs = []
        for a in range(1,17):
            best = np.full(257,np.inf)
            for lam in (.00016,.0002,.00025,.00032,.0004,.0005,.00064):
                moments = log_moments(*regions(*transfers(data,lam,a)))
                best = np.minimum(best,np.minimum(0,moments+lam*209715))
            upper = log(512*comb(16,a))+logsumexp([log(spectrum[w])+best[w]
                                                     for w in range(1,257) if spectrum[w]])
            total_logs.append(upper)
            print(f'BINARY64 ONLY rank1 active rows={a}: log2 union={upper/log(2):.3f}',flush=True)
        print(f'BINARY64 ONLY rank1 all row subsets: log2 union={logsumexp(total_logs)/log(2):.3f}',flush=True)
        return
    caps = joint.support_caps(spectrum,dimensions=joint.dimension_caps())
    best = np.full((16,257),np.inf)
    for lam in (.0001,.0002,.0004,.0008,.0016,.0032,.0064,.0128,.0256):
        moments = log_moments(*regions(*transfers(data,lam)))
        coefficients = np.minimum(0,moments+lam*(209715))
        # A separate tilt may be chosen for every fixed support class.
        best = np.minimum(best,coefficients)
        by_rank = [log(512)+log_weighted_caps(row,best[h]) for h,row in enumerate(caps)]
        print(f'BINARY64 ONLY tilt={lam:.5g} log2 union={logsumexp(by_rank)/log(2):.3f} '
              f'rank1={by_rank[0]/log(2):.3f} rank16={by_rank[-1]/log(2):.3f}',flush=True)
    print('Diagnostic only: no outward replay or multi-group coverage.')


if __name__ == '__main__':
    main()
