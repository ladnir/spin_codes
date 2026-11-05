"""Find support histories adapted to a sampled route, not codewords.

Fixed seed, bounded memory (<128MiB arrays), serial <=30s, no files.
Reports an exact shortened RS dimension; FLOAT count floors are diagnostic.
"""
import time
import numpy as np
from diagnose_random_pairing_confined import choose2,crossing_logs
from diagnose_random_pairing_projection import softplus2


def evaluate(route,masks,block):
    states=np.repeat(masks,block,axis=0)
    unions=states.copy();unions[1:] |= states[:-1]
    counts=np.bincount(route[unions],minlength=8192)
    B=int(np.count_nonzero(states[1:] & ~states[:-1]))
    dimension=int(np.maximum(counts-128,0).sum())
    T=len(route)-1;w=256;W=460*(1<<21)//1000
    I=T*w-B;z=max(0,(T+1)*w+B-W)
    ks,law=crossing_logs(512,w)
    logm=max(float(p)-choose2(512-w,int(k)) for k,p in zip(ks,law))
    # Z_outer >= q^dimension from words entirely supported on these unions.
    # Thus this is a LOWER bound on our UPPER-bound expression, not on the
    # actual failure probability. No inner maps or state values are tested.
    floor=1+T*logm+B-32*(B+z)+choose2(I,z)+32*dimension
    envelope=sum(softplus2(256-4096+31*int(b)) for b in counts)
    return counts,dict(block_updates=block,swaps=B,shortened_dimension=dimension,
                       relaxation_log_floor=round(floor,2),
                       route_enumerator_log_upper=round(envelope,2),
                       chosen_history_log_upper=round(floor-32*dimension+envelope,2),
                       not_a_codeword=True)


def run():
    began=time.monotonic();rng=np.random.default_rng(20261026)
    route=np.stack([rng.permutation(8192).astype(np.uint16) for _ in range(256)]).reshape(4096,512)
    for block in (1024,256,128,32,16):
        best=None
        for start in range(2):
            selected=rng.random(8192)<.5
            masks=np.zeros((4096//block,512),dtype=bool)
            for iteration in range(12):
                for j in range(len(masks)):
                    score=selected[route[j*block:(j+1)*block]].sum(axis=0)
                    # Deterministic random tie ordering, no seed search.
                    order=np.lexsort((rng.random(512),score))
                    masks[j]=False;masks[j,order[-256:]]=True
                counts,result=evaluate(route,masks,block)
                if best is None or result['shortened_dimension']>best['shortened_dimension']:
                    best=dict(result,iteration=iteration,start=start)
                selected=counts>128
                if time.monotonic()-began>25:raise TimeoutError('route adaptation budget')
        print(best,flush=True)
    print(dict(seconds=time.monotonic()-began,seed=20261026,no_inner_maps=True))


if __name__=='__main__':run()
