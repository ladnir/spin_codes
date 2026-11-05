"""Exact checks of swap deltas before full-size structural search."""
import shutil
import tempfile
import time
import numpy as np
from diagnose_random_pairing_rank_search import block_marginals,improve_block


def cost(target,partners,masks,block):
    states=np.repeat(masks,block,axis=0)
    previous=np.zeros_like(states);previous[1:]=states[:-1]
    visible=previous|states|np.take_along_axis(previous,partners,axis=1)
    return int(np.count_nonzero(target&visible)+np.count_nonzero(previous&~states))


def main():
    assert shutil.disk_usage(tempfile.gettempdir()).free>=20*(1<<30),'disk guard'
    began=time.monotonic();rng=np.random.default_rng(20261027);checks=steps=0
    for block in (1,2,4):
        L=8;n=8
        for fixture in range(10):
            partners=np.empty((L,n),dtype=np.int32)
            for t in range(L):
                order=rng.permutation(n);partners[t,order]=order[np.arange(n)^1]
            target=rng.random((L,n))<.3
            masks=np.zeros((L//block,n),dtype=bool)
            for H in masks:H[rng.permutation(n)[:n//2]]=True
            for j in range(len(masks)):
                baseline=cost(target,partners,masks,block)
                flip,ps,ws,lastp,lastw=block_marginals(target,partners,masks,j,block)
                for a in np.flatnonzero(masks[j]):
                    for b in np.flatnonzero(~masks[j]):
                        correction=int(np.sum((ps[:,a]==b)*ws[:,a]))
                        if lastp is not None:correction+=int(lastp[a]==b)*int(lastw[a])
                        trial=masks.copy();trial[j,a]=False;trial[j,b]=True
                        assert cost(target,partners,trial,block)-baseline==flip[a]+flip[b]+correction
                        checks+=1
                changed=improve_block(target,partners,masks,j,block)
                actual=cost(target,partners,masks,block)
                assert actual<baseline if changed else actual==baseline
                assert np.all(masks.sum(axis=1)==n//2)
                steps+=1
    assert time.monotonic()-began<25
    print(dict(exact_swap_checks=checks,descent_checks=steps,seconds=time.monotonic()-began))


if __name__=='__main__':main()
