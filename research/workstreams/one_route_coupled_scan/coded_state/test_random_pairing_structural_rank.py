"""Exact small checks of the support-restricted rank obstruction.

Uses repetition outers and local invertible GF2 gates, not the accepted
SL2 ensemble. Tests the coefficient-independent dimension argument only.
Serial <=30s, no output files.
"""
from itertools import combinations,product
import random
import shutil
import tempfile
import time


def main():
    assert shutil.disk_usage(tempfile.gettempdir()).free>=20*(1<<30),'disk guard'
    began=time.monotonic();rng=random.Random(20261026)
    masks=[sum(1<<j for j in c) for c in combinations(range(4),2)]
    pairings=(((0,1),(2,3)),((0,2),(1,3)),((0,3),(1,2)))
    identities=witness_checks=0
    for fixture in range(3):
        route=list(range(8))*2
        if fixture:
            left=list(range(8));right=left.copy();rng.shuffle(left);rng.shuffle(right)
            route=left+right
        gates=[pairings[0] if not fixture else rng.choice(pairings) for _ in range(4)]
        words=[]
        for message in range(1,256):
            Q=0;states=[];outputs=[]
            for step,pairs in enumerate(gates):
                mixed=0
                for a,b in pairs:
                    mixed|=(((Q>>a)^(Q>>b))&1)<<a
                    mixed|=((Q>>b)&1)<<b
                X=sum(((message>>route[step*4+j])&1)<<j for j in range(4))
                next_state=mixed^X;outputs.append(next_state^Q)
                states.append(next_state);Q=next_state
            words.append((states,outputs))
        for states in product(masks,repeat=4):
            previous=0;coverage=[0]*8;J=D=E=0;intersections=[]
            for step,(H,pairs) in enumerate(zip(states,gates)):
                live=0
                for a,b in pairs:
                    pair=(1<<a)|(1<<b)
                    if previous&pair:live|=pair
                union=previous|H;visible=live|H
                J+=(live&~union).bit_count();D+=(previous&~H).bit_count()
                E+=visible.bit_count()
                for j in range(4):
                    if visible>>j&1:coverage[route[4*step+j]]+=1
                    if (previous&H)>>j&1:intersections.append((step,j))
                previous=H
            A=sum(H.bit_count() for H in states)
            assert E==A+J+D;identities+=1
            # Each repetition row shortened to its two available positions.
            d=sum(max(c-1,0) for c in coverage)
            if d-J-D<=0:continue
            z=min(len(intersections),d-J-D-1)
            selected=intersections[:z]
            assert any(all(not (Q&~H) for Q,H in zip(qs,states))
                       and all(not ((ys[t]>>j)&1) for t,j in selected)
                       for qs,ys in words),(fixture,states,d,J,D,z)
            witness_checks+=1
        assert time.monotonic()-began<25
    assert witness_checks>0
    print(dict(visible_equation_identities=identities,positive_gap_witness_checks=witness_checks,
               seconds=time.monotonic()-began,full_size_obstruction_found=False))


if __name__=='__main__':main()
