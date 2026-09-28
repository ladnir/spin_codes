"""Different distribution: uniformly shuffle all 16 positions in each bundle.

The four source columns still select the same macroregion and epoch window.
Only the independent permutation inside that window changes.
"""
from collections import Counter,defaultdict
from math import comb


def lift(data):
    spectrum,allowed,moments,atoms,cancel=data
    by_weight=defaultdict(list)
    for s in allowed:
        by_weight[sum(s)].append(s)
    masks={};mom={};atom={};can={}
    for a,shapes in by_weight.items():
        masks[a]=[x for s in shapes for x in allowed[s]]
        assert len(masks[a])==len(set(masks[a]))==comb(16,a)
        atom[a]=Counter();mom[a]=defaultdict(Counter);can[a]=defaultdict(Counter)
        for s in shapes:
            atom[a].update(atoms[s])
            for v,row in moments[s].items():
                mom[a][v].update(row)
            for v,row in cancel[s].items():
                can[a][v].update(row)
        choices=8*comb(16,a)
        assert sum(atom[a].values())==choices
        assert sum(sum(row.values()) for row in can[a].values())+atom[a][0]==choices
        for v,count in spectrum.items():
            assert sum(mom[a][v].values())==choices*count
    assert [atom[a][0] for a in range(1,17)]==[0]*7+[2]+[0]*8
    print('Full-window ensemble: all 16-subset orbit counts reconstructed exactly',flush=True)
    return (spectrum,{s:masks[sum(s)] for s in allowed},
            {s:mom[sum(s)] for s in allowed},{s:atom[sum(s)] for s in allowed},
            {s:can[sum(s)] for s in allowed})
