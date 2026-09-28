"""Exact local orbit counts for four-row bundles of one or two columns."""
from collections import Counter,defaultdict
from group_moment import maps
from rank_two_moment import intersection_counts


def census(columns_per_bundle=2):
    assert columns_per_bundle in (1,2)
    width=4*columns_per_bundle
    mask_count=1<<width
    window_count=128//width
    images,columns,spectrum=maps()
    shapes=[tuple(sorted(((mask>>(4*j))&15).bit_count() for j in range(columns_per_bundle))) for mask in range(mask_count)]
    allowed=defaultdict(list)
    for mask,s in enumerate(shapes):
        if mask:
            allowed[s].append(mask)
    windows=defaultdict(Counter)
    for image in images[1:]:
        v=image.bit_count()
        for start in range(0,128,width):
            windows[v][shapes[(image>>start)&(mask_count-1)]]+=1
    moments={s:defaultdict(Counter) for s in allowed}
    for state_shape in set(shapes):
        for s in allowed:
            counts=intersection_counts(state_shape,s)
            assert sum(counts.values())==len(allowed[s])
            for v in spectrum:
                frequency=windows[v][state_shape]
                for overlap,count in counts.items():
                    moments[s][v][v+sum(s)-2*overlap]+=frequency*count
    atoms=defaultdict(Counter);cancel={s:defaultdict(Counter) for s in allowed}
    for start in range(0,128,width):
        syndromes=[0]*mask_count
        for mask in range(1,mask_count):
            bit=mask&-mask
            q=syndromes[mask^bit]^columns[start+bit.bit_length()-1]
            syndromes[mask]=q;s=shapes[mask]
            atoms[s][q]+=1
            if q:
                image=images[q]
                cancel[s][image.bit_count()][(image^(mask<<start)).bit_count()]+=1
    for s in allowed:
        choices=window_count*len(allowed[s])
        assert sum(atoms[s].values())==choices
        assert sum(sum(row.values()) for row in cancel[s].values())==choices-atoms[s][0]
        for v,count in spectrum.items():
            assert sum(moments[s][v].values())==choices*count
    for mask in range(mask_count):
        for s in allowed:
            assert Counter((mask&x).bit_count() for x in allowed[s])==intersection_counts(shapes[mask],s)
    print('All',len(allowed),'local orbit moments and feedback counts checked; columns',columns_per_bundle,'zero-feedback counts:',
          {s:atoms[s][0] for s in allowed if atoms[s][0]},flush=True)
    return spectrum,allowed,moments,atoms,cancel
