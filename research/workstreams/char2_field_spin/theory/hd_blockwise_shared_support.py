"""Exact support constraints for the shared FFT8 / four-state inner.

Ranks are algebraic facts about the stated joint variable space. No random
amplitude assumption or whole-code probability claim is made.
"""

from functools import lru_cache
from itertools import combinations

from hd_blockwise_maps import matmul,matvec
from hd_blockwise_shared_maps import (matrices,identity,rank,power,
                                      multiplication_table)


def rref(matrix,columns=None):
    if not matrix:
        return (),()
    work=[list(row) for row in matrix]
    table=multiplication_table()
    end=len(work[0]) if columns is None else columns
    pivots=[]
    for column in range(end):
        pivot=next((i for i in range(len(pivots),len(work)) if work[i][column]),None)
        if pivot is None:
            continue
        row=len(pivots)
        work[row],work[pivot]=work[pivot],work[row]
        scale=power(work[row][column],254)
        work[row]=[table[scale][x] for x in work[row]]
        for i in range(len(work)):
            scale=work[i][column]
            if i!=row and scale:
                work[i]=[x^table[scale][y] for x,y in zip(work[i],work[row])]
        pivots.append(column)
        if len(pivots)==len(work):
            break
    return tuple(tuple(row) for row in work),tuple(pivots)


def kernel(matrix,ncolumns):
    """Columns form a basis of the right kernel; supports empty matrices."""
    if not matrix:
        return identity(ncolumns)
    reduced,pivots=rref(matrix)
    free=[j for j in range(ncolumns) if j not in pivots]
    vectors=[]
    for j in free:
        vector=[0]*ncolumns
        vector[j]=1
        for i,pivot in enumerate(pivots):
            vector[pivot]=reduced[i][j]
        vectors.append(vector)
    return tuple(tuple(vector[i] for vector in vectors) for i in range(ncolumns))


def _product(left,right):
    if not left:
        return ()
    if not right or not right[0]:
        return tuple(() for _ in left)
    return matmul(left,right)


def _support(values,limit=64):
    result=tuple(sorted(values))
    if len(set(result))!=len(result) or any(not 0<=v<limit for v in result):
        raise ValueError('invalid coordinate set')
    return result


def output_constraints(input_support,zero_outputs,translation_offset=8):
    support=_support(input_support)
    zeros=_support(zero_outputs)
    f,a,_,_,_=matrices(translation_offset)
    return tuple(tuple(f[i][j] for j in support)+a[i] for i in zeros)


def state_residual(tile,input_support,zero_outputs,translation_offset=8):
    """Local coordinates; eliminate input columns, retain four state columns."""
    if not 0<=tile<8:
        raise ValueError('invalid tile')
    support=_support(input_support,8)
    zeros=_support(zero_outputs,8)
    block=output_constraints(tuple(8*tile+j for j in support),
                             tuple(8*tile+j for j in zeros),translation_offset)
    if not block:
        return dict(input_pivots=0,residual=(),rank=0)
    reduced,pivots=rref(block,len(support))
    residual=tuple(tuple(row[len(support):]) for row in reduced[len(pivots):])
    return dict(input_pivots=len(pivots),residual=residual,rank=rank(residual))


def stacked_output_rank(input_support,zero_outputs,translation_offset=8):
    """Exact output rank = local input pivots + joint state residual rank."""
    support=_support(input_support)
    zeros=_support(zero_outputs)
    residual=[]
    pivot_count=0
    local_ranks=[]
    for tile in range(8):
        local=state_residual(tile,[x%8 for x in support if x//8==tile],
                             [x%8 for x in zeros if x//8==tile],translation_offset)
        pivot_count+=local['input_pivots']
        residual.extend(local['residual'])
        local_ranks.append(local['rank'])
    state_rank=rank(residual)
    return dict(input_pivots=pivot_count,state_rank=state_rank,
                local_state_ranks=local_ranks,output_rank=pivot_count+state_rank)


def rank_profile(input_support,zero_outputs,translation_offset=8,
                 source_embedding=None,mixer=None,require_live_state=True):
    """Impose output zeros on exact joint source forms, then inspect images.

    Variables are ordered as supported X coordinates followed by four Q
    coordinates. Optional source_embedding has that many rows and linearly
    independent columns. It retains any prescribed input/state dependencies.
    mixer is a fixed 4x4 map; its randomness is not analyzed by this helper.
    """
    support=_support(input_support)
    zeros=_support(zero_outputs)
    w=len(support)
    _,_,c,_,_=matrices(translation_offset)
    embedding=identity(w+4) if source_embedding is None else tuple(tuple(row) for row in source_embedding)
    if len(embedding)!=w+4 or len({len(row) for row in embedding})!=1:
        raise ValueError('invalid source embedding')
    dimension=len(embedding[0])
    if rank(embedding)!=dimension:
        raise ValueError('source embedding columns must be independent')
    mix=identity(4) if mixer is None else tuple(tuple(row) for row in mixer)
    if len(mix)!=4 or any(len(row)!=4 for row in mix):
        raise ValueError('invalid fixed mixer')
    outputs=_product(output_constraints(support,zeros,translation_offset),embedding)
    surviving=kernel(outputs,dimension)
    ambient_survivors=_product(embedding,surviving)
    feedback_map=tuple(tuple(c[i][j] for j in support)+(0,)*4 for i in range(4))
    successor_map=tuple(tuple(c[i][j] for j in support)+mix[i] for i in range(4))
    feedback=_product(feedback_map,ambient_survivors)
    successor=_product(successor_map,ambient_survivors)
    state=ambient_survivors[w:]

    def live_possible(vectors):
        inputs_live=all(any(row) for row in vectors[:w])
        state_live=any(x for row in vectors[w:] for x in row)
        return inputs_live and (state_live or not require_live_state)

    returns=_product(successor_map,embedding)
    return_kernel=kernel(outputs+returns,dimension)
    return_survivors=_product(embedding,return_kernel)
    return dict(input_weight=w,zero_outputs=len(zeros),source_dimension=dimension,
                output_rank=rank(outputs),output_nullity=dimension-rank(outputs),
                state_image_rank=rank(state),feedback_image_rank=rank(feedback),
                successor_image_rank=rank(successor),
                joint_output_feedback_rank=rank(outputs)+rank(feedback),
                joint_output_return_rank=rank(outputs)+rank(successor),
                full_input_support_possible=live_possible(ambient_survivors),
                live_input_and_return_possible=live_possible(return_survivors),
                support_feasibility_field_size=256,
                whole_certificate=False)


def paired_masks(support,pair_start=0):
    """Complementary swapped-output masks of the beta8 defect."""
    support=_support(support)
    if pair_start%16 or any(not pair_start<=x<pair_start+16 for x in support):
        raise ValueError('support must lie in one aligned pair')
    zero_outputs=tuple(i for i in range(pair_start,pair_start+16) if (i^8) not in support)
    return support,zero_outputs


def tile_generator(tile=0,translation_offset=8):
    """Generator rows are the 16 observed coordinates, columns are X8,Q4."""
    output=output_constraints(tuple(range(8*tile,8*tile+8)),
                              tuple(range(8*tile,8*tile+8)),translation_offset)
    return tuple(row+(0,)*4 for row in identity(8))+output


def defect_receipts():
    support,zeros=paired_masks((0,1,2,8,9,10))
    return {str(beta):dict(rank_profile(support,zeros,beta),
                            elimination=stacked_output_rank(support,zeros,beta))
            for beta in (8,64)}


@lru_cache(None)
def residual_plane_catalog(tile,translation_offset):
    """All w=3,z=5 masks on one tile; retain one mask per residual plane."""
    catalog={}
    for support in combinations(range(8),3):
        for zeros in combinations(range(8),5):
            local=state_residual(tile,support,zeros,translation_offset)
            reduced,pivots=rref(local['residual'])
            key=reduced[:len(pivots)]
            assert len(key)==2
            catalog.setdefault(key,(support,zeros))
    return tuple(catalog.items())


def plane_comparison():
    """Bounded exhaustive comparison of the first two tiles, not a law."""
    inventories={}
    for beta in (8,64):
        left=dict(residual_plane_catalog(0,beta))
        right=dict(residual_plane_catalog(1,beta))
        inventories[str(beta)]=dict(unique_planes=[len(left),len(right)],
                                   identical_planes=len(left.keys()&right.keys()))
    table=multiplication_table()
    minors=tuple(combinations(range(4),2))
    plucker=[]
    for tile in (0,1):
        rows=[]
        for basis,masks in residual_plane_catalog(tile,64):
            coordinates=tuple(table[basis[0][i]][basis[1][j]]
                              ^table[basis[0][j]][basis[1][i]] for i,j in minors)
            rows.append((coordinates,basis,masks))
        plucker.append(rows)
    deficient=0
    first=None
    for p,left,s in plucker[0]:
        products=tuple(table[x] for x in p)
        for q,right,t in plucker[1]:
            determinant=(products[0][q[5]]^products[1][q[4]]^products[2][q[3]]
                         ^products[3][q[2]]^products[4][q[1]]^products[5][q[0]])
            if determinant==0:
                deficient+=1
                if first is None:
                    first=dict(first_tile_masks=s,second_tile_masks=t,
                               stacked_rank=rank(left+right))
    inventories['64'].update(plane_pairs=len(plucker[0])*len(plucker[1]),
                             rank_three_pairs=deficient,first_rank_three_example=first)
    return dict(inventories=inventories,tiles=[0,1],input_weight_per_tile=3,
                zero_outputs_per_tile=5,uniform_probability_claim=False,
                whole_certificate=False)


if __name__=='__main__':
    import json
    print(json.dumps(dict(defects=defect_receipts(),planes=plane_comparison()),indent=2))
