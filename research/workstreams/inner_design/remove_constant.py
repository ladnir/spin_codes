"""Construct a same-size RM2 subspace with no all-one expansion state.

Extend the retained quadratic space by one dimension, then choose a hyperplane
excluding the constant function. Avoid singleton/triple evaluation syndromes
when choosing that hyperplane. Exact full spectra validate the final map; the
rank-two quadratic screen is only a construction aid, not trusted evidence.
"""
import hashlib
import itertools
import json
from pathlib import Path

import general_occupancies as g

HERE=Path(__file__).resolve().parent


def basis(words):
    pivots={}
    for word in words:
        while word:
            p=word.bit_length()-1
            if p not in pivots:pivots[p]=word;break
            word^=pivots[p]
    return [pivots[p] for p in sorted(pivots,reverse=True)]


def reduce(word,rows):
    for row in rows:
        if word>>(row.bit_length()-1)&1:word^=row
    return word


def main():
    _,old_kernel,sources=g.tv.fixed.load_inner()
    selection_path=g.tv.fixed.BRIDGE/'generated/larger_state_inputs_v1/t128_s19_selection.json'
    selected=json.loads(selection_path.read_text())['selected']
    masks=[int(v,16) for v in selected['quadratic_masks_hex']]
    columns=[int(v,16) for v in selected['B_columns_hex']]
    pairs=list(itertools.combinations(range(7),2))
    evaluations=[sum((((x>>i)&1)*((x>>j)&1))<<k for k,(i,j) in enumerate(pairs)) for x in range(128)]
    reconstructed=[1|(x<<1)|sum(((evaluations[x]&q).bit_count()&1)<<(8+k) for k,q in enumerate(masks)) for x in range(128)]
    assert reconstructed==columns
    rows=basis(masks);assert len(rows)==11
    rank_two=set()
    for u,v in itertools.combinations(range(1,128),2):
        rank_two.add(sum(((((u>>i)&1)*((v>>j)&1))^(((u>>j)&1)*((v>>i)&1)))<<k for k,(i,j) in enumerate(pairs)))
    assert len(rank_two)==2667 and not any(reduce(q,rows)==0 for q in rank_two)
    bad={reduce(q,rows) for q in rank_two}
    pivot_positions={row.bit_length()-1 for row in rows}
    free=[i for i in range(21) if i not in pivot_positions]
    representatives=[sum(((mask>>j)&1)<<p for j,p in enumerate(free)) for mask in range(1<<len(free))]
    extensions=[q for q in representatives if q and q not in bad]
    assert extensions,'retained quadratic space has no screened extension'
    extra=extensions[0]
    parent=[col|(((ev&extra).bit_count()&1)<<19) for col,ev in zip(columns,evaluations)]
    parent_rows=g.tv.fixed.maps.generators(parent,20)
    parent_spectrum=g.tv.fixed.maps.spectrum(parent_rows)
    assert min(w for w in parent_spectrum if w)>=48
    projected=[col>>1 for col in parent]
    assert len(set(projected))==128
    forbidden=bytearray(1<<19)
    for col in projected:forbidden[col]=1
    for a,b,c in itertools.combinations(projected,3):forbidden[a^b^c]=1
    ell=next(i for i,flag in enumerate(forbidden) if not flag)
    new_columns=[col^ell for col in projected]
    new_rows=g.tv.fixed.maps.generators(new_columns,19)
    spectrum=g.tv.fixed.maps.spectrum(new_rows)
    kernel=g.tv.fixed.maps.dual_spectrum(spectrum,128,19)
    assert spectrum.get(128,0)==0 and min(w for w in spectrum if w)>=48
    assert not any((a&b).bit_count()&1 for a in new_rows for b in new_rows)
    assert min(j for j in kernel if j)>=5
    assert all(kernel.get(j,0)<=old_kernel[j] for j in range(0,129,2))
    payload=dict(status='EXACT_MAP_AUDIT_NOT_SPIN_CERTIFICATE',t=128,s=19,
                 construction='Extend quadratics, then delete the constant direction by a shifted evaluation hyperplane',
                 original_quadratic_masks_hex=[hex(q) for q in masks],extra_quadratic_mask_hex=hex(extra),
                 screened_extensions=len(extensions),hyperplane_shift_hex=hex(ell),
                 forbidden_singleton_or_triple_syndromes=sum(forbidden),
                 columns=new_columns,generator_rows_hex=[hex(row) for row in new_rows],
                 spectrum=spectrum,kernel={j:str(n) for j,n in kernel.items()},
                 parent_spectrum=parent_spectrum,
                 minimum_image_weight=min(w for w in spectrum if w),minimum_kernel_weight=min(j for j in kernel if j),
                 all_one_in_image=False,ba_zero=True)
    sources += [Path(__file__),selection_path,Path(g.tv.fixed.maps.__file__)]
    payload['source_sha256']={p.relative_to(g.tv.ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
    (HERE/'NO_CONSTANT_MAP.json').write_text(json.dumps(payload,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in payload.items() if k not in ('columns','generator_rows_hex','kernel','source_sha256')},indent=2))


if __name__=='__main__':main()
