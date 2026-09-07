"""Audit a complete partial rank cover and the distance-32 F7-kernel subcode.

The complete P dual still has only a distance-30 certificate. Embedded exact
witnesses make replay independent of the discovery beam and legacy folder.
"""
import argparse
import math
from pathlib import Path
import bridge as base
import pdual_rank_refinement as single
import pdual_multi_pivot as multi
from certify_bch_shift_rank import orbit, check_path
from affine_wambach import gf_mul, gf_pow


def forest(legacy):
    def leaf(name):
        return dict(leaf=base.read(base.HERE/'generated/pdual_refinement_v1'/name))
    roots={str(p):dict(retained=base.read(legacy/f'case_{p:03d}.json'))
           for p in base.read(legacy/'summary.json')['case_representatives']}
    roots['7']=dict(split=15,
        zero=dict(split=23,zero=leaf('p7_f15z_f23z.json'),nonzero=leaf('p7_f15z_f23nz_qseed_r26.json')),
        nonzero=leaf('p7_f15nz_qseed_r26.json'))
    roots['15']=dict(split=23,zero=leaf('p15_f23_zero_seeded.json'),
                      nonzero=leaf('p15_f23nz_qseed_w5000.json'))
    return roots


def audit_tree(node,parent,zeros,known,zflags,nzflags,leaves):
    if 'split' in node:
        assert set(node)=={'split','zero','nonzero'}
        e=node['split'];coset=orbit(e)
        assert e==min(coset) and coset.isdisjoint(zeros|known)
        return min(audit_tree(node['zero'],parent,zeros|coset,known,zflags+[e],nzflags,leaves),
                   audit_tree(node['nonzero'],parent,zeros,known|coset,zflags,nzflags+[e],leaves))
    assert set(node) in ({'leaf'},{'retained'})
    if 'leaf' in node:
        record=node['leaf']
        assert record['parent']==parent and record['zero_cosets']==zflags and record['nonzero_cosets']==nzflags
        assert record['zero_indices']==sorted(zeros) and record['known_nonzero_indices']==sorted(known)
        rank=multi.verify(record) if 'initial_pivot' in record['witness'] else single.verify(record)
    else:
        assert not zflags and not nzflags
        record=node['retained']
        assert record['pivot']==parent and record['zero_indices']==sorted(zeros)
        if record.get('proof_type')=='BCH_root_run':
            start,step,count=(record[k] for k in ('start','step','root_count'))
            assert 0<=start<255 and 0<step<255 and math.gcd(step,255)==1 and 0<count<255
            assert all((start+j*step)%255 in zeros for j in range(count))
            rank=count+1
        else:
            state=check_path(zeros,parent,record['steps']);rank=len(state)
            assert record['relative_independent_set']==sorted((e-parent)%255 for e in state)
        assert rank==record['rank']
    leaves.append(dict(parent=parent,zero_cosets=zflags,nonzero_cosets=nzflags,rank=rank))
    return rank


def kernel_extension(data,rank):
    g=int(data['generator_hex'],16);gd=int(data['dual_generator_hex'],16)
    p_rows=[g<<i for i in range(data['dimension'])]
    p_rows=[w|((w.bit_count()%2)<<255) for w in p_rows]
    dual=[gd<<i for i in range(data['dual_dimension'])]+[(1<<256)-1]
    assert len(dual)+len(p_rows)==256
    assert all((w&v).bit_count()%2==0 for w in dual for v in p_rows)
    # Distinct leading degrees give independence of the first 124 rows;
    # the all-one vector has the only nonzero parity coordinate.
    assert len({w.bit_length() for w in dual})==len(dual)
    assert all(w.bit_count()%2==0 for w in dual)
    coords=[gf_pow(2,i) for i in range(255)]+[0]
    pos={x:i for i,x in enumerate(coords)}
    table=[gf_pow(x,7) if x else 0 for x in coords]
    def phi(word):
        value=0
        for i in multi.positions(word):value ^= table[i]
        return value
    images=[phi(w) for w in dual]
    pivots={}
    for v in images:
        while v:
            j=v.bit_length()-1
            if j not in pivots:pivots[j]=v;break
            v ^= pivots[j]
    assert len(pivots)==8
    for a,b in ((1,1),(2,0)):
        permutation=[pos[gf_mul(a,x)^b] for x in coords]
        assert len(set(permutation))==256
        for w in dual:
            moved=sum(1<<permutation[i] for i in multi.positions(w))
            assert phi(moved)==gf_mul(gf_pow(a,7),phi(w))
    assert data['extended_affine_generators_verified'] and rank>=31
    # The two affine generators preserve P dual and ker(phi). Their generated
    # group contains every translation, allowing any zero coordinate to be
    # moved to the parity position before the checked punctured case cover.
    return dict(length=256,dimension=len(dual)-8,dual_basis_rows=len(dual),
        F7_image_rank=8,affine_kernel_invariance_generators_checked=True,
        punctured_kernel_rank_lower=rank,minimum_distance_lower=rank+rank%2,
        conclusion='Every weight-30 word of the extended P dual has nonzero F7.')


def audit(roots):
    data=single.code_data(37)
    zeros=set(data['dual_zero_indices'])
    reps=sorted({min(orbit(e)) for e in set(range(255))-zeros})
    assert set(roots)=={str(p) for p in reps} and reps[0]==7
    leaves=[];ranks={}
    for p in reps:
        assert zeros.isdisjoint(orbit(p))
        ranks[p]=audit_tree(roots[str(p)],p,zeros,orbit(p),[],[],leaves)
        zeros |= orbit(p)
    assert zeros==set(range(255))
    full_rank=min(ranks.values());kernel_rank=min(r for p,r in ranks.items() if p!=7)
    return dict(status='COMPLETE_PARTITION_AUDIT_WITH_PARTIAL_DISTANCE32_PROGRESS',
        full_Pdual_rank_lower=full_rank,full_Pdual_distance_lower=full_rank+full_rank%2,
        full_Pdual_distance32_proved=full_rank>=31,leaf_count=len(leaves),leaves=leaves,
        case_minimum_ranks={str(p):r for p,r in ranks.items()},
        extended_F7_kernel=kernel_extension(data,kernel_rank))


def run(args):
    if args.verify:
        saved=base.read(args.output)
        for name,digest in saved['source_sha256'].items():assert base.sha(base.ROOT/name)==digest
        result=audit(saved['forest'])
        assert result==saved['audit']
    else:
        assert not args.output.exists() and args.legacy is not None
        roots=forest(args.legacy);result=audit(roots)
        sources=[Path(__file__),Path(single.__file__),Path(multi.__file__),
                 base.BCH/'code/certify_bch_shift_rank.py',base.BCH/'code/bch_quotient.py',
                 base.BCH/'code/affine_wambach.py']
        base.write_new(args.output,dict(forest=roots,audit=result,
            source_sha256={p.relative_to(base.ROOT).as_posix():base.sha(p) for p in sources},
            legacy_provenance={p.name:base.sha(p) for p in sorted(args.legacy.glob('*.json'))}))
    print({k:v for k,v in result.items() if k!='leaves'},flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--legacy',type=Path)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--verify',action='store_true')
    run(p.parse_args())
