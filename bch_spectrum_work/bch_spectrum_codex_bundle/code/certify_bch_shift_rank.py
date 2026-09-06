"""Independent exact audit of complete first-nonzero Fourier rank certificates."""
import json
import math
import sys
from pathlib import Path
from bch_quotient import generator_polynomial,binary_poly_divmod
from affine_wambach import gf_mul,gf_pow
from run_higher_endpoint_preflight import write_new,sha
ROOT=Path(__file__).resolve().parents[1]


def orbit(e,n=255):
    result=[]; x=e%n
    while x not in result:
        result.append(x); x=x*2%n
    assert x==e%n
    return set(result)


def evaluate(poly,x):
    value=0
    for j in range(poly.bit_length()-1,-1,-1):
        value=gf_mul(value,x)^((poly>>j)&1)
    return value


def in_span(word,basis):
    pivots={}
    for v in basis:
        while v:
            j=v.bit_length()-1
            if j not in pivots:
                pivots[j]=v; break
            v^=pivots[j]
    while word:
        j=word.bit_length()-1
        if j not in pivots:
            return False
        word^=pivots[j]
    return True


def check_path(zeros,pivot,steps,n=255):
    assert pivot not in zeros and 0<=pivot<n
    state={pivot}
    for k,shift in steps:
        assert isinstance(k,int) and 0<=k<8 and isinstance(shift,int) and 0<=shift<n
        actual_shift=(shift+(1-(1<<k))*pivot)%n
        moved={((1<<k)*e+actual_shift)%n for e in state}
        assert len(moved)==len(state) and moved<=zeros
        assert pivot not in moved
        state=moved|{pivot}
    return state


def code_data(distance):
    g=generator_polynomial(distance)
    k=255-(g.bit_length()-1)
    quotient,remainder=binary_poly_divmod((1<<255)|1,g)
    assert remainder==0
    degree=quotient.bit_length()-1
    gd=sum(((quotient>>i)&1)<<(degree-i) for i in range(degree+1))
    basis=[g<<i for i in range(k)]
    dual_basis=[gd<<i for i in range(255-k)]
    assert all((x&y).bit_count()%2==0 for x in basis for y in dual_basis)
    # Distinct leading degrees give the asserted dimensions.
    roots={e for e in range(255) if evaluate(gd,gf_pow(2,e))==0}
    original_roots={e for e in range(255) if evaluate(g,gf_pow(2,e))==0}
    assert len(roots)==degree==k and len(original_roots)==255-k
    assert roots=={-e%255 for e in set(range(255))-original_roots}
    assert 0 in roots and all(w.bit_count()%2==0 for w in dual_basis)
    assert all(orbit(e)<=roots for e in roots)
    extended=[w|((w.bit_count()%2)<<255) for w in basis]
    coords=[gf_pow(2,i) for i in range(255)]+[0]
    pos={x:i for i,x in enumerate(coords)}
    assert len(pos)==256
    for a,b in ((1,1),(2,0)):
        permutation=[pos[gf_mul(a,x)^b] for x in coords]
        assert len(set(permutation))==256
        for w in extended:
            moved=sum(1<<permutation[i] for i in range(256) if w>>i&1)
            assert in_span(moved,extended)
    return dict(dimension=k,generator_hex=hex(g),dual_generator_hex=hex(gd),
        dual_dimension=255-k,dual_zero_indices=sorted(roots),
        dual_evenness_verified=True,extended_affine_generators_verified=True)


def audit_folder(distance,folder_name):
    assert folder_name.isidentifier()
    folder=ROOT/'generated'/folder_name
    data=code_data(distance)
    roots=set(data['dual_zero_indices'])
    reps=sorted({min(orbit(e)) for e in set(range(255))-roots})
    summary=json.loads((folder/'summary.json').read_text())
    assert summary['original_BCH_designed_distance']==distance and summary['case_representatives']==reps
    zeros=set(roots); cases=[]; paths=[folder/'summary.json']
    for pivot in reps:
        path=folder/f'case_{pivot:03d}.json'
        record=json.loads(path.read_text())
        assert record['pivot']==pivot and record['zero_indices']==sorted(zeros)
        assert orbit(pivot).isdisjoint(zeros)
        if record.get('proof_type')=='BCH_root_run':
            start,step,count=(record[k] for k in ('start','step','root_count'))
            assert 0<=start<255 and 0<step<255 and math.gcd(step,255)==1
            assert 0<count<255 and all((start+j*step)%255 in zeros for j in range(count))
            assert record['rank']==count+1
            cases.append(dict(pivot=pivot,rank=count+1,proof_type='BCH_root_run',
                start=start,step=step,root_count=count))
        else:
            independent=check_path(zeros,pivot,record['steps'])
            assert len(independent)==record['rank']
            assert sorted((e-pivot)%255 for e in independent)==record['relative_independent_set']
            cases.append(dict(pivot=pivot,rank=len(independent),steps_checked=len(record['steps']),proof_type='independent_set'))
        zeros|=orbit(pivot); paths.append(path)
    assert zeros==set(range(255))
    rank=min(c['rank'] for c in cases)
    even=rank+rank%2
    data.update(original_BCH_designed_distance=distance,folder_name=folder_name,
        cases=cases,minimum_certified_rank=rank,all_nonzero_Fourier_cases_covered=True,
        punctured_dual_distance_lower=even,extended_dual_distance_lower=even,
        proof_inputs_from_literature_table=False,
        witness_sha256={str(p.relative_to(ROOT)):sha(p) for p in paths})
    return data


def toy_check():
    # Exhaust every nonzero binary word of length five. This subgroup lies in F256.
    n=5; root=gf_pow(2,51); tested=0
    for word in range(1,1<<n):
        support=[j for j in range(n) if word>>j&1]
        fourier=[0]*n
        for e in range(n):
            for j in support:
                fourier[e]^=gf_pow(root,e*j)
        zeros={e for e,f in enumerate(fourier) if f==0}
        assert len(zeros)<n
        pivot=next(e for e in range(n) if e not in zeros)
        states={(pivot,):[]}
        for _ in range(n):
            updated={}
            for state,path in states.items():
                assert len(state)<=len(support)
                for k in range(4):
                    for shift in range(n):
                        actual=(shift+(1-(1<<k))*pivot)%n
                        moved={((1<<k)*e+actual)%n for e in state}
                        if moved<=zeros:
                            new=tuple(sorted(moved|{pivot}))
                            steps=path+[(k,shift)]
                            assert check_path(zeros,pivot,steps,n)==set(new)
                            assert len(new)<=len(support)
                            updated[new]=steps; tested+=1
            if not updated:
                break
            states=updated
    return dict(all_31_nonzero_length5_words_tested=True,valid_transitions_checked=tested)


def build():
    q=audit_folder(39,'shift_rank_q24_complete')
    p=audit_folder(37,'shift_rank_p26_complete')
    assert q['extended_dual_distance_lower']>=24 and p['extended_dual_distance_lower']>=26
    toy=toy_check()
    return dict(classification='Independent exact complete Fourier-case rank certificate for both BCH dual distances',
        P=p,Q=q,toy_regression=toy,published_SchaubPlus_table_used=False,
        Q_and_cosets_OA_strength_at_least=23,P_OA_strength_at_least=25,
        original_M22_target_closed=False,source_sha256={str(path.relative_to(ROOT)):sha(path) for path in
            (Path(__file__),ROOT/'code/bch_quotient.py',ROOT/'code/affine_wambach.py')})


if __name__=='__main__':
    value=build()
    path=ROOT/'generated/bch256_shift_rank_dual_distances.json'
    if '--verify' in sys.argv:
        assert value==json.loads(path.read_text())
    else:
        write_new(path,value)
    print(json.dumps(dict(classification=value['classification'],
        Pdual_distance_lower=value['P']['extended_dual_distance_lower'],
        Qdual_distance_lower=value['Q']['extended_dual_distance_lower'],
        P_cases=len(value['P']['cases']),Q_cases=len(value['Q']['cases']),
        toy_regression=value['toy_regression'],published_table_used=False),indent=2))
