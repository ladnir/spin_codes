"""Construct a confined nonzero BCH message for a cache-local permutation.

The general obstruction is rank-nullity, not probabilistic evidence. This
exact instance demonstrates it while retaining the asymmetric inner's maps
and a global state that is not reset between chunks.
"""
import hashlib
import json
from pathlib import Path
import random
import screen_macro as macro
search=macro.search
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]


def dependency(columns):
    pivots={}
    for i,value in enumerate(columns):
        mask=1<<i
        while value:
            bit=value.bit_length()-1
            if bit not in pivots:pivots[bit]=(value,mask);break
            v,m=pivots[bit];value^=v;mask^=m
        if value==0:return mask
    raise AssertionError('more message coordinates than interface rank required')


def main():
    a_path=HERE.parent/'NO_CONSTANT_MAP.json';a=json.loads(a_path.read_text())['columns']
    b_path=HERE.parent/'asymmetric/FEEDBACK_SCREEN.json';bank=json.loads(b_path.read_text())
    b=next(r['b_columns'] for r in bank['candidates'] if r['name']=='greedy3_2')
    outer=search.g.tv.smaller_outer.construction();rows=[int(v,16) for v in outer['generator_rows_hex']]
    assert len(rows)==32
    n=1<<22;local=n//8;epochs=local//128;rng=random.Random(771)
    positions=rng.sample(range(local),128)
    mixers=[]
    for _ in range(epochs):
        u=rng.randrange(1,1<<19);v=rng.randrange(1<<19)
        if (u&v).bit_count()&1:v^=u&-u
        mixers.append((u,v))
    def inputs(word):
        result={}
        for p,position in enumerate(positions):
            if word>>p&1:result[position//128]=result.get(position//128,0)^(1<<(position%128))
        return result
    def terminal(word,measure=False):
        sparse=inputs(word);q=0;weight=0
        for e,(u,v) in enumerate(mixers):
            x=sparse.get(e,0)
            if measure:weight+=(x^search.image(a,q)).bit_count()
            q=search.g.tv.update(q,u,v)^search.inject(b,x)
        return q,weight
    boundary=[terminal(word)[0] for word in rows]
    message=dependency(boundary);word=0
    for i,row in enumerate(rows):
        if message>>i&1:word^=row
    assert message and word
    final,weight=terminal(word,True)
    assert final==0 and 0<weight<=local
    assert weight*200<33*n and weight*100<19*n
    sources=[Path(__file__),Path(macro.__file__),Path(search.__file__),Path(search.g.tv.smaller_outer.__file__),a_path,b_path]
    result=dict(status='EXACT_COUNTEREXAMPLE_TO_CHUNK_LOCAL_PERMUTATIONS_NOT_TO_SPIN',
        output_bits=n,chunk_output_bits=local,chunk_count=8,inner_state_bits=19,
        local_message_dimension=32,boundary_rank=search.rank(boundary),
        message_hex=hex(message),outer_word_hex=hex(word),outer_weight=word.bit_count(),
        chunk_permutation_seed=771,selected_row_positions=positions,
        outgoing_state=final,full_output_weight=weight,relative_weight=weight/n,
        boundary_columns=boundary,
        source_sha256={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in sources})
    (HERE/'LOCALITY_COUNTEREXAMPLE.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in result.items() if k not in ('selected_row_positions','source_sha256','boundary_columns')},indent=2))


if __name__=='__main__':main()
