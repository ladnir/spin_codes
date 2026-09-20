"""Remove inequalities redundant modulo the equality space (safe weakening).

Modular membership is a selection heuristic only. Deleting constraints cannot
invalidate an eventual upper certificate, even if dependence fails over Q.
"""
import json
import sys
import textwrap
from pathlib import Path
sys.dont_write_bytecode=True
from flint import nmod_mat
from prepare_bch_split38_probe import export
from run_higher_endpoint_preflight import sha,write_new
ROOT=Path(__file__).resolve().parents[1]
FOLDER=ROOT/'generated/split38_clean_probe'


def main():
    FOLDER.mkdir(exist_ok=False)
    source=ROOT/'generated/split38_balanced_probe/model.json'
    model=json.loads(source.read_text())
    variables=model['variables']
    index={n:i for i,n in enumerate(variables)}
    equations=[r for r in model['constraints'] if r['sense']=='eq']
    inequalities=[r for r in model['constraints'] if r['sense']!='eq']
    prime=2147483647
    def dense(rows):
        output=[]
        for row in rows:
            values=[0]*(len(variables)+1)
            for n,v in row['coeffs'].items():
                values[index[n]]=int(v)%prime
            values[-1]=int(row['rhs'])%prime
            output.append(values)
        return output
    e,rank=nmod_mat(dense(equations),prime).rref()
    assert rank==len(equations)
    pivots=[next(j for j in range(len(variables)+1) if e[i,j]) for i in range(rank)]
    raw=dense(inequalities)
    im=nmod_mat(raw,prime)
    factors=nmod_mat([[row[j] for j in pivots] for row in raw],prime)
    residual=im-factors*e
    removed={}
    for i,row in enumerate(inequalities):
        if all(residual[i,j]==0 for j in range(len(variables)+1)):
            removed[row['name']]='augmented row lies in equality span modulo prime'
        elif len(row['coeffs'])==1 and int(row['rhs'])==0:
            value=int(next(iter(row['coeffs'].values())))
            if (row['sense']=='ge' and value>0) or (row['sense']=='le' and value<0):
                removed[row['name']]='duplicate nonnegative variable bound'
    model['constraints']=[r for r in model['constraints'] if r['name'] not in removed]
    model['classification']='Safe weakening of the split relaxation with duplicate/selected equality-dependent inequalities removed'
    model['removed_inequalities']=removed
    model['source_sha256']={str(p.relative_to(ROOT)):sha(p) for p in (Path(__file__),source)}
    write_new(FOLDER/'model.json',model)
    write_new(FOLDER/'all_rows_after_substitution.json',json.loads((ROOT/'generated/split38_local_probe/all_rows_after_substitution.json').read_text()))
    export(model,FOLDER/'h38_unwrapped.lp')
    raw=(FOLDER/'h38_unwrapped.lp').read_text()
    wrapped='\n'.join('\n'.join(textwrap.wrap(line,width=2000,break_long_words=False,break_on_hyphens=False,
                                             subsequent_indent=' ')) for line in raw.splitlines())+'\n'
    assert raw.split()==wrapped.split()
    with (FOLDER/'h38.lp').open('x',encoding='ascii') as stream:
        stream.write(wrapped)
    print(json.dumps(dict(variables=len(variables),equalities=len(equations),
                          original_inequalities=len(inequalities),removed_inequalities=len(removed),
                          remaining_rows=len(model['constraints'])),indent=2))


if __name__=='__main__':
    main()
