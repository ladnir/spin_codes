"""Change only exact column scaling in the saved local split relaxation.

Square-root binomial scales balance the Krawtchouk transform more evenly.
Every scale is an integer; the feasible physical spectra are unchanged.
"""
import json
import math
import sys
import textwrap
from pathlib import Path
sys.dont_write_bytecode=True
from prepare_bch_split38_probe import export
from run_higher_endpoint_preflight import sha,write_new
ROOT=Path(__file__).resolve().parents[1]
FOLDER=ROOT/'generated/split38_balanced_probe'


def main():
    FOLDER.mkdir(exist_ok=False)
    source=ROOT/'generated/split38_local_probe/model.json'
    model=json.loads(source.read_text())
    for name in model['variables']:
        if name.startswith('b_'):
            _,i,j=name.split('_')
            volume=math.comb(38,int(i))*math.comb(218,int(j))
        else:
            volume=math.comb(256,int(name.split('_')[1]))
        model['scales'][name]=max(1,math.isqrt(volume))
    model['objective']['physical_multiplier']=model['scales']['h_38']
    model['classification']='Same local split relaxation with exact integer square-root-binomial column scales'
    model['source_sha256']={str(p.relative_to(ROOT)):sha(p) for p in (Path(__file__),source)}
    write_new(FOLDER/'model.json',model)
    all_rows=json.loads((ROOT/'generated/split38_local_probe/all_rows_after_substitution.json').read_text())
    write_new(FOLDER/'all_rows_after_substitution.json',all_rows)
    export(model,FOLDER/'h38_unwrapped.lp')
    raw=(FOLDER/'h38_unwrapped.lp').read_text()
    wrapped='\n'.join('\n'.join(textwrap.wrap(line,width=2000,break_long_words=False,break_on_hyphens=False,
                                             subsequent_indent=' ')) for line in raw.splitlines())+'\n'
    assert raw.split()==wrapped.split()
    with (FOLDER/'h38.lp').open('x',encoding='ascii') as stream:
        stream.write(wrapped)
    print('Prepared algebraically identical split LP with new positive integer scales')


if __name__=='__main__':
    main()
