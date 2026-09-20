"""Read-only search for stronger BCH root runs under cyclic multipliers."""
import json
import math
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]


def build():
    source=json.loads((ROOT/'generated/bch256_hull_structure.json').read_text())
    units=[m for m in range(1,255) if math.gcd(m,255)==1]
    result={}
    for label in ('HP','HQ','HPdual','HQdual'):
        roots=set(source['codes'][label]['roots'])
        best=(0,0,0)
        for m in units:
            transformed={(m*r)%255 for r in roots}
            for b in transformed:
                if (b-1)%255 in transformed:
                    continue
                k=0
                while (b+k)%255 in transformed:
                    k+=1
                best=max(best,(k,m,b))
        length,m,b=best
        bound=length+1
        result[label]=dict(consecutive_root_count=length,multiplier=m,start=b,
                           extended_even_distance_bound=bound+bound%2)
    return result


if __name__=='__main__':
    print(json.dumps(build(),indent=2))
