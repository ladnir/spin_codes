"""Dense gate for cheap updates with one full refresh every R epochs.

No implementation or certificate claim. Regions have 256 epochs, so every
tested period divides the region length and there is no boundary phase drift.
"""
from functools import lru_cache
import json
import math
from pathlib import Path

import numpy as np
import general_occupancies as g

HERE=Path(__file__).resolve().parent


def power(matrix,n):
    size=matrix.shape[0];out=np.full((1,size,size),-np.inf)
    for i in range(size):out[0,i,i]=0
    value=matrix[None]
    while n:
        if n&1:out=g.wm.product(out,value)
        n>>=1
        if n:value=g.wm.product(value,value)
    return out[0]


def main():
    spectrum,kernel,_=g.tv.fixed.load_inner()
    columns=json.loads((g.tv.ROOT/'workstreams/bare_bch_rm2sub/generated/manifest.json').read_text())['t128_s19']['columns']
    caps=g.fiber_caps(128,19,spectrum,kernel);low=g.low_cancellation(columns,19)
    witnesses=json.loads((g.tv.fixed.HERE/'SMALLER_OUTWARD_WITNESSES.json').read_text())
    counts={w:n for w,n in enumerate(g.tv.smaller_outer.spectrum()) if w and n}
    @lru_cache(maxsize=128)
    def epoch(z,full):return g.epochs(spectrum,kernel,columns,caps,low,math.exp(z),None if full else 1,128)
    results=[]
    for row in witnesses['results']:
        dense=row['dense'];groups={period:[] for period in (1,2,4,8,16)}
        for i,box in enumerate(dense['selected_boxes']):
            witness=box['witness'];z=witness['log_surprisal'];lam=math.exp(z)
            ps=np.array(dense['probability_banks'][witness['probability_bank']]['probabilities'])
            p=np.array(witness['proposal']);p/=p.sum();theta=float(ps@p)
            costs=[0.]
            for band,b in zip(dense['bands'],ps[1:]):
                costs.append(math.log(counts[128]) if band==[128] else max(math.log(counts[w])-math.log(math.comb(128,w))-w*math.log(b)-(128-w)*math.log1p(-b) for w in band))
            probs=np.array([math.log(math.comb(128,j))+j*math.log(theta)+(128-j)*math.log1p(-theta) for j in range(129)])
            ordinary=np.logaddexp.reduce(epoch(z,False)+probs[:,None,None],axis=0)
            refresh=np.logaddexp.reduce(epoch(z,True)+probs[:,None,None],axis=0)
            corners=g.typed.vertices(box['lower'],box['upper'],32768)
            common=max(g.typed.point_logs(corners,32768,128,costs,p,0.,row['bad_weight'],lam))+g.typed.lattice_log_count(box['lower'],box['upper'])
            for period in groups:
                block=g.wm.product(power(ordinary,period-1)[None],refresh[None])[0]
                groups[period].append(common+g.terminal(block,32768//period))
            if i%64==0:print(row['distance_target'],i+1,flush=True)
        for period,values in groups.items():
            result=dict(distance_target=row['distance_target'],refresh_period=period,
                        dense_margin_bits=-float(np.logaddexp.reduce(values))/math.log(2),
                        failing_boxes=int(sum(v>-60*math.log(2) for v in values)),box_log_bounds=[float(v) for v in values])
            results.append(result);print(json.dumps({k:v for k,v in result.items() if k!='box_log_bounds'}),flush=True)
    (HERE/'PERIODIC_DENSE_PROBE.json').write_text(json.dumps(dict(status='BINARY64_DENSE_ONLY_NOT_CERTIFICATE',results=results),indent=2)+'\n',encoding='utf-8')


if __name__=='__main__':main()
