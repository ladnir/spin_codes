"""Sharper independent-map Fourier cap from exhaustively audited supports."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parent.parent))
import bch_model as model
import dense_bch as original
from flint import arb

HERE=Path(__file__).resolve().parent


def load_groups():
    record=model.base.read(HERE/'PACKING_AUDIT.json')
    model.authenticate(record)
    assert record['status']=='EXACT_ALL_DUAL_SUPPORT_BASIS_PACKING'
    assert (record['t'],record['s'])==(128,19)
    assert record['audit']['states']==(1<<19)-1
    return record['audit']['groups']


def bernoulli(engine,theta,lam,groups):
    assert 0<theta<1 and lam>0
    t=len(engine.columns);m=engine.m;n=engine.n
    z=(-lam).exp();g0=1-theta+theta*z;g1=theta+(1-theta)*z
    rho0=min(arb(1),model.up(abs(1-2*theta*z/g0)))
    rho1=min(arb(1),model.up(abs(1-2*(1-theta)*z/g1)))
    entries=[]
    for v in engine.levels:
        cap=arb(1)
        for w,d,e,count in groups:
            lo=max(0,v+w-t,d)
            hi=min(v,w,v-e)
            assert lo<=hi
            cap+=count*max(model.up(rho0**(w-h)*rho1**h) for h in (lo,hi))
        moment=g0**(t-v)*g1**v
        entries.append(model.up(moment*(cap/(2*(m+1))+arb(1)/(2*m))))
    # Zero-source transitions depend only on the unchanged B kernel.
    result=[arb(0)]*(n*n)
    import math
    result[0]=sum((engine.kernel[j]*(theta*z)**j*(1-theta)**(t-j) for j in range(t+1)),arb(0))
    result[1]=sum(((math.comb(t,j)-engine.kernel[j])*(theta*z)**j*(1-theta)**(t-j) for j in range(t+1)),arb(0))
    for i,entry in [(1,max(entries))]+[(j+2,v) for j,v in enumerate(entries)]:
        result[i*n]=entry
        for k,w in enumerate(engine.levels):result[i*n+k+2]=entry*engine.spectrum[w]
    return tuple(map(model.up,result))


class Checker(original.Checker):
    def __init__(self,exponent,ps=None):
        super().__init__(exponent,ps)
        self.groups=load_groups()

    def _fixed_moment(self,index,r):
        old=super()._fixed_moment(index,r)
        lam=(arb(index)/40).exp()
        matrix=bernoulli(self.engine,model.number(r),lam,self.groups)
        moment=model.independent.terminal(matrix,self.engine.n,1<<self.power)
        return min(old,model.up(self.cutoff*lam+moment.log()))
