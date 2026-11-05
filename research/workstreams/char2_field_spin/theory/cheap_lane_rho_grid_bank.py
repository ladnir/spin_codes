"""One complete FLOAT bank with hierarchical N8 or N16 rho subdivision.

Reuse the pinned N4 fine-bank orchestration. Compute its completed N4 cap
first, then enter grid2^e*(1+j/N), e=-16..15,j=0..N-1,plus2^16. For N16,
retain the completed N8 cap as well. Each stage intersects its finished
upper cap with the prior stage, so reassignment cannot cause regression.
The raw/powered local table is built once and every fit/marginal uses the
final stage. No extra source/count factor or alpha is introduced.

Changing GRID alone would leave nearest_bucket's old captured default.
Here production assignment passes the new grid explicitly. Frozen N4
tests run on its original globals; independent exact tests cover N8/N16.
All substitutions, argv and private reporting hooks are restored.
One --subdivisions8|16 witness per run, default8; other CLI arguments
are forwarded to the frozen fine-bank wrapper. Its300s/512MiB/20GiB
guards, complete profile coverage and explicit fallbacks remain intact.
FLOAT only, stdout only, no generated files, no certificate claim.
"""

import argparse
import builtins
from collections import defaultdict
from fractions import Fraction as F
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import sys
from types import SimpleNamespace

sys.dont_write_bytecode=True
BANK_PIN="a82044f15e0cfdf20dd85bd92818f49d42c5a24be707ad55eedec077c8f83625"
PATHS=((0,0,1,2),(0,1,1,3),(2,0,3,2),(2,1,3,3),(0,4,1,5),(2,4,3,5))


def load_bank():
    path=Path(__file__).with_name("cheap_lane_rho_fine_bank.py")
    if hashlib.sha256(path.read_bytes()).hexdigest()!=BANK_PIN:
        raise RuntimeError("frozen fine bank mismatch")
    spec=importlib.util.spec_from_file_location("rho_grid_frozen_bank",path)
    module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module


def grid_for(n):
    if n not in (8,16): raise ValueError("choose8 or16 subdivisions")
    return tuple(F(2)**e*F(n+j,n) for e in range(-16,16) for j in range(n))+(F(1<<16),)


def exact_tests(n,nearest,pair,correction,np,guard):
    grid=grid_for(n); checks=0
    assert len(grid)==32*n+1 and tuple(sorted(set(grid)))==grid
    assert max(b/a for a,b in zip(grid,grid[1:]))==F(n+1,n)
    for i,value in enumerate(grid):
        assert nearest(F(1),value,grid)==i; checks+=1
    assert nearest(F(1),F(2),(F(1),F(4),F(16)))==0
    assert nearest(F(1),F(8),(F(1),F(4),F(16)))==1
    assert nearest(F(1),grid[0]/4,grid)==0 and nearest(F(1),grid[-1]*4,grid)==len(grid)-1
    assert nearest(F(0),F(1),grid) is None
    atoms=(((0,0),F(2),(F(1),F(9,8),F(2),F(15,4),F(3),F(7))),
           ((1,0),F(3),(F(3),F(1),F(0),F(5),F(2),F(11))),
           ((0,1),F(1),(F(2),F(3),F(5),F(0),F(13),F(3))))
    assigned={degree:tuple(nearest(row[2*i],row[2*i+1],grid) for i in range(2))
              for degree,w,row in atoms}
    expected=defaultdict(lambda:[F(0)]*6); ideal=defaultdict(lambda:[F(0)]*6)
    for da,wa,a in atoms:
        for db,wb,b in atoms:
            degree=da[0]+db[0],da[1]+db[1]
            for j,(ai,bi,ci,di) in enumerate(PATHS):
                target=min(a[ai]*b[bi],a[ci]*b[di]); ideal[degree][j]+=wa*wb*target
                index=assigned[da][ai//2]
                if index is not None:
                    rho=grid[index]; witness=min(rho*a[ai],a[ci])*min(b[bi]/rho,b[di])
                    assert 0<=witness<=target
                    if grid[0]<=a[ci]/a[ai]<=grid[-1]: assert (n+1)*witness*witness>=n*target*target
                    expected[degree][j]+=wa*wb*witness; checks+=2
    computed=np.full((3,3,6),-math.inf)
    occupied=sorted({i for row in assigned.values() for i in row if i is not None})
    for index in occupied:
        guard(); rho=grid[index]; left=np.zeros((2,2,2)); right=np.zeros((2,2,3))
        for (f,s),weight,row in atoms:
            for i in range(2):
                if assigned[f,s][i]==index: left[f,s,i]=float(weight*min(rho*row[2*i],row[2*i+1]))
            right[f,s]=list(map(float,(weight*min(row[0]/rho,row[2]),
                weight*min(row[1]/rho,row[3]),weight*min(row[4]/rho,row[5]))))
        lower=correction.cross_convolution(pair.as_logs(left,np),pair.as_logs(right,np),np,guard)
        np.logaddexp(computed,lower,out=computed)
    for degree,row in expected.items():
        for j,value in enumerate(row):
            assert value<=ideal[degree][j]
            got=math.exp(float(computed[degree][j])) if math.isfinite(computed[degree][j]) else 0.0
            assert abs(got-float(value))<=3e-12*max(1,float(value)); checks+=2
    return {"subdivisions":n,"grid_size":len(grid),"exact_checks":checks,
        "explicit_grid_argument":True,"noncommuting_six_channels":True,
        "inside_grid_squared_loss_control":True,"zero_clip_and_rational_tie_controls":True}


def main():
    parser=argparse.ArgumentParser(description=__doc__,add_help=False)
    parser.add_argument("--subdivisions",type=int,choices=(8,16),default=8)
    args,remaining=parser.parse_known_args()
    levels=(8,) if args.subdivisions==8 else (8,16)
    bank=load_bank(); restorations=[]
    def replace(obj,name,value):
        restorations.append((obj,name,getattr(obj,name))); setattr(obj,name,value)
    original_load=bank.load
    def intercepted_load(name,pin,label):
        module=original_load(name,pin,label)
        if name=="cheap_lane_rho_fine_buckets.py":
            original_builder,original_tests,nearest=module.fine_bucket_pair_logs,module.exact_tests,module.nearest_bucket
            original_grid=module.GRID
            def builder(*values,**kwargs):
                # N4 must finish BEFORE changing module globals. Later callbacks
                # return retained coefficients; they do not reconstruct N4.
                tables,report=original_builder(*values,**kwargs)
                for n in levels:
                    previous_tables,previous_report=tables,report
                    def retained(*unused,**ignored):
                        kept=dict(previous_tables)
                        kept["dense33_bucket_sum"]=previous_tables["fine129_bucket_sum"]
                        return kept,previous_report
                    grid=grid_for(n)
                    def explicit_nearest(a,c,chosen=grid): return nearest(a,c,chosen)
                    try:
                        module.GRID=grid; module.nearest_bucket=explicit_nearest
                        stage_kwargs=dict(kwargs); stage_kwargs["coarse_builder"]=retained
                        tables,report=original_builder(*values,**stage_kwargs)
                    finally:
                        module.GRID=original_grid; module.nearest_bucket=nearest
                    report["grid"]=f"2^e*(1+j/{n}),e=-16..15,j=0..{n-1},plus2^16"
                    report["subdivisions"]=n
                    report["retained_levels"]=[4]+[x for x in levels if x<=n]
                    report["legacy_fine129_key_contains_current_finished_grid_cap"]=True
                return tables,report
            def tests(pair,correction,np,guard):
                return {"frozen_N4":original_tests(pair,correction,np,guard),
                    "higher_grids":[exact_tests(n,nearest,pair,correction,np,guard) for n in levels]}
            replace(module,"fine_bucket_pair_logs",builder); replace(module,"exact_tests",tests)
        return module
    def emit(*values,**kwargs):
        if len(values)==1 and isinstance(values[0],str):
            row=json.loads(values[0])
            if isinstance(row,dict) and "result" in row and "pair_bound" in row:
                row["N4_wrapper_pair_bound_label_overridden"]=row["pair_bound"]
                row["pair_bound"]=f"hierarchical disjoint-left rational rho buckets N4 through N{args.subdivisions}; all prior finished upper caps retained"
                row["parameters"]["rho_subdivisions"]=args.subdivisions
                row["parameters"]["retained_rho_levels"]=[4,*levels]
                row["grid_wrapper_sha256"]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
                row["frozen_N4_bank_pin"]=BANK_PIN
                values=(json.dumps(row,indent=2),)
        builtins.print(*values,**kwargs)
    original_argv=sys.argv[:]
    try:
        replace(bank,"load",intercepted_load)
        # Replace only this private module binding, never builtins.print.
        replace(bank,"builtins",SimpleNamespace(print=emit))
        sys.argv=[sys.argv[0],*remaining]
        bank.main()
    finally:
        sys.argv=original_argv
        for obj,name,old in reversed(restorations): setattr(obj,name,old)


if __name__=="__main__":
    main()
