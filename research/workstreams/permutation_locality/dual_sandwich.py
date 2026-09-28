"""Exact repaired LP witnesses for dual shells of the fixed BCH sandwich.

HiGHS proposes multipliers. Exact residuals are covered using authenticated
variable upper bounds. Empirical orbit-search lower constraints are omitted.
Only the checked rational upper, never the numerical optimum, is retained.
"""
import argparse
import json
from fractions import Fraction as Q
from math import comb, log2
from pathlib import Path
import sys
import numpy as np
from scipy.optimize import linprog
from bch_joint_support import authenticated_caps
from dual_moments import dual_shell_caps
from dual_shortening import verify_bch_premise
from shortened_bound import krawtchouk


def raw_model():
    repo=Path(__file__).resolve().parents[3]
    bundle=repo/'research/bch_spectrum_work/bch_spectrum_codex_bundle'
    sys.path.insert(0,str(bundle/'code'))
    from audit_bch_closure_envelope import audit
    from export_lp import independent_oa_constraints
    # This re-derives every retained original row and the quotient algebra.
    checked=audit()
    assert checked==json.loads((bundle/'generated/bch256_closure_deterministic_envelope.json').read_text())
    verify_bch_premise()
    caps=authenticated_caps()
    original=json.loads((bundle/'generated/coupled_lp_exact.json').read_text())
    raw=independent_oa_constraints(original)
    raw=[r for r in raw if not r['name'].startswith(('Wambach_','orbit_search_','q_nonneg_','h_nonneg_'))]
    weights=list(range(0,129,2))
    table={w:krawtchouk(256,w) for w in weights}
    for prefix in ('q','h'):
        for degree in range(16,30,2):
            raw.append(dict(name=f'{prefix}_OA29_{degree}',sense='eq',rhs=0,
                            coeffs={f'{prefix}_{w}':table[w][degree]*(2 if w<128 else 1) for w in weights}))
    for w in range(38,129,2):
        raw.append(dict(name=f'authenticated_C_cap_{w}',sense='le',rhs=caps[w],coeffs={f'q_{w}':1,f'h_{w}':31}))
    fixed={f'{p}_{w}':int(p=='q' and w==0) for p in ('q','h') for w in weights
           if w<(40 if p=='q' else 38)}
    return original['metadata']['variables'],raw,fixed,caps,table


def model():
    original_variables,raw,fixed,caps,table=raw_model()
    variables=[x for x in original_variables if x not in fixed]
    scales=[max(1,comb(256,int(x.split('_')[1]))//(1<<132)) for x in variables]
    upper=[Q(caps[int(x.split('_')[1])],s*(31 if x.startswith('h_') else 1)) for x,s in zip(variables,scales)]
    inequalities=[];equalities=[]
    for raw_row in raw:
        rhs=Q(int(raw_row['rhs']))-sum(int(raw_row['coeffs'].get(x,0))*v for x,v in fixed.items())
        row=[Q(int(raw_row['coeffs'].get(x,0))*s) for x,s in zip(variables,scales)]
        sense=raw_row['sense']
        if sense=='ge':row=[-x for x in row];rhs=-rhs;sense='le'
        if not any(row):
            assert rhs==0 if sense=='eq' else rhs>=0
            continue
        scale=max([abs(rhs)]+[abs(x) for x in row])
        normalized=([x/scale for x in row],rhs/scale,raw_row['name'])
        (equalities if sense=='eq' else inequalities).append(normalized)
    print('Replayed sandwich model:',len(variables),'variables,',len(inequalities),'inequalities,',len(equalities),'equalities; no empirical lower constraints',flush=True)
    return variables,scales,upper,inequalities,equalities,fixed,table


def repaired_upper(objective,constant,upper,inequalities,equalities,result):
    """All arguments except proposed floating multipliers are exact rationals."""
    prices=[max(Q(0),-Q(float(x))) for x in result.ineqlin.marginals]
    eqprices=[-Q(float(x)) for x in result.eqlin.marginals]
    boundprices=[max(Q(0),-Q(float(x))) for x in result.upper.marginals]
    return price_upper(objective,constant,upper,inequalities,equalities,prices,eqprices,boundprices)


def price_upper(objective,constant,upper,inequalities,equalities,prices,eqprices,boundprices=None):
    boundprices=boundprices if boundprices is not None else [Q(0)]*len(upper)
    assert all(x>=0 for x in prices+boundprices)
    columns=boundprices[:]
    bound=constant+sum((a*b for a,b in zip(boundprices,upper)),Q(0))
    for price,(row,rhs,_) in zip(prices+eqprices,inequalities+equalities):
        if not price:continue
        bound+=price*rhs
        for i,x in enumerate(row):columns[i]+=price*x
    residual=[max(Q(0),a-b) for a,b in zip(objective,columns)]
    correction=sum((a*b for a,b in zip(residual,upper)),Q(0))
    assert all(b+r>=a for a,b,r in zip(objective,columns,residual))
    assert all(x>=0 for x in prices+boundprices+residual)
    return bound+correction,correction


def reconstructed_prices(rows,objective,fit,ni,limit):
    """Propose rational prices by solving selected active columns exactly.

    This need not preserve feasibility: price_upper checks signs and repairs
    every column afterwards. Missing rank simply rejects the proposal.
    """
    from flint import fmpq_mat,fmpq
    nr=len(rows);prices=[Q(float(x)) for x in fit.x[:nr]]
    unknown=[]
    for j,x in enumerate(prices):
        if abs(float(x))<1e-10:prices[j]=Q(0)
        elif abs(abs(float(x))-limit)<max(1e-8,limit*1e-9):prices[j]=Q(limit if x>0 else -limit)
        else:unknown.append(j)
    active=[i for i,x in enumerate(fit.x[nr:]) if x<1e-8]
    if not unknown or len(active)<len(unknown):return None
    matrix=fmpq_mat([[fmpq(rows[j][0][i].numerator,rows[j][0][i].denominator) for j in unknown] for i in active])
    reduced,rank=matrix.transpose().rref()
    if rank!=len(unknown):return None
    selected=[active[next(i for i in range(len(active)) if reduced[j,i])] for j in range(rank)]
    square=fmpq_mat([[fmpq(rows[j][0][i].numerator,rows[j][0][i].denominator) for j in unknown] for i in selected])
    fixed=[j for j in range(nr) if j not in unknown]
    rhs=[objective[i]-sum((prices[j]*rows[j][0][i] for j in fixed),Q(0)) for i in selected]
    answer=square.solve(fmpq_mat([[fmpq(x.numerator,x.denominator)] for x in rhs]))
    for i,j in enumerate(unknown):prices[j]=Q(int(answer[i,0].numerator),int(answer[i,0].denominator))
    prices[:ni]=[max(Q(0),x) for x in prices[:ni]]
    print('Reconstructed',rank,'prices from exact active-column equations',flush=True)
    return prices


def exact_solver(objective,constant,upper,inequalities,equalities,variables,tight=False):
    import subprocess,tempfile
    runtime=Path('C:/Users/peter/.codex/worktrees/ba80/permute_conv/bch_spectrum_work/bch_spectrum_codex_bundle/generated/split38_soplex_tools/runtime')
    assert (runtime/'usr/bin/soplex').is_file()
    linux=lambda p:'/mnt/'+str(p.resolve()).replace('\\','/')[0].lower()+str(p.resolve()).replace('\\','/')[2:]
    def expr(row):
        # SoPlex's LP reader has a fixed 8190-character line limit.
        terms=[('+' if x>=0 else '-')+' '+str(abs(x))+' '+name for name,x in zip(variables,row) if x]
        return '\n  '.join(' '.join(terms[i:i+8]) for i in range(0,len(terms),8)) or '0'
    lines=['Maximize',' obj: '+expr(objective),'Subject To']
    for prefix,rows,sign in (('i',inequalities,'<='),('e',equalities,'=')):
        lines.extend(f' {prefix}{j}: {expr(row)} {sign} {rhs}' for j,(row,rhs,_) in enumerate(rows))
    lines+=['Bounds']+[f' 0 <= {name} <= {u}' for name,u in zip(variables,upper)]+['End']
    # This is a generated solver input, not an edit to any source/model file.
    with tempfile.TemporaryDirectory(prefix='spin-dual-') as folder:
        path=Path(folder)/'shell.lp';path.write_text('\n'.join(lines)+'\n',encoding='ascii')
        command=['wsl.exe','-d','Ubuntu-24.04','--','env',
                 'LD_LIBRARY_PATH='+linux(runtime/'usr/lib/x86_64-linux-gnu'),linux(runtime/'usr/bin/soplex'),
                 '--arithmetic=2','--precision=100','--readmode=1','--solvemode=2',
                 '--int:syncmode=1','--int:checkmode=2','-f0','-o0','-t60','-s0','-Y','-c']
        if tight:
            # The runtime's saved settings retain double-scale pivot/zero
            # tolerances even with 100-decimal-digit arithmetic.
            command += [f'--real:{name}=1e-80' for name in
                        ('epsilon_zero','epsilon_factorization','epsilon_update','epsilon_pivot')]
            command += ['--real:fpfeastol=1e-12','--real:fpopttol=1e-12']
        command.append(linux(path))
        print('Starting bounded 60-second multiprecision/rational solver',flush=True)
        process=subprocess.run(command,capture_output=True,text=True,timeout=90)
    output=process.stdout+process.stderr
    print(output[-16000:],flush=True)
    return None


def shell(data,w,old,direct=False,limit=1e6,soplex=False):
    variables,scales,upper,inequalities,equalities,fixed,table=data
    unit=1<<(old.bit_length()-1)
    physical={f'{p}_{v}':Q(table[v][w]*(2 if v<128 else 1)*(31 if p=='h' else 1),1<<128)
              for p in ('q','h') for v in table}
    objective=[physical[x]*s/unit for x,s in zip(variables,scales)]
    constant=sum((physical[x]*v for x,v in fixed.items()),Q(0))/unit
    if soplex:
        exact_solver(objective,constant,upper,inequalities,equalities,variables)
        return old
    if direct:
        from types import SimpleNamespace as NS
        # A bounded search for dual multipliers is feasible at zero prices
        # with positive column residuals. It avoids a delicate primal start.
        rows=inequalities+equalities;ni=len(inequalities);ne=len(equalities);nv=len(variables)
        matrix=np.array([[float(x) for x in row] for row,_,_ in rows])
        cost=np.array([float(rhs) for _,rhs,_ in rows]+list(map(float,upper)))
        fit=linprog(cost,A_ub=np.hstack((-matrix.transpose(),-np.eye(nv))),
                    b_ub=-np.array(list(map(float,objective))),
                    bounds=[(0,limit)]*ni+[(-limit,limit)]*ne+[(0,None)]*nv,
                    method='highs',options={'time_limit':30})
        if not fit.success:
            print('No bounded dual witness',w,fit.message,flush=True);return old
        result=NS(success=True,ineqlin=NS(marginals=-fit.x[:ni]),
                  eqlin=NS(marginals=-fit.x[ni:ni+ne]),upper=NS(marginals=np.zeros(nv)),fun=-fit.fun)
    else:
        result=linprog(-np.array(list(map(float,objective))),
                   A_ub=np.array([[float(x) for x in row] for row,_,_ in inequalities]),
                   b_ub=np.array([float(rhs) for _,rhs,_ in inequalities]),
                   A_eq=np.array([[float(x) for x in row] for row,_,_ in equalities]),
                   b_eq=np.array([float(rhs) for _,rhs,_ in equalities]),
                   bounds=[(0,float(x)) for x in upper],method='highs',
                   options={'time_limit':30,'dual_feasibility_tolerance':1e-9,'primal_feasibility_tolerance':1e-9})
    if not result.success:
        print('No witness',w,result.message,flush=True);return old
    bound,correction=repaired_upper(objective,constant,upper,inequalities,equalities,result)
    if direct:
        prices=reconstructed_prices(rows,objective,fit,ni,limit)
        if prices is not None:
            exact,repair=price_upper(objective,constant,upper,inequalities,equalities,prices[:ni],prices[ni:])
            if exact<bound:bound,correction=exact,repair
    physical_bound=bound*unit
    assert physical_bound>=0
    checked=min(old,physical_bound.numerator//physical_bound.denominator)
    print('DUAL SHELL',w,'log2 old/checked',log2(old),log2(checked) if checked else '-inf',
          'numerical scaled upper',float(constant)-result.fun,
          'exact repair scaled',float(correction),'checked scaled',float(bound),flush=True)
    return checked


def self_test():
    # Arbitrarily rounded proposals must still upper-bound every feasible point.
    from types import SimpleNamespace as NS
    constraints=[([Q(1),Q(2)],Q(3),'toy')]
    equalities=[([Q(1),Q(-1)],Q(0),'eq')]
    objective=[Q(3),Q(1)];upper=[Q(2),Q(2)]
    for p in (-2.,-.3,0.,.7):
        for e in (-1.1,0.,2.3):
            result=NS(ineqlin=NS(marginals=[p]),eqlin=NS(marginals=[e]),upper=NS(marginals=[-.03,.04]))
            bound,_=repaired_upper(objective,Q(2),upper,constraints,equalities,result)
            assert bound>=6  # maximum at x=y=1, plus constant 2
    print('Exact residual-repair tests: 12 deliberately imprecise witnesses passed',flush=True)


def propagate(data):
    from shortened_bound import dimension_caps
    from shortening_polynomial import improve_dimensions
    from dual_shortening import improve_dimensions as dual_dimensions
    from bch_joint_support import support_caps
    from basis_lattice import improve_caps
    from shortening_moments import improve
    from dual_moments import refine_bch
    new=refined_spectrum(data)
    spectrum=authenticated_caps()
    dimensions=dual_dimensions(improve_dimensions(dimension_caps()))
    primal=improve_caps(support_caps(spectrum,g=4,dimensions=dimensions),spectrum)
    primal,_=improve(primal,dimensions)
    before=refine_bch(primal,dimensions)
    after=refine_bch(primal,dimensions,dual_spectrum=new)
    assert all(b<=a for rowa,rowb in zip(before,after) for a,b in zip(rowa,rowb))
    print('Propagated support CDF log2 before/after:',
          [(u,log2(sum(row[u] for row in before)),log2(sum(row[u] for row in after)))
           for u in (96,112,128,144,160,176,192)],flush=True)
    print('No inner probability or full occupancy certificate.',flush=True)


def refined_spectrum(data=None):
    data=model() if data is None else data
    old=dual_shell_caps();new=old[:]
    for w,limit in ((30,100),(38,100000),(56,1000000),(64,1000000)):
        new[w]=new[256-w]=shell(data,w,old[w],True,limit)
    return new


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--weights',nargs='+',type=int,default=[30,32,38,48,56,64])
    parser.add_argument('--direct-dual',action='store_true')
    parser.add_argument('--multiplier-limit',type=float,default=1e6)
    parser.add_argument('--soplex',action='store_true')
    parser.add_argument('--propagate',action='store_true')
    args=parser.parse_args();self_test();data=model();old=dual_shell_caps()
    if args.propagate:propagate(data)
    else:
        for w in args.weights:
            assert 30<=w<=226 and w%2==0
            shell(data,w,old[w],args.direct_dual,args.multiplier_limit,args.soplex)
