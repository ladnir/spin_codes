"""One bounded exact LP payoff test for D=E dual containing P containing Q.

D is the even extension of the designed-distance-31 BCH code. Its nonzero
P-cosets have a common spectrum u; E=D dual has certified distance >=32.
"""
import argparse
import copy
import math
from pathlib import Path
import bridge as base
import low_shell_roi as roi
import low_shell_constraint_probe as probe
import pdual_case_audit as cover
from bch_quotient import generator_polynomial,binary_poly_divmod,multiply_by_x_mod
from affine_wambach import gf_pow
from prepare_bch_hull_probe import export
from export_scaled_rational_lp import expression
from bch_hull_cut_iteration import warm_basis,solve


def structure():
    saved=base.read(base.HERE/'generated/pdual_partial_cover_v1.json')
    assert cover.audit(saved['forest'])==saved['audit']
    assert saved['audit']['extended_F7_kernel']['minimum_distance_lower']>=32
    gp=generator_polynomial(37);gd=generator_polynomial(31)
    assert 255-(gp.bit_length()-1)==131 and 255-(gd.bit_length()-1)==139
    factor,rem=binary_poly_divmod(gp,gd)
    assert rem==0 and factor.bit_length()==9
    seen=set();x=1
    while x not in seen:seen.add(x);x=multiply_by_x_mod(x,factor)
    assert seen==set(range(1,256)) and x==1
    # P plus the eight F7 check vectors is E dual, and equals D by dimension.
    coords=[gf_pow(2,i) for i in range(255)]+[0]
    values=[gf_pow(x,7) if x else 0 for x in coords]
    rows=[gp<<i for i in range(131)]
    rows=[w|((w.bit_count()%2)<<255) for w in rows]
    checks=[sum(1<<i for i,v in enumerate(values) if v>>j&1) for j in range(8)]
    pivots={}
    for w in rows+checks:
        assert w.bit_count()%2==0
        assert binary_poly_divmod(w&((1<<255)-1),gd)[1]==0
        while w:
            j=w.bit_length()-1
            if j not in pivots:pivots[j]=w;break
            w^=pivots[j]
    assert len(pivots)==139
    return dict(D_dimension=139,D_designed_distance=31,D_extended_distance_lower=32,
        E_dimension=117,E_distance_lower=32,quotient_factor_hex=hex(factor),
        nonzero_P_coset_orbit_size=255,E_dual_equals_D_checked=True)


def prepare():
    facts=structure();ref=base.BCH/'generated/shift_rank_oa29_joint'
    old=base.read(ref/'model.json');model=copy.deepcopy(old);scales=base.read(ref/'scales.json')
    weights=range(0,129,2);kt=roi.kraw_table()
    model['metadata']['variables'] += [f'u_{w}' for w in weights]
    scales.update({f'u_{w}':max(1,math.comb(256,w)//(1<<126)) for w in weights})
    def transform(terms,j):
        return {f'{a}_{w}':mult*(kt[j][w]+(kt[j][256-w] if w!=128 else 0))
                for a,mult in terms for w in weights}
    def add(name,co,sense='ge',rhs=0):
        model['constraints'].append(dict(name='KERNEL_'+name,coeffs={n:str(v) for n,v in co.items() if v},sense=sense,rhs=str(rhs)))
    add('u_mass',transform([('u',1)],0),'eq',1<<131)
    for w in range(0,32,2):add(f'u_support_{w}',{f'u_{w}':1},'eq')
    for j in weights:
        e=transform([('q',1),('h',255),('u',255)],j)
        h=transform([('q',1),('h',255),('u',-1)],j)
        add(f'E_spectrum_{j}',e,'eq' if 0<j<32 else 'ge')
        add(f'nonzero_E_coset_spectrum_{j}',h,'eq' if 0<j<30 else 'ge')
        # E is contained in P dual and contains HP=P intersect P dual:
        # P has F7 zero, hence HP does too. Thus D is contained in HP dual.
        if j%4==0:
            c=e.copy();c[f'r_{j}']=-(1<<139);add(f'HP_subset_E_{j}',c)
        c={f'r_{w}':kt[j][w]+(kt[j][256-w] if w!=128 else 0) for w in range(0,129,4)}
        c.update({f'q_{j}':-(1<<85),f'h_{j}':-255*(1<<85),f'u_{j}':-255*(1<<85)})
        add(f'D_subset_HPdual_{j}',c)
    # Borrow only the existing rounded current objective, not its extra hypothesis.
    _,_,objective,norm,_,_,_=probe.prepare('pdual30_zero')
    lp=export(model,scales)[0].replace(' obj: h_38\n',' obj: '+expression([(v,n) for n,v in objective.items()])+'\n')
    warm,mapping=warm_basis(old,base.read(ref/'scales.json'),model,scales,ref/'h_38.bas')
    return model,scales,objective,norm,lp,warm,mapping,facts


def run(args):
    folder=base.HERE/'generated/kernel_coset_payoff_v1'
    model,scales,objective,norm,lp,warm,mapping,facts=prepare()
    if not args.verify:
        folder.mkdir(exist_ok=False)
        for name,text in (('h_38.lp',lp),('warm.bas',warm)):
            with (folder/name).open('x',encoding='ascii') as stream:stream.write(text)
        base.write_new(folder/'mapping.json',mapping)
        print(solve(folder,args.seconds,True),flush=True)
    assert (folder/'h_38.lp').read_text()==lp
    if not (folder/'h_38.sol').exists() or 'status = OPTIMAL' not in (folder/'h_38.sol').read_text():
        print('Bounded exact solve inconclusive; no claimed improvement.',flush=True);return
    result=probe.audit(folder,model,scales,objective,norm)
    result['status']='EXACT_KERNEL_COSET_LP_PAYOFF_WITH_CHECKED_STRUCTURE'
    result['structure']=facts
    result['source_sha256'].update({p.relative_to(base.ROOT).as_posix():base.sha(p) for p in
        (Path(__file__),base.HERE/'generated/pdual_partial_cover_v1.json',Path(cover.__file__))})
    if args.verify:assert result==base.read(folder/'audit.json')
    else:base.write_new(folder/'audit.json',result)
    print({k:result[k] for k in ('status','conditional_margin_bits','improvement_bits','fixed_coefficients_gain_ceiling_bits','structure')},flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--seconds',type=int,default=60)
    p.add_argument('--verify',action='store_true');run(p.parse_args())
