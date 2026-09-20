"""Bounded exact LP probes of explicitly hypothetical additional BCH facts.

The LP witnesses prove implications, not the proposed additional code facts.
Run serially; every attempt and receipt has a new directory.
"""
import argparse
import copy
from fractions import Fraction as F
from pathlib import Path
from flint import fmpq
import bridge as base
import low_shell_roi as roi
from prepare_bch_hull_probe import export
from export_scaled_rational_lp import expression
from bch_hull_cut_iteration import warm_basis, solve
from verify_scaled_rational_solution import normalized_rows, parse_assignments


def prepare(candidate):
    reference = base.BCH/'generated/shift_rank_oa29_joint'
    old = base.read(reference/'model.json')
    model = copy.deepcopy(old)
    scales = base.read(reference/'scales.json')
    kt = roi.kraw_table()
    def add(terms, degree, sense='eq', rhs=0):
        coefficients = {n: str(mult*(kt[degree][w]+(kt[degree][256-w] if w != 128 else 0)))
                        for prefix, mult in terms for n in scales if n.startswith(prefix+'_')
                        for w in [int(n.split('_')[1])]}
        model['constraints'].append(dict(name=f'ROI_{candidate}_{len(model["constraints"])}',
                                         coeffs=coefficients, sense=sense, rhs=str(rhs)))
    if candidate == 'pdual30_zero':
        add([('q', 1), ('h', 255)], 30)
    elif candidate == 'qdual30_zero':
        add([('q', 1)], 30)
        add([('h', 1)], 30)
    elif candidate in ('hqdual20_zero', 'hqdual20_partial'):
        add([('r', 1), ('s', 255)], 20, 'le',
            (1 << 93)*(832972800 if candidate == 'hqdual20_partial' else 0))
    elif candidate == 'hqdual20_22_24_zero':
        for j in (20, 22, 24):
            add([('r', 1), ('s', 255)], j)
    else:
        raise ValueError(candidate)
    q1 = base.read(base.HERE/'generated/frontier_k28_q1_v1.json')
    co = {int(w): base.decode(v) for w, v in q1['coefficient_upper'].items()}
    physical = {f'{a}_{w}': (co[w]+co[256-w])*mult*scales[f'{a}_{w}']
                for w in (38,40,42) for a,mult in (('q',1),('h',31))}
    norm = max(physical.values())
    objective = {}
    for n, v in physical.items():
        ratio = v/norm
        rounded = F(-(-(ratio.numerator << 64)//ratio.denominator),1 << 64)
        assert ratio <= rounded < ratio+F(1,1 << 64)
        objective[n] = rounded
    lp = export(model,scales)[0].replace(' obj: h_38\n',' obj: '+expression([(v,n) for n,v in objective.items()])+'\n')
    warm, mapping = warm_basis(old,scales,model,scales,reference/'h_38.bas')
    return model,scales,objective,norm,lp,warm,mapping


def audit(folder, model, scales, objective, norm):
    solution = (folder/'h_38.sol').read_text()
    assert 'status = OPTIMAL' in solution
    variables = model['metadata']['variables']
    raw = parse_assignments(solution,'VARS:','REDUCED COST:')
    prices = parse_assignments(solution,'PI:','SLACK:')
    rows = normalized_rows(model,scales)
    assert set(raw) <= set(variables) and set(prices) <= {f'c{i}' for i in range(1,len(rows)+1)}
    x = {n: roi.fq(raw.get(n,F(0))) for n in variables}
    assert all(v >= 0 for v in x.values())
    combined = {n: fmpq(0) for n in variables}
    bound = fmpq(0)
    for i,row in enumerate(rows,1):
        lhs = sum((roi.fq(a)*x[n] for n,a in row['coeffs'].items()),fmpq(0))
        rhs = roi.fq(row['rhs'])
        assert {'eq':lhs == rhs,'le':lhs <= rhs,'ge':lhs >= rhs}[row['sense']],row['name']
        price = roi.fq(prices.get(f'c{i}',F(0)))
        assert row['sense'] == 'eq' or (price >= 0 if row['sense'] == 'le' else price <= 0)
        bound += price*rhs
        for n,a in row['coeffs'].items():
            combined[n] += price*roi.fq(a)
    assert all(combined[n] >= roi.fq(objective.get(n,F(0))) for n in variables)
    assert bound == sum((roi.fq(a)*x[n] for n,a in objective.items()),fmpq(0))
    q1path = base.HERE/'generated/frontier_k28_q1_v1.json'
    fullpath = base.HERE/'generated/curve_k28_full_retained_v1.json'
    co = {int(w):base.decode(v) for w,v in base.read(q1path)['coefficient_upper'].items()}
    original,_,rest = base.bch_bound(co)
    higher = base.decode(base.read(fullpath)['higher_occupancy_upper'])
    conditional = min(original,roi.fraction(bound)*norm+rest)+higher
    # A feasible point of the stronger model limits its possible full-objective gain.
    full_score = sum(((roi.fraction(x[f'q_{min(w,256-w)}'])*scales[f'q_{min(w,256-w)}']+
                       31*roi.fraction(x[f'h_{min(w,256-w)}'])*scales[f'h_{min(w,256-w)}'])*c
                      for w,c in co.items()),F(0))
    return dict(status='EXACT_LP_IMPLICATION_ONLY_ADDITIONAL_CODE_FACT_NOT_PROVED',
        candidate=folder.name,rows_checked=len(rows),variables_checked=len(variables),
        rational_primal_dual_checks_passed=True,paired_upper=base.encode(roi.fraction(bound)*norm),
        conditional_full_upper=base.encode(conditional),conditional_margin_bits=roi.margin(conditional),
        baseline_margin_bits=roi.margin(original+higher),
        improvement_bits=roi.margin(conditional)-roi.margin(original+higher),
        feasible_full_score=base.encode(full_score),
        fixed_coefficients_gain_ceiling_bits=roi.margin(full_score)-roi.margin(original),
        primal_shells={str(w):base.encode(roi.fraction(x[f'q_{w}'])*scales[f'q_{w}']+
                          31*roi.fraction(x[f'h_{w}'])*scales[f'h_{w}']) for w in (38,40,42)},
        source_sha256={p.relative_to(base.ROOT).as_posix():base.sha(p) for p in
            (Path(__file__),Path(roi.__file__),q1path,fullpath,folder/'h_38.lp',folder/'h_38.sol')})


def run(candidate, seconds, verify):
    folder = base.HERE/'generated/low_shell_constraints_v1'/candidate
    model,scales,objective,norm,lp,warm,mapping = prepare(candidate)
    if not verify:
        folder.mkdir(parents=True,exist_ok=False)
        for name,text in (('h_38.lp',lp),('warm.bas',warm)):
            with (folder/name).open('x',encoding='ascii') as stream:
                stream.write(text)
        base.write_new(folder/'mapping.json',mapping)
        print(solve(folder,seconds,True),flush=True)
    assert (folder/'h_38.lp').read_text() == lp
    if not (folder/'h_38.sol').exists():
        print('No solution: inconclusive',candidate,flush=True)
        return
    result = audit(folder,model,scales,objective,norm)
    if verify:
        assert result == base.read(folder/'audit.json')
    else:
        base.write_new(folder/'audit.json',result)
    print({k:result[k] for k in ('candidate','conditional_margin_bits','improvement_bits','fixed_coefficients_gain_ceiling_bits')},flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('candidate',choices=['pdual30_zero','qdual30_zero','hqdual20_zero','hqdual20_partial','hqdual20_22_24_zero'])
    parser.add_argument('--seconds',type=int,default=45)
    parser.add_argument('--verify',action='store_true')
    args = parser.parse_args()
    run(args.candidate,args.seconds,args.verify)
