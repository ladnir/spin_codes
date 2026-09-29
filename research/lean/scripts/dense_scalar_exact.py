"""Generate an untrusted rational scalar witness and a kernel-checked local box.

Python chooses z and an upper radius. Lean proves six exact integer inequalities
and four logarithmic vertex checks; no Python numerical answer is trusted.
"""
from decimal import Decimal, localcontext
from fractions import Fraction as Q
import argparse
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT.parent/'workstreams/inner_design/imt_asymptotic/d11'


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('index', type=int)
    args = parser.parse_args()
    leaves = json.loads((SOURCE/'DENSE_REPLAY.json').read_text())['leaves']
    boxes = [b for b in leaves if b['rational_witness']['family'] == 'statefree']
    box = boxes[args.index]
    w = box['rational_witness']
    p, y = Q(w['p']), Q(w['y'])
    q = p*y
    with localcontext() as ctx:
        ctx.prec = 100
        ell = Q(w['log_lam'])
        value = (-(Decimal(ell.numerator)/Decimal(ell.denominator)).exp()).exp()
        z = Q(int(value*10**30), 10**30)
    assert 0 < q < 1 and 0 < z < 1
    g0, g1 = 1-q+q*z, q+(1-q)*z
    exact = max(g0**(128-d)*g1**d for d in [0,48,56,64,72,80])
    r = Q(-((-exact.numerator*10**100)//exact.denominator), 10**100)
    assert 0 < exact <= r < 1
    outer = json.loads((SOURCE/'OUTER_REFINED.json').read_text())
    seg = box['segment']
    if seg < 19:
        a = outer['left_segments'][seg]; m, c = Q(a['slope']), Q(a['intercept'])
    elif seg == 19:
        m, c = Q(0), Q(outer['central_constant_upper'])
    else:
        a = outer['left_segments'][38-seg]; m, c = -Q(a['slope']), Q(a['slope'])+Q(a['intercept'])
    tag = f'B{args.index:03d}'
    values = dict(m=m,c=c,p=p,y=y,radius=r,z=z,a0=Q(box['alpha'][0]),a1=Q(box['alpha'][1]),
                  x0=Q(box['row_density'][0]),x1=Q(box['row_density'][1]))
    data = ['import SpinCodes.Structured.DenseScalarExactDefs',
            'import SpinCodes.Structured.DenseOccupationFixedVertexDefs', '',
            f'namespace Spin.Structured.DenseScalarExact.{tag}',
            'open Spin.Structured.DenseOccupationFixed',
            'set_option maxHeartbeats 0', 'set_option maxRecDepth 100000']
    for name, value in values.items():
        data.append(f'def {name} : QInput := ⟨{value.numerator}, {value.denominator}⟩')
    data.extend([f'def qn : Int := {q.numerator}', f'def qd : Int := {q.denominator}'])
    for d in [0,48,56,64,72,80]:
        data.append(f'theorem check{d} : check qn qd z.num z.den radius.num radius.den {d} = true := by decide')
    data.append('def upper : Int := -400000000000000000000000')
    for a in [0,1]:
        for x in [0,1]:
            data.append(f'theorem check_{a}{x} : vertexCheck m c p y radius z a{a} x{x} 50 upper = true := by decide')
    data.append(f'end Spin.Structured.DenseScalarExact.{tag}')
    (ROOT/f'SpinCodes/Structured/DenseScalarExact{tag}Data.lean').write_text('\n'.join(data)+'\n',encoding='utf-8')
    template = (ROOT/'SpinCodes/Structured/DenseOccupationFixedFirstBox.lean').read_text(encoding='utf-8-sig')
    start = template.index('lemma vertex00')
    end = template.index('theorem parameters_match')
    body = template[start:end].replace('first retained occupation box', f'scalar box {args.index}')
    head = [f'import SpinCodes.Structured.DenseScalarExact{tag}Data',
            'import SpinCodes.Structured.DenseScalarExactSound',
            'import SpinCodes.Structured.DenseOccupationFixedVertex', '', 'noncomputable section',
            f'namespace Spin.Structured.DenseScalarExact.{tag}',
            'open Spin.Numeric Spin.Structured.DenseOccupationFixed Set', '', body,
            'theorem scalar_bound : ConcreteScalar.scalarBound (p.real*y.real) z.real ≤ radius.real := by',
            '  have h := checked_scalar (qn := qn) (qd := qd) (zn := z.num) (zd := z.den)',
            '    (rn := radius.num) (rd := radius.den) (by decide) (by decide) (by decide) check0',
            '    (by intro i\n        fin_cases i\n        · exact check48\n        · exact check56\n        · exact check64\n        · exact check72\n        · exact check80)',
            '  have he : p.real*y.real = (qn:ℝ)/qd := by norm_num [p,y,QInput.real,qn,qd]',
            '  rw [he]', '  exact h', '', '#print axioms scalar_bound', '#print axioms exponent_bound',
            f'end Spin.Structured.DenseScalarExact.{tag}']
    (ROOT/f'SpinCodes/Structured/DenseScalarExact{tag}.lean').write_text('\n'.join(head)+'\n',encoding='utf-8')
    (ROOT/f'scripts/map_data/dense_scalar_exact_{tag}_candidate.json').write_text(json.dumps({
        'status':'UNTRUSTED_CANDIDATE', 'source_box':box, 'q':str(q), 'z':str(z), 'radius':str(r),
        'scope':'Local scalar transfer and affine exponent box; outer semantics and global coverage separate.',
    },indent=2)+'\n',encoding='utf-8')
    print(tag, 'generated; radius', float(r), 'z', float(z), 'segment', seg)


if __name__ == '__main__':
    main()
