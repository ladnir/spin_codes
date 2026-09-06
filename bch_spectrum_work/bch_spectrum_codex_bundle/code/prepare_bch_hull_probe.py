"""Couple P/Q to their nested doubly-even hulls, including signed transforms."""
import json
import math
import sys
from fractions import Fraction
from pathlib import Path
sys.dont_write_bytecode = True
from audit_bch_closure_envelope import kraw_table
from audit_bch_hulls import build as hull_audit
from bch_quotient import binary_poly_divmod, multiply_by_x_mod
from verify_scaled_rational_solution import normalized_rows
from export_scaled_rational_lp import expression
from run_higher_endpoint_preflight import write_new, sha

ROOT = Path(__file__).resolve().parents[1]
FOLDER = ROOT/'generated/oa21_hull_probe'


def build():
    structure = hull_audit()
    assert structure == json.loads((ROOT/'generated/bch256_hull_structure.json').read_text())
    assert structure['containments']['HP']['HQ']
    assert structure['intersection_dimensions']['P']['HQdual'] == 123
    g1 = int(structure['codes']['HP']['generator_hex'], 16)
    g2 = int(structure['codes']['HQ']['generator_hex'], 16)
    factor, rem = binary_poly_divmod(g1, g2)
    assert rem == 0 and factor.bit_length() == 9
    seen, x = set(), 1
    while x not in seen:
        seen.add(x)
        x = multiply_by_x_mod(x, factor)
    assert seen == set(range(1, 256)) and x == 1
    model = json.loads((ROOT/'generated/oa21_mod4_probe/model.json').read_text())
    variables = model['metadata']['variables']
    variables += [f'{a}_{w}' for a in ('r', 's') for w in range(0,129,4)]
    kt = kraw_table()

    def add(name, terms, sense='ge', rhs=0):
        out = {}
        for scalar, values in terms:
            for n, v in values.items():
                out[n] = out.get(n, 0) + scalar*v
        model['constraints'].append(dict(name=name, coeffs={n:str(v) for n,v in out.items() if v},
                                        sense=sense, rhs=str(rhs)))

    def point(a, w):
        return {f'{a}_{w}': 1} if f'{a}_{w}' in variables else {}

    def transform(a, j, signed=False):
        return {n: (kt[j][w]+(kt[j][256-w] if w != 128 else 0)) *
                   ((-1)**(w//2) if signed else 1)
                for n in variables if n.startswith(a+'_') for w in [int(n.split('_')[1])]}

    for a in ('r','s'):
        add(a+'_mass', [(1,transform(a,0))], 'eq', 1 << 85)
        for w in range(0,40,4):
            add(a+'_support_'+str(w), [(1,point(a,w))], 'eq', int(a=='r' and w==0))
    for j in range(0,129,2):
        kr, ks, kq, kh = (transform(a,j) for a in ('r','s','q','h'))
        q, h, r, s = (point(a,j) for a in ('q','h','r','s'))
        # HQ subset Q; HP subset Pdual; HQ subset Qdual.
        add('HQ_le_Q_'+str(j), [(1,q),(-1,r),(-255,s)])
        add('HP_le_Pdual_'+str(j), [(1,kq),(255,kh),(-(1<<131),r)])
        add('HQ_le_Qdual_'+str(j), [(1,kq),(-(1<<123),r),(-255*(1<<123),s)])
        # Q subset HQdual; matching nonzero P/Q cosets subset HPdual/HQdual cosets.
        add('Q_le_HQdual_'+str(j), [(1,kr),(255,ks),(-(1<<93),q)])
        add('H_le_Hdualcoset_'+str(j), [(1,kr),(-1,ks),(-(1<<93),h)])
        # Hull-dual distance >= 8, established by the cyclic root audit.
        if 0 < j < 8:
            add('HPdual_zero_'+str(j), [(1,kr)], 'eq')
            add('HQdual_zero_'+str(j), [(1,kr),(255,ks)], 'eq')
        for label, dim, hd, spectrum, dual, hullspec, hull_transform in (
            ('P',131,85,[(1,q),(255,h)],[(1,kq),(255,kh)],[(1,r)],[(1,kr)]),
            ('Q',123,93,[(1,q)],[(1,kq)],[(1,r),(255,s)],[(1,kr),(255,ks)])):
            signed = [(1,transform('q',j,True))]
            if label == 'P':
                signed.append((255,transform('h',j,True)))
            # G_j = signed Krawtchouk transform. Fourier values on hull-dual
            # are +/- 2^108, and zero outside. Count each sign separately.
            # Positive sign includes all weight-2(mod4) words of E.
            pos = [(c*(1<<(dim+108-hd)),v) for c,v in hull_transform]
            pos += [(c*(1<<dim),v) for c,v in signed]
            if j % 4 == 2:
                pos += [(-c*(1<<(dim+109)),v) for c,v in spectrum]
            add('Fourier_positive_'+label+'_'+str(j), pos)
            # Negative sign contains Edual union the 0(mod4) part of E.
            # Their intersection is the hull, whose spectrum is subtracted once.
            neg = [(c*(1<<(dim+108-hd)),v) for c,v in hull_transform]
            neg += [(-c*(1<<dim),v) for c,v in signed]
            neg += [(-c*(1<<109),v) for c,v in dual]
            neg += [(c*(1<<(dim+109)),v) for c,v in hullspec]
            if j % 4 == 0:
                neg += [(-c*(1<<(dim+109)),v) for c,v in spectrum]
            add('Fourier_negative_'+label+'_'+str(j), neg)
    model['metadata']['hull_meaning'] = 'r=HP spectrum, s=common nonzero HQ/HP coset spectrum; both half spectra at multiples of four'
    model['metadata']['hull_quotient_factor_hex'] = hex(factor)
    model['metadata']['hull_sources_sha256'] = {str(p.relative_to(ROOT)):sha(p) for p in
        (Path(__file__),ROOT/'generated/bch256_hull_structure.json',ROOT/'generated/oa21_mod4_probe/model.json')}
    scales = {n:max(1,math.comb(256,int(n.split('_')[1]))//(1<<(170 if n[0] in 'rs' else 132)))
              for n in variables}
    return model, scales


def export(model, scales):
    rows = normalized_rows(model,scales)
    lines = ['Maximize',' obj: h_38','Subject To']
    for i,row in enumerate(rows,1):
        terms = list(row['coeffs'].items())
        # Wrap only between complete coefficient-variable terms.
        pieces = [expression([(v,n)]) for n,v in terms if v]
        chunks = [' '.join(pieces[j:j+8]) for j in range(0,len(pieces),8)]
        if not chunks:
            chunks = ['0 h_38']
        chunks[0] = f' c{i}: '+chunks[0]
        chunks[-1] += ' '+{'eq':'=','ge':'>=','le':'<='}[row['sense']]+' '+str(row['rhs'])
        lines.extend(chunks)
    lines += ['Bounds']+[f' 0 <= {n}' for n in model['metadata']['variables']]+['End']
    return '\n'.join(lines)+'\n',len(rows)


if __name__ == '__main__':
    model,scales = build()
    lp,count = export(model,scales)
    FOLDER.mkdir(exist_ok=False)
    write_new(FOLDER/'model.json',model)
    write_new(FOLDER/'scales.json',scales)
    with (FOLDER/'h_38.lp').open('x') as f:
        f.write(lp)
    print(json.dumps(dict(variables=len(scales),rows=count)))
