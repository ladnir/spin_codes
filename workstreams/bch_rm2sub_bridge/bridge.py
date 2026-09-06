"""Isolated BCH/RM2Sub occupation-one bridge; never writes another worktree.

Snapshot selected maps once, verify their algebra, screen per-weight tilts,
then recompute chosen witnesses with Arb. The existing BCH weighted bound
is reused by coefficient domination, not by substituting the M22 inner law.
All outputs are write-once. No occupation >=2 claim is made.
"""
import argparse
import hashlib
import json
import math
import sys
from collections import Counter
from fractions import Fraction as F
from pathlib import Path

import numpy as np

sys.dont_write_bytecode = True
sys.set_int_max_str_digits(0)
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
BCH = ROOT / 'bch_spectrum_work/bch_spectrum_codex_bundle'
CONFIGS = ('t64_s16', 't128_s15', 't256_s14')
WEIGHTS = tuple(range(38, 220, 2)) + (256,)
ROWS, LENGTH, CUTOFF = 8192, 256, 209716


def read(path):
    return json.loads(path.read_text())


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_new(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('x', encoding='utf-8') as out:
        json.dump(value, out, indent=2)
        out.write('\n')


def encode(x):
    return {'numerator': str(x.numerator), 'denominator': str(x.denominator)}


def decode(x):
    return F(int(x['numerator']), int(x['denominator']))


def snapshot(source):
    directory = source / 'workstreams/finite_asymptotic_theory/small_k_replay/rm2sub_calibration_constituents'
    target = HERE / 'inputs'
    assert not target.exists(), 'Inputs are immutable; use a new snapshot directory for another version'
    payloads, hashes = {}, {}
    for name in CONFIGS:
        for suffix in ('selection', 'a_spectrum', 'b_kernel_spectrum'):
            filename = f'{name}_{suffix}.json'
            path = directory / filename
            payloads[filename] = read(path)
            hashes[str(path.relative_to(source))] = sha(path)
    target.mkdir()
    for name, payload in payloads.items():
        write_new(target / name, payload)
    reference = source / 'workstreams/finite_asymptotic_theory/small_k_replay/certify_rm2sub_rm49_q1_outward.py'
    write_new(target / 'manifest.json', dict(source_worktree=str(source.resolve()),
        source_sha256=hashes, reference_checker=dict(path=str(reference.relative_to(source)), sha256=sha(reference)),
        snapshot_sha256={name: sha(target / name) for name in payloads}))
    print('Snapshot created: nine map/spectrum JSON inputs; no source-worktree writes', flush=True)


def load_map(name):
    manifest = read(HERE / 'inputs/manifest.json')
    for filename, expected in manifest['snapshot_sha256'].items():
        assert sha(HERE / 'inputs' / filename) == expected
    selection = read(HERE / 'inputs' / f'{name}_selection.json')
    t, s = (selection['parameters'][key] for key in ('step_bits', 'state_bits'))
    assert name == f't{t}_s{s}' and ROWS % t == 0
    selected = selection['selected']
    generators = [int(v, 16) for v in selected['A_generator_words_hex']]
    columns = [int(v, 16) for v in selected['B_columns_hex']]
    assert len(generators) == s and len(columns) == t
    assert len(set(columns)) == t and 0 not in columns
    assert all(0 < c < 1 << s for c in columns)
    assert generators == [sum(((c >> j) & 1) << i for i, c in enumerate(columns)) for j in range(s)]
    assert all((a & b).bit_count() % 2 == 0 for a in generators for b in generators)
    pivots = {}
    for v in generators:
        assert v < 1 << t
        while v:
            p = v.bit_length() - 1
            if p not in pivots:
                pivots[p] = v
                break
            v ^= pivots[p]
    assert len(pivots) == s
    # Gray-code enumeration changes one generator per word; no word list stored.
    spectrum = Counter({0: 1})
    word = 0
    for i in range(1, 1 << s):
        word ^= generators[(i & -i).bit_length() - 1]
        spectrum[word.bit_count()] += 1
    claimed = {row['weight']: row['count'] for row in read(HERE / 'inputs' / f'{name}_a_spectrum.json')['spectrum'] if row['count']}
    assert dict(spectrum) == claimed
    kernel = {row['total_weight']: row['kernel_words'] for row in read(HERE / 'inputs' / f'{name}_b_kernel_spectrum.json')['by_total_weight']}
    for degree in range(t + 1):
        total = 0
        for w, count in spectrum.items():
            kraw = sum((-1) ** j * math.comb(w, j) * math.comb(t-w, degree-j)
                       for j in range(max(0, degree-(t-w)), min(degree, w)+1))
            total += count * kraw
        value, remainder = divmod(total, 1 << s)
        assert remainder == 0 and value == kernel.get(degree, 0)
    assert spectrum[0] == 1 and kernel.get(1, 0) == 0
    return t, s, {w: count for w, count in sorted(spectrum.items()) if w}


def lse(values):
    values = np.asarray(values)
    return float(np.logaddexp.reduce(values))


def log_regions(t, s, spectrum, lam):
    den = (1 << s) - 1
    log_m0 = lse([math.log(c / den) - lam*w for w, c in spectrum.items()])
    terms = []
    for w, c in spectrum.items():
        if w:
            terms.append(math.log(c / den * w / t) - lam*(w-1))
        if w < t:
            terms.append(math.log(c / den * (t-w) / t) - lam*(w+1))
    log_m1 = lse(terms)
    log_kappa = math.log1p(1 / (den-1))
    a, b = log_kappa+log_m0, log_kappa+log_m1
    epochs = ROWS // t
    geom = lse(np.arange(epochs)*a) - math.log(epochs)
    return epochs*a, (-lam+geom, b-math.log(den)+geom, b+(epochs-1)*a)


def log_coefficients(z11, active):
    a01, a10, a11 = active
    current = np.full((LENGTH+1, 2), -np.inf)
    current[0, 0] = 0.
    for n in range(LENGTH):
        updated = np.full_like(current, -np.inf)
        old = current[:n+1]
        updated[:n+1, 0] = old[:, 0]
        updated[:n+1, 1] = old[:, 1] + z11
        updated[1:n+2, 0] = np.logaddexp(updated[1:n+2, 0], old[:, 1]+a10)
        updated[1:n+2, 1] = np.logaddexp(updated[1:n+2, 1],
            np.logaddexp(old[:, 0]+a01, old[:, 1]+a11))
        current = updated
    return np.logaddexp(current[:, 0], current[:, 1]) - np.array([math.log(math.comb(LENGTH, w)) for w in range(LENGTH+1)])


def bch_bound(coefficients):
    # Bound only Q1. No M22 transfer coefficient is treated as an RM2Sub coefficient.
    sys.path.insert(0, str(BCH / 'code'))
    from certify_bch_m25_closure import deterministic_caps
    old = read(BCH / 'generated/shift_rank_oa29_joint/objective.json')
    audit = read(BCH / 'generated/shift_rank_oa29_joint/audit.json')
    pairs = {w: coefficients[w]+coefficients[256-w] for w in (38, 40, 42)}
    factor = max(pairs[w]/decode(old['dyadic_upper_pair_coefficients'][str(w)]) for w in pairs)
    bound = factor*decode(audit['paired_shells_upper'])
    caps = deterministic_caps()
    excluded = {38, 40, 42, 214, 216, 218}
    rest = sum((caps[w]*coefficients[w] for w in WEIGHTS if w not in excluded), F(0))
    return bound+rest, factor, rest


def screen(name):
    t, s, spectrum = load_map(name)
    best = np.full(LENGTH+1, np.inf)
    witness = np.zeros(LENGTH+1, dtype=int)
    for tenth in range(-120, 1):
        lam = math.exp(tenth/10)
        values = log_coefficients(*log_regions(t, s, spectrum, lam)) + CUTOFF*lam
        values = math.log(ROWS) + np.minimum(0., values)
        improved = values < best
        best = np.minimum(values, best)
        witness[improved] = tenth
    # Diagnostic dyadics only. Underflow is rounded to a positive placeholder;
    # certification recomputes the selected witnesses independently with Arb.
    coefficients = {w: F.from_float(math.exp(max(-700., float(best[w])))) for w in WEIGHTS}
    upper, factor, rest = bch_bound(coefficients)
    payload = dict(status='BINARY64_Q1_SCREEN_ONLY', configuration=name,
        parameters=dict(message_bits=1 << 20, output_bits=1 << 21, outer_length=LENGTH, outer_rows=ROWS,
                        distance_cutoff=CUTOFF, step_bits=t, state_bits=s),
        coefficient_rows={str(w): dict(log_coefficient=float(best[w]), witness_tenth=int(witness[w])) for w in WEIGHTS},
        diagnostic_margin_bits=math.log2(upper.denominator)-math.log2(upper.numerator),
        dual_domination_factor=float(factor), remaining_shell_fraction_diagnostic=float(rest/upper),
        all_occupations_certified=False, source_sha256={'bridge.py':sha(Path(__file__)), 'inputs/manifest.json':sha(HERE/'inputs/manifest.json')})
    write_new(HERE / 'generated' / f'{name}_screen.json', payload)
    print(name, 'Q1 screen margin', payload['diagnostic_margin_bits'], 'dual factor',float(factor), flush=True)


def arb_coefficients(t, s, spectrum, tenth):
    from flint import arb
    lam = (arb(tenth)/10).exp()
    z = (-lam).exp()
    den = (1 << s)-1
    m0 = sum((arb(c)*z**w/den for w, c in spectrum.items()), arb(0))
    m1 = sum((arb(c)*(w*z**(w-1)+(t-w)*z**(w+1))/(t*den) for w, c in spectrum.items()), arb(0))
    kappa = arb(den)/(den-1)
    zero = (arb(1), arb(0), arb(0), kappa*m0)
    active = (arb(0), z, kappa*m1/den, kappa*m1)
    # Independent positive matrix products rather than the screening formulas.
    def mul(a,b):
        return (a[0]*b[0]+a[1]*b[2], a[0]*b[1]+a[1]*b[3],
                a[2]*b[0]+a[3]*b[2], a[2]*b[1]+a[3]*b[3])
    rz, ra = (arb(1),arb(0),arb(0),arb(1)), (arb(0),)*4
    for _ in range(ROWS//t):
        first, second = mul(ra,zero), mul(rz,active)
        ra = tuple(a+b for a,b in zip(first,second))
        rz = mul(rz,zero)
    ra = tuple(a/(ROWS//t) for a in ra)
    current = [(arb(1),arb(0))]
    for n in range(LENGTH):
        updated = []
        for w in range(n+2):
            a,b = arb(0),arb(0)
            if w<=n:
                x,y = current[w]
                a += (x*rz[0]+y*rz[2])*(n+1-w)/(n+1)
                b += (x*rz[1]+y*rz[3])*(n+1-w)/(n+1)
            if w:
                x,y = current[w-1]
                a += (x*ra[0]+y*ra[2])*w/(n+1)
                b += (x*ra[1]+y*ra[3])*w/(n+1)
            updated.append((a,b))
        current = updated
    correction = ROWS*(CUTOFF*lam).exp()
    return [(a+b)*correction for a,b in current]


def certify(name):
    raise RuntimeError('Two-state activation is not justified. Use activation_bridge.py; old screens are historical diagnostics only.')
    from flint import ctx
    sys.path.insert(0, str(BCH/'code'))
    from audit_bch_q1_full_arb import rational
    t,s,spectrum = load_map(name)
    discovery_path = HERE/'generated'/f'{name}_screen.json'
    discovery = read(discovery_path)
    assert discovery['source_sha256']['bridge.py'] == sha(Path(__file__))
    groups = {}
    for w in WEIGHTS:
        tenth = discovery['coefficient_rows'][str(w)]['witness_tenth']
        assert isinstance(tenth,int) and -120<=tenth<=0
        groups.setdefault(tenth,[]).append(w)
    ctx.prec = 256
    co = {}
    for tenth, weights in sorted(groups.items()):
        values = arb_coefficients(t,s,spectrum,tenth)
        for w in weights:
            co[w] = min(F(ROWS),rational(values[w].upper()))
            assert co[w]>0
        print(name,'Arb weights',weights,flush=True)
    upper, factor, rest = bch_bound(co)
    payload = dict(status='OUTWARD_Q1_BOUND_USING_EXISTING_BCH_DUAL', configuration=name,
        parameters=discovery['parameters'], coefficient_upper={str(w):encode(v) for w,v in co.items()},
        Q1_upper=encode(upper), dual_domination_factor=encode(factor), remaining_shell_upper=encode(rest),
        margin_bits_diagnostic=math.log2(upper.denominator)-math.log2(upper.numerator),
        Q1_below_2_to_minus_40=upper<F(1,1<<40), all_occupations_certified=False,
        assumption_scope='Uses the reviewed two-state RM2Sub envelope and retained BCH exact-dual certificate; not an M22 probability comparison',
        source_sha256={'bridge.py':sha(Path(__file__)), 'inputs/manifest.json':sha(HERE/'inputs/manifest.json'),
            f'generated/{name}_screen.json':sha(discovery_path),
            'BCH_joint_audit':sha(BCH/'generated/shift_rank_oa29_joint/audit.json'),
            'BCH_joint_objective':sha(BCH/'generated/shift_rank_oa29_joint/objective.json')})
    write_new(HERE/'generated'/f'{name}_q1_outward.json',payload)
    print(name,'OUTWARD Q1 margin',payload['margin_bits_diagnostic'],flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('mode',choices=('snapshot','screen','certify','audit-map'))
    parser.add_argument('--source',type=Path)
    parser.add_argument('--configuration',choices=CONFIGS)
    args=parser.parse_args()
    if args.mode=='snapshot':
        assert args.source is not None
        snapshot(args.source)
    else:
        assert args.configuration is not None
        if args.mode=='screen': screen(args.configuration)
        elif args.mode=='certify': certify(args.configuration)
        else:
            t,s,spectrum=load_map(args.configuration)
            print('Exact A/B, full A spectrum, and MacWilliams kernel checks passed', t,s)
