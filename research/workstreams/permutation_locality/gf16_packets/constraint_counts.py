"""Joint shortening/dual CDF constraints with exact LP-dual verification.

Unknowns count binary subspaces by union support, for a code and its dual.
Each CDF is divided by an authenticated integer upper bound. Thus every
variable lies in [0,1]. Floating LP solutions merely propose multipliers;
dyadic arithmetic and a box-residual correction verify the final upper.
"""
import argparse
import json
from fractions import Fraction as Q
from math import comb, isfinite, log2
from pathlib import Path

import numpy as np
from scipy.optimize import linprog
from scipy.sparse import csr_matrix, vstack
import shared_support  # Establish the existing proof-module import paths.
from bch_joint_support import rank_total
from shortening_moments import gaussian
from dual_moments import split_coefficient, dual_shell_caps
from count_refinements import refine_pair


def normalize(row, rhs):
    """Power-of-two scaling preserves exact, cheap dyadic coefficients."""
    row = {i: x for i, x in row.items() if x}
    scale = 1 << max([abs(rhs).bit_length(), *(abs(x).bit_length() for x in row.values()), 0])
    return {i: Q(x, scale) for i, x in row.items()}, Q(rhs, scale)


def dual_residual(inequalities, equalities, objective, y, z, variables):
    """Return exact y.b+z.d and c-yA-zE, checking multiplier signs."""
    if (type(variables) is not int or variables < 1
            or len(y) != len(inequalities) or len(z) != len(equalities) or any(v < 0 for v in y)
            or any(type(i) is not int or not 0 <= i < variables
                   for row in [objective, *(row for row, _ in inequalities+equalities)] for i in row)):
        raise ValueError('matching dual multipliers with nonnegative inequality weights required')
    residual = [Q(objective.get(i, 0)) for i in range(variables)]
    budget = Q(0)
    for constraints, multipliers in ((inequalities, y), (equalities, z)):
        for (row, rhs), value in zip(constraints, multipliers):
            value = Q(value)
            if not value:
                continue
            budget += value*rhs
            for i, coefficient in row.items():
                residual[i] -= value*coefficient
    return budget, residual


def verify_dual(inequalities, equalities, objective, y, z, variables):
    """For Ax<=b, Ex=d, 0<=x<=1, verify c.x <= y.b+z.d+sum(r_+).

    Here y>=0 and r=c-yA-zE. No approximate feasibility is trusted.
    All inputs are exact rationals, including the rounded multipliers.
    """
    budget, residual = dual_residual(inequalities, equalities, objective, y, z, variables)
    correction = sum((max(Q(0), r) for r in residual), Q(0))
    return budget+correction, correction


class Constraints:
    def __init__(self, primal, dual, dimensions, dual_dimensions, k, last,
                 spectra=None, containment=True, complements=True, ones=(False, False)):
        n = len(dimensions)-1
        g = len(primal)
        if not 0 <= k <= n or not (n+1)//2 <= last <= n or not g or len(dual) != g:
            raise ValueError('matching ranks and a support limit covering at least half the code required')
        if (len(dual_dimensions) != n+1 or any(len(row) != n+1 for row in primal+dual)
                or any(type(v) is not int or v < 0 for row in primal+dual for v in row)
                or type(containment) is not bool or type(complements) is not bool
                or len(ones) != 2 or any(type(v) is not bool for v in ones)):
            raise ValueError('valid integer rank CDF caps and boolean constraint flags required')
        self.n = n
        self.g = g
        self.last = last
        self.caps = [[ [value//rank_total(h, g, h) for value in row]
                       for h, row in enumerate(code, 1)] for code in (primal, dual)]
        self.keys = [(side, h, u) for side in range(2) for h in range(1, g+1)
                     for u in range(last+1) if self.caps[side][h-1][u]]
        self.indices = {key: i for i, key in enumerate(self.keys)}
        self.inequalities = []
        self.equalities = []

        def add_term(row, key, coefficient):
            if key in self.indices:
                i = self.indices[key]
                side, h, u = key
                row[i] = row.get(i, 0)+coefficient*self.caps[side][h-1][u]

        def moment(row, side, h, t, coefficient):
            # Summation by parts: M_h(t)=sum_v F_h(v) C(n-v-1,t-v).
            # At t=n only the final CDF survives.
            if t == n:
                add_term(row, (side, h, n), coefficient)
            else:
                for v in range(t+1):
                    add_term(row, (side, h, v), coefficient*comb(n-v-1, t-v))

        for side, (dim, dimension) in enumerate(((dimensions, k), (dual_dimensions, n-k))):
            if ones[side]:
                # Caller establishes that the code contains the all-ones
                # word. Rank-one subcodes are its nonzero binary words.
                for u in range(max(0, n-last-1), min(last, (n-1)//2)+1):
                    row = {}
                    add_term(row, (side, 1, u), 1)
                    add_term(row, (side, 1, n-u-1), 1)
                    self.equalities.append(normalize(row, (1 << dimension)-2))
            for h in range(1, g+1):
                for u in range(last):
                    row = {}
                    add_term(row, (side, h, u), 1)
                    add_term(row, (side, h, u+1), -1)
                    if row:
                        self.inequalities.append(normalize(row, 0))
                if containment:
                    for r in range(1, h):
                        for t in range(last+1):
                            if dim[t] < r:
                                continue
                            row = {}
                            moment(row, side, h, t, gaussian(dim[t], r))
                            moment(row, side, r, t, -gaussian(dim[t], h))
                            if row:
                                self.inequalities.append(normalize(row, 0))
                if complements:
                    for t in range(max(n-dimension, n-last), last+1):
                        row = {}
                        moment(row, side, h, t, 1)
                        a = t-(n-dimension)
                        for r in range(1, h+1):
                            moment(row, 1-side, r, n-t, -split_coefficient(a, h, r))
                        constant = split_coefficient(a, h, 0)*comb(n, n-t)
                        self.equalities.append(normalize(row, constant))
            if spectra is not None:
                if len(spectra) != 2 or len(spectra[side]) != n+1:
                    raise ValueError('two matching binary spectrum caps required')
                for u in range(1, last+1):
                    row = {}
                    add_term(row, (side, 1, u), 1)
                    add_term(row, (side, 1, u-1), -1)
                    if row:
                        self.inequalities.append(normalize(row, int(spectra[side][u])))

    def verify_witness(self, witness):
        """Apply saved multipliers to fresh constraints, never trust a saved cap."""
        rank, support = witness['rank'], witness['support']
        if (type(rank) is not int or not 1 <= rank <= self.g
                or type(support) is not int or not 0 <= support <= self.last
                or witness['inequality_count'] != len(self.inequalities)
                or witness['equality_count'] != len(self.equalities)):
            raise ValueError('witness dimensions do not match regenerated model')
        def expand(entries, size):
            result = [Q(0)]*size
            seen = set()
            for i, value in entries:
                if type(i) is not int or not 0 <= i < size or i in seen:
                    raise ValueError('distinct valid multiplier indices required')
                seen.add(i)
                result[i] = Q(value)
            return result
        y = expand(witness['y'], len(self.inequalities))
        z = expand(witness['z'], len(self.equalities))
        if any(v < 0 for v in y):
            raise ValueError('inequality multipliers must be nonnegative')
        key = (0, rank, support)
        if key not in self.indices:
            return 0
        upper, _ = verify_dual(self.inequalities, self.equalities, {self.indices[key]: Q(1)},
                              y, z, len(self.keys))
        if upper < 0:
            raise ArithmeticError('saved dual contradicts the counting premises')
        upper = min(Q(1), upper)
        return self.caps[0][rank-1][support]*upper.numerator//upper.denominator

    def solve(self, rank, support, repair_rounds=0, time_limit=5.):
        if (type(rank) is not int or not 1 <= rank <= self.g
                or type(support) is not int or not 0 <= support <= self.last
                or type(repair_rounds) is not int or not 0 <= repair_rounds <= 3
                or not isfinite(time_limit) or not 0 < time_limit <= 60):
            raise ValueError('valid rank and support within the modeled range required')
        key = (0, rank, support)
        if key not in self.indices:
            return dict(rank=rank, support=support, cap=0, proposal=0., correction='0', status='zero prior cap')
        size = len(self.keys)
        objective = {self.indices[key]: Q(1)}
        target = np.zeros(size)
        target[self.indices[key]] = -1.
        def matrix(rows):
            values, ri, ci = [], [], []
            for j, (row, _) in enumerate(rows):
                for i, value in row.items():
                    values.append(float(value)); ri.append(j); ci.append(i)
            return csr_matrix((values, (ri, ci)), shape=(len(rows), size))
        if not hasattr(self, '_arrays'):
            self._arrays = dict(A_ub=matrix(self.inequalities),
                b_ub=np.array([float(b) for _, b in self.inequalities]),
                A_eq=matrix(self.equalities) if self.equalities else None,
                b_eq=np.array([float(b) for _, b in self.equalities]) if self.equalities else None)
        prior = self.caps[0][rank-1][support]
        def dyadic(value):
            if not isfinite(value):
                raise ArithmeticError('nonfinite LP dual proposal')
            return Q(round(float(value)*(1 << 40)), 1 << 40)
        def decode(fit, name):
            multipliers = [max(Q(0), dyadic(-v)) for v in fit.ineqlin.marginals]
            y = multipliers[:len(self.inequalities)]
            if name in ('relaxed-equality-proposal', 'coarse-safe-proposal'):
                q = len(self.inequalities)
                count = len(self.equalities)
                z = [multipliers[q+i]-multipliers[q+count+i] for i in range(count)]
            else:
                z = [dyadic(-v) for v in fit.eqlin.marginals]
            return y, z
        best = dict(rank=rank, support=support, cap=prior, prior=prior, proposal=None,
                    status='retained prior bound', attempts=[])
        # A stricter float model can cycle or become numerically infeasible.
        # Retain both independently checked candidates, with bounded effort.
        for name, options in (('default', {'time_limit': time_limit}),
                              ('small-coefficients', {'time_limit': time_limit,
                               'dual_feasibility_tolerance': 1e-9, 'primal_feasibility_tolerance': 1e-9,
                               'small_matrix_value': 1e-12}),
                              ('relaxed-equality-proposal', {'time_limit': time_limit}),
                              ('coarse-safe-proposal', {'time_limit': time_limit})):
            if name == 'relaxed-equality-proposal':
                if not self.equalities or best['cap']*8 <= prior:
                    continue
                if not hasattr(self, '_relaxed_arrays'):
                    a = self._arrays
                    self._relaxed_arrays = dict(A_ub=vstack((a['A_ub'], a['A_eq'], -a['A_eq'])),
                        b_ub=np.concatenate((a['b_ub'], a['b_eq']+1e-7, -a['b_eq']+1e-7)))
                arrays = self._relaxed_arrays
            elif name == 'coarse-safe-proposal':
                if best['cap']*8 <= prior:
                    continue
                if not hasattr(self, '_coarse_arrays'):
                    a = self._arrays
                    if self.equalities:
                        matrix = vstack((a['A_ub'], a['A_eq'], -a['A_eq'])).tocsr()
                        rhs = np.concatenate((a['b_ub'], a['b_eq'], -a['b_eq']))
                    else:
                        matrix = a['A_ub'].copy()
                        rhs = a['b_ub'].copy()
                    # The solver discards small matrix entries anyway.
                    # Explicitly budget removed negative coefficients
                    # using 0<=x<=1 before asking for multipliers. This
                    # float relaxation is NOT used by exact verification.
                    tiny = np.abs(matrix.data) < 1e-8
                    negative = matrix.copy()
                    negative.data = np.where(tiny, np.minimum(matrix.data, 0.), 0.)
                    rhs = rhs-np.asarray(negative.sum(axis=1)).ravel()+1e-7
                    matrix.data[tiny] = 0.
                    matrix.eliminate_zeros()
                    self._coarse_arrays = dict(A_ub=matrix, b_ub=rhs)
                arrays = self._coarse_arrays
            else:
                arrays = self._arrays
            fit = linprog(target, **arrays, bounds=(0., 1.), method='highs', options=options)
            attempt = dict(solver=name, status=str(fit.message))
            best['attempts'].append(attempt)
            if not fit.success:
                continue
            y, z = decode(fit, name)
            # Even for the relaxed search, verification uses the ORIGINAL
            # exact equalities. The +/-1e-7 bands are never proof premises.
            upper, correction = verify_dual(self.inequalities, self.equalities, objective, y, z, size)
            repairs = []
            for iteration in range(repair_rounds):
                # Re-optimize the exact residual at unit scale. Summing
                # its proposed multipliers with the original dual gives
                # another dual for the ORIGINAL objective and equations.
                _, residual = dual_residual(self.inequalities, self.equalities, objective, y, z, size)
                scale = max(map(abs, residual))
                if not scale or not correction:
                    break
                repair = linprog(-np.array([float(v/scale) for v in residual]),
                    **arrays, bounds=(0., 1.), method='highs', options=options)
                row = dict(round=iteration+1, status=str(repair.message))
                repairs.append(row)
                if not repair.success:
                    break
                dy, dz = decode(repair, name)
                yy = [a+scale*b for a, b in zip(y, dy)]
                zz = [a+scale*b for a, b in zip(z, dz)]
                candidate, residual_correction = verify_dual(
                    self.inequalities, self.equalities, objective, yy, zz, size)
                row.update(verified_ratio=str(candidate), correction=str(residual_correction))
                if candidate >= upper:
                    break
                y, z, upper, correction = yy, zz, candidate, residual_correction
            if upper < 0:
                raise ArithmeticError('LP dual contradicts the counting premises')
            upper = min(Q(1), upper)
            cap = prior*upper.numerator//upper.denominator
            attempt.update(proposal=-float(fit.fun), verified_ratio=str(upper), correction=str(correction), repairs=repairs)
            if cap < best['cap'] or (cap == best['cap'] and best['proposal'] is None):
                best.update(cap=cap, proposal=-float(fit.fun), verified_ratio=str(upper),
                            correction=str(correction), status='exact dual checked', solver=name,
                            active_inequalities=sum(bool(v) for v in y), active_equalities=sum(bool(v) for v in z))
                best['dual_witness'] = dict(rank=rank, support=support,
                    inequality_count=len(self.inequalities), equality_count=len(self.equalities),
                    y=[[i, str(v)] for i, v in enumerate(y) if v],
                    z=[[i, str(v)] for i, v in enumerate(z) if v])
        return best


def refine(primal, dual, dimensions, dual_dimensions, spectrum, last=144,
           supports=range(96, 137), record=None, symmetry=False, repair_rounds=0, time_limit=5., witnesses=None):
    """Tighten rank four only; all preceding valid caps remain available."""
    system = Constraints(primal, dual, dimensions, dual_dimensions, 128, last,
                         spectra=(spectrum, dual_shell_caps()), ones=(symmetry, symmetry))
    result = [row[:] for row in primal]
    multiplicity = rank_total(4, 4, 4)
    print('JOINT COUNTS', len(system.keys), 'variables', len(system.inequalities),
          'inequalities', len(system.equalities), 'equalities', flush=True)
    for u in supports:
        prior = system.caps[0][3][u]
        if witnesses is None:
            row = system.solve(4, u, repair_rounds, time_limit)
        else:
            witness = witnesses[u]
            if witness['rank'] != 4 or witness['support'] != u:
                raise ValueError('saved witness targets a different rank or support')
            row = dict(rank=4, support=u, cap=system.verify_witness(witness), prior=prior,
                       status='exact dual replayed', dual_witness=witness)
        gain = log2(prior)-log2(row['cap']) if row['cap'] else None
        print('JOINT COUNTS support', u, 'status', row['status'], 'proposal ratio', row.get('proposal'),
              'verified gain bits', gain, 'residual correction', float(Q(row.get('correction', '0'))), flush=True)
        result[3][u] = min(result[3][u], row['cap']*multiplicity)
        if record:
            record(dict(row, cap=str(row['cap']), prior=str(prior), gain_bits=gain))
    for u in range(255, -1, -1):
        result[3][u] = min(result[3][u], result[3][u+1])
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--last', type=int, default=144)
    parser.add_argument('--supports', type=int, nargs='+', default=[104, 108, 112, 116, 120, 128])
    parser.add_argument('--output', type=Path)
    parser.add_argument('--complement-symmetry', action='store_true')
    parser.add_argument('--repair-rounds', type=int, choices=range(4), default=0)
    parser.add_argument('--time-limit', type=float, default=5., help='Seconds per numerical proposal, at most 60; exact verification remains required')
    parser.add_argument('--replay', type=Path, help='Recheck saved dual multipliers against freshly generated constraints; no saved caps trusted')
    args = parser.parse_args()
    if not 128 <= args.last <= 256 or any(not 38 <= u <= args.last for u in args.supports):
        parser.error('support limit 128..256 and targets 38..last required')
    witnesses = None
    if args.replay:
        saved = json.loads(args.replay.read_text())
        if (saved['schema'] != 'shared-gf16-constraint-counts-1' or saved['last'] != args.last
                or saved['complement_symmetry'] != args.complement_symmetry):
            parser.error('saved model settings do not match')
        witnesses = {}
        for row in saved['rows']:
            if 'dual_witness' in row:
                if row['support'] in witnesses:
                    parser.error('duplicate saved support')
                witnesses[row['support']] = row['dual_witness']
        if any(u not in witnesses for u in args.supports):
            parser.error('every requested support needs a saved dual witness')
    _, baseline = shared_support.shared_counts(True)
    spectrum = shared_support.authenticated_caps()
    primal, dual, dims, ddims = refine_pair(baseline, spectrum)
    rows = []
    def record(row):
        rows.append(row)
        if args.output:
            args.output.parent.mkdir(parents=True, exist_ok=True)
            args.output.write_text(json.dumps(dict(schema='shared-gf16-constraint-counts-1',
                last=args.last, complement_symmetry=args.complement_symmetry,
                repair_rounds=args.repair_rounds, time_limit=args.time_limit, rows=rows), indent=2)+'\n')
    # The regenerated BCH premises check all-ones containment and even
    # generator rows. The latter establishes all-ones containment in its dual.
    refine(primal, dual, dims, ddims, spectrum, args.last, args.supports, record,
           args.complement_symmetry, args.repair_rounds, args.time_limit, witnesses)


if __name__ == '__main__':
    main()
