"""FLOAT diagnostic: retain nonuniform flat-state refresh as residual mass.

This changes the proof decomposition, not the encoder. Unlike the previous
flat-density transfer, its B->B entry contains only the uniform refresh.
The nonuniform part contributes B->R with mass epsilon V. Complete old and
new operators may be compared, but their decompositions must not be mixed
by taking unrelated entrywise minima.
"""

from fractions import Fraction as F
from functools import lru_cache
from math import comb, exp2, log2

from hd_randomized_outer_flat_moment import flat_moments
from hd_randomized_outer_joint_input import fractional_logs
from hd_randomized_outer_refined_profile import refined_profiles
from hd_randomized_outer_state_profile import power, total


@lru_cache(None)
def local_mass_transfer(log_z):
    original = refined_profiles(log_z)
    uniform = flat_moments(F(2.0**float(log_z)))
    result = []
    for j, packet in enumerate(original):
        rows = [list(row) for row in packet]
        rows[1][1] = log2(1-1/65536)+uniform[j]
        rows[1][2] = -16+uniform[j]
        result.append(tuple(tuple(row) for row in rows))
    return tuple(result)


@lru_cache(maxsize=32)
def fractional_packets(mark, alpha, mode='mass'):
    mark, alpha = F(mark), F(alpha)
    if not 0 < mark < 1 or not 0 < alpha <= 1 or mode not in ('mass', 'density'):
        raise ValueError('Invalid mark, exponent, or decomposition')
    log_z = log2(mark)
    profiles = local_mass_transfer(log_z) if mode == 'mass' else refined_profiles(log_z)
    quiet, _ = fractional_logs(mark, alpha)
    result = []
    for j, packet in enumerate(profiles):
        rows = [[float(alpha)*x for x in row] for row in packet]
        rows[0][0] = quiet[j]
        result.append(tuple(tuple(row) for row in rows))
    return tuple(result)


def screen(rows, mark, alpha, mode='mass'):
    from scipy.optimize import minimize_scalar
    mark, alpha = F(mark), F(alpha)
    packets = fractional_packets(mark, alpha, mode)
    choose = log2(comb(8192, rows))
    beta = 4096-128*log2(65535)
    prefix = choose+float(alpha)*(rows*beta-log2(65535)-838860*log2(mark))

    def objective(odds):
        matrix = tuple(tuple(total(log2(comb(64, j))+j*odds+packets[j][a][b]
                                   for j in range(65)) for b in range(3)) for a in range(3))
        run = power(matrix, 32768)
        return prefix+total(run[0])-256*(choose+rows*odds)

    selected = minimize_scalar(objective, bounds=(-20, 24), method='bounded',
                               options={'xatol':1e-7})
    return dict(active_rows=rows, mark=str(mark), alpha=str(alpha), mode=mode,
                marker_log_odds=float(selected.x), diagnostic_log2=float(selected.fun),
                directed=False, whole_code_certificate=False)


if __name__ == '__main__':
    import json
    for mode in ('density', 'mass'):
        for rows, mark, alpha in ((16,F(63,64),F(1,8)),(2048,F(1,16),F(1,16)),
                                  (3072,F(1,64),F(1,16)),(4096,F(1,256),F(1,8)),
                                  (5120,F(1,1024),F(1,8)),(6144,F(1,4096),F(1,8))):
            print(json.dumps(screen(rows,mark,alpha,mode)),flush=True)
