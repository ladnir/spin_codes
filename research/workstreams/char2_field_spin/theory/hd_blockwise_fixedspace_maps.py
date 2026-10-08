"""Fixed-space state maps for the beta64 shared packet.

Exact GF256 matrix receipt; coefficients extend to GF65536. This changes
the state maps, not the beta64 feedthrough, and is not a distance certificate.
"""

from functools import lru_cache

from hd_blockwise_maps import mul, power, matmul, matvec
from hd_blockwise_shared_maps import matrices, rank, identity, inverse
from hd_blockwise_shared_support import kernel


def z_value(x):
    return mul(x, x) ^ mul(64, x)


@lru_cache(None)
def fixedspace_maps():
    """Return F64,A64x4,C4x64,K4 in monomial z-state coordinates."""
    feed = matrices(64)[0]
    a = tuple(tuple(power(z_value(x), j) for j in range(4)) for x in range(64))
    first = tuple(feed[0][j] ^ int(j == 0) for j in range(4))
    g = matvec(inverse(a[:4]), first)
    k = ((g[0], g[1], g[2], g[3]),
         (g[1], 0, g[3], 0),
         (g[2], g[3], 0, 0),
         (g[3], 0, 0, 0))
    c = matmul(k, tuple(zip(*a)))
    return feed, a, c, k


def local_output_rank(support, zero_outputs, tile=0):
    """Exact rank on free supported input coordinates and shared Q4."""
    if not 0 <= tile < 8:
        raise ValueError('tile must lie in 0..7')
    sets = tuple(tuple(sorted(s)) for s in (support, zero_outputs))
    if any(len(set(s)) != len(s) or any(not 0 <= i < 8 for i in s) for s in sets):
        raise ValueError('local coordinate set required')
    support, zeros = sets
    f, a, _, _ = fixedspace_maps()
    offset = 8 * tile
    block = tuple(tuple(f[offset+i][offset+j] for j in support)+a[offset+i]
                  for i in zeros)
    formula = len(set(support) & set(zeros)) + min(4, len(set(zeros)-set(support)))
    return dict(exact_rank=rank(block), formula=formula)


def weight(values):
    return sum(value != 0 for value in values)


def singleton_cancellation(position=0):
    f, a, c, _ = fixedspace_maps()
    x = tuple(int(i == position) for i in range(64))
    q = tuple(row[position] for row in c)
    fx, aq = matvec(f, x), matvec(a, q)
    y = tuple(u ^ v for u, v in zip(fx, aq))
    tile = position // 8
    return dict(position=position, state=q,
                same_tile_output=y[8*tile:8*tile+8],
                local_joint_weight=1+weight(y[8*tile:8*tile+8]),
                packet_joint_weight=1+weight(y))


def global_quiet_not_fixed():
    """A five-symbol global-kernel word with nonzero local feedback."""
    f, _, c, _ = fixedspace_maps()
    support = (0, 8, 9, 10, 11)
    selected = tuple(tuple(row[i] for i in support) for row in c)
    basis = kernel(selected, len(support))
    assert len(basis[0]) == 1
    x = [0] * 64
    for i, row in zip(support, basis):
        x[i] = row[0]
    y = matvec(f, x)
    return dict(support=support, input=tuple(x), output=y,
                global_feedback=matvec(c, x),
                local_feedback=tuple(matvec(tuple(row[b:b+8] for row in c), x[b:b+8])
                                     for b in range(0, 64, 8)),
                input_weight=weight(x), output_weight=weight(y))


def receipt():
    f, a, c, k = fixedspace_maps()
    residual = tuple(tuple(v ^ int(i == j) for j, v in enumerate(row))
                     for i, row in enumerate(f))
    return dict(construction='beta64-fixedspace-z-monomial-state4-width64',
                k=k, distinct_z_nodes=len({z_value(i) for i in range(64)}),
                emission_rank=rank(a), feedback_rank=rank(c),
                residual_rank=rank(residual), stacked_fixed_rank=rank(c+residual),
                nonzero_state_joint_weight_lower=37,
                singleton=singleton_cancellation(),
                quiet=global_quiet_not_fixed(), whole_certificate=False)


if __name__ == '__main__':
    import json
    print(json.dumps(receipt(), indent=2))
