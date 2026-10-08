"""Proof-side GF16/GF65536 sum-rank inner; not a native implementation.

Eight local Cauchy tiles share their eight-point evaluation circuit. A
dimension-four linearized Reed--Solomon generator separates their residual
subspaces. The field representation matches the frozen native candidate.
"""

from functools import lru_cache

from hd_blockwise_shared_transfer import dot, mm, multiply, power, rank, rref

Q = 65536


def transpose(matrix):
    return tuple(zip(*matrix))


def inverse(matrix):
    n = len(matrix)
    augmented = tuple(tuple(row)+tuple(int(i == j) for j in range(n))
                      for i, row in enumerate(matrix))
    reduced, pivots = rref(augmented, Q)
    if pivots[:n] != tuple(range(n)):
        raise ValueError('Singular matrix')
    return tuple(row[n:] for row in reduced)


@lru_cache(None)
def construction(state=4,tower_basis=False,fft_weights=False):
    if not 1 <= state <= 32:
        raise ValueError('LRS dimension must lie in 1..32')
    g = power(3, 17, Q)
    basis16 = tuple(power(g, j, Q) for j in range(4))
    subfield = tuple(_xor(basis16[j] for j in range(4) if i >> j & 1)
                     for i in range(16))
    assert len(set(subfield)) == 16
    assert all(power(x, 16, Q) == x for x in subfield)
    inputs, outputs = subfield[:8], subfield[8:]
    ein = tuple(tuple(power(x, j, Q) for j in range(8)) for x in inputs)
    eout = tuple(tuple(power(x, j, Q) for j in range(8)) for x in outputs)
    ei = inverse(ein)
    local_f = mm(eout, ei, Q)
    local_a = tuple(row[:4] for row in eout)
    local_c = ei[4:]
    if tower_basis:
        # The normalized additive-polynomial basis is graded, so its low
        # four polynomials span exactly the degree-below-four subspace.
        def novel(x,j):
            value=1
            for bit in range(3):
                if j>>bit&1:
                    top=bottom=1
                    for h in range(1<<bit):
                        top=multiply(top,x^subfield[h],Q)
                        bottom=multiply(bottom,subfield[1<<bit]^subfield[h],Q)
                    value=multiply(value,multiply(top,power(bottom,Q-2,Q),Q),Q)
            return value
        n0=tuple(tuple(novel(x,j) for j in range(8)) for x in inputs)
        n1=tuple(tuple(novel(x,j) for j in range(8)) for x in outputs)
        local_a=tuple(row[:4] for row in n1)
        local_c=inverse(n0)[4:]
        v=256;u=power(v,16,Q)^v
        assert u==79 and power(u,16,Q)^u==1
        basis=(1,u,v,multiply(u,v,Q))
    else:
        basis = tuple(power(256, j, Q) for j in range(4))
    moore = tuple(tuple(power(b, 16**j, Q) for b in basis) for j in range(state))
    assert rank(moore, Q) == min(4,state)
    if fft_weights and (state!=16 or not tower_basis):
        raise ValueError('FFT weights require the sixteen-state tower basis')
    representatives = outputs if fft_weights else tuple(power(g, i, Q) for i in range(8))
    assert len({power(a, 4, Q) for a in representatives}) == 8
    blocks = tuple(tuple(tuple(multiply(power(a, j, Q), moore[j][l], Q)
                               for l in range(4)) for j in range(state))
                   for a in representatives)
    if fft_weights:
        # n_g(a)^4 is the graded novel basis evaluated at a^4, on the
        # fourth-powered additive basis. This is an invertible row change
        # within each j mod 4 residue of the original Moore generator.
        blocks=tuple(tuple(tuple(multiply(multiply(power(a,j%4,Q),
                    power(novel(a,j//4),4,Q),Q),moore[j][l],Q)
                    for l in range(4)) for j in range(state))
                    for a in representatives)
    f = tuple(tuple(local_f[i % 8][j % 8] if i // 8 == j // 8 else 0
                    for j in range(64)) for i in range(64))
    c_blocks = tuple(mm(block, local_c, Q) for block in blocks)
    a_blocks = tuple(mm(local_a, transpose(block), Q) for block in blocks)
    c = tuple(tuple(x for block in c_blocks for x in block[j]) for j in range(state))
    a = tuple(row for block in a_blocks for row in block)
    return dict(subfield=subfield, inputs=inputs, outputs=outputs, basis=basis,
                moore=moore, representatives=representatives, blocks=blocks,
                local_f=local_f, local_a=local_a, local_c=local_c, F=f, A=a, C=c)


def _xor(values):
    result = 0
    for value in values:
        result ^= value
    return result


def mask_ranks(potential, output_zeros, state=4,tower_basis=False,fft_weights=False):
    """Return actual and proposed ranks for one concrete pair of masks."""
    maps = construction(state,tower_basis,fft_weights)
    f, a, c = (maps[key] for key in ('F', 'A', 'C'))
    restricted = tuple(tuple(f[i][j] for j in potential) for i in output_zeros)
    quiet = restricted+tuple(tuple(row[j] for j in potential) for row in c)
    flat = tuple(tuple(f[i][j] for j in potential)+a[i] for i in output_zeros)
    occupancies = [sum(j // 8 == i for j in potential) for i in range(8)]
    zeros = [sum(j // 8 == i for j in output_zeros) for i in range(8)]
    r0 = sum(min(x, y) for x, y in zip(occupancies, zeros))
    residual = sum(min(4,max(0, x-y)) for x, y in zip(occupancies, zeros))
    constraints = sum(min(4,max(0, y-x)) for x, y in zip(occupancies, zeros))
    return dict(quiet_actual=rank(quiet, Q), quiet_expected=r0+min(state, residual),
                flat_actual=rank(flat, Q), flat_expected=r0+min(state, constraints),
                input_actual=rank(restricted, Q), input_expected=r0)
