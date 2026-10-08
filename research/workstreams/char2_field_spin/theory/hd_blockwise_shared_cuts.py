"""Source/route cuts for quiet words of the shared-transform construction.

The sparse quiet-word receipt depends only on the outer distance and the
packet feedback-kernel distance. It applies to beta8 and any beta64 redesign
that retains those parameters. It is NOT a whole-code distance certificate.
"""

from fractions import Fraction as F
from itertools import product
from math import comb


def crowded_count(positions, rows, width, threshold):
    """Count selected positions in packets with at least threshold selections."""
    if rows % width or not 1 <= threshold <= width:
        raise ValueError('Invalid packet partition or threshold')
    occupancy = [0]*(rows//width)
    for position in positions:
        if not 0 <= position < rows:
            raise ValueError('Invalid selected position')
        occupancy[position//width] += 1
    return sum(value for value in occupancy if value >= threshold)


def packet_polynomial(mark, width=64, threshold=5):
    """Coefficients of sum_j C(width,j) x^j mark^(j*[j>=threshold])."""
    mark = F(mark)
    if mark < 1 or not 1 <= threshold <= width:
        raise ValueError('Need an upper-tail mark and a feasible threshold')
    return tuple(comb(width, j)*mark**(j if j >= threshold else 0)
                 for j in range(width+1))


def regional_moment_exact(rows, active, mark, width=4, threshold=3):
    """Small exact coefficient oracle; no K20 high-dimensional DP is used."""
    if rows % width or not 0 <= active <= rows or rows > 64:
        raise ValueError('This exact oracle is restricted to at most 64 positions')
    coefficients = [F(1)]
    packet = packet_polynomial(mark, width, threshold)
    for _ in range(rows//width):
        updated = [F(0)]*min(active+1, len(coefficients)+width)
        for i, value in enumerate(coefficients):
            for j in range(min(width, active-i)+1):
                updated[i+j] += value*packet[j]
        coefficients = updated
    return coefficients[active]/comb(rows, active)


def fixed_marks(active):
    """A fixed six-piece rational schedule, not an optimization routine."""
    if not 5 <= active <= 256:
        raise ValueError('The directed receipt covers active row counts 5..256')
    mark = (44 if active <= 8 else 22 if active <= 16 else
            12 if active <= 32 else 7 if active <= 64 else
            4 if active <= 128 else 2)
    return F(active, 16384), F(mark)


def sparse_quiet_receipt(precision=160):
    """Directed union over all nonzero quiet words using at most256 source rows.

    The event is only a necessary support condition, so no field-value or
    message union is needed. A free final packet may contain up to min(a,64)
    nonzero symbols without satisfying the feedback equation.
    """
    from flint import arb, ctx

    previous = ctx.prec
    ctx.prec = precision
    try:
        rows, regions, width, outer_distance = 8192, 256, 64, 129
        packets_per_region = rows//width
        log_two = arb(2).log()
        total, worst = arb(0), (-float('inf'), None)
        for active in range(5, 257):
            x, mark = fixed_marks(active)
            polynomial = packet_polynomial(mark)
            value = sum((coefficient*x**j for j, coefficient in enumerate(polynomial)), F(0))
            g = arb(value.numerator)/arb(value.denominator)
            ax = arb(x.numerator)/arb(x.denominator)
            au = arb(mark.numerator)/arb(mark.denominator)
            target = outer_distance*active-min(active, width)
            logarithm = (regions*packets_per_region*g.log()
                         -regions*active*ax.log()
                         -(regions-1)*arb(comb(rows, active)).log()
                         -target*au.log())
            total += logarithm.exp()
            displayed = float(logarithm/log_two)
            if displayed > worst[0]:
                worst = (displayed, dict(active_rows=active, coefficient_mark=str(x),
                                        crowded_mark=str(mark), crowded_target=target,
                                        class_log2_upper=str(logarithm/log_two)))
        upper = total.log()/log_two
        return dict(outer=[256, 128], rows=rows, packet_width=width,
                    feedback_kernel_minimum=5, state=4, terminal_free=True,
                    active_row_classes=[5, 256], class_count=252,
                    at_most_four_active_rows_deterministically_impossible=True,
                    quiet_event_log2_upper=str(upper), below_2_to_minus_40=bool(upper < -40),
                    worst_class=worst[1],
                    scope='any nonzero all-emitted-states-zero word with <=256 active outer rows',
                    output_weight_unrestricted=True, feedthrough_independent=True,
                    nonquiet_words_not_covered=True, whole_code_certificate=False)
    finally:
        ctx.prec = previous


def tile_support_rank(input_support, output_support, offset=0):
    """Exact support-space and feedback ranks for one tile of the saved beta8 map."""
    from hd_blockwise_shared_maps import matrices, rank

    source, target = set(input_support), set(output_support)
    if (not source <= set(range(8)) or not target <= set(range(8))
            or offset not in range(0, 64, 8)):
        raise ValueError('Invalid tile or tile-local coordinate set')
    full, _, feedback, _, _ = matrices()
    constraints = [tuple(int(i == j) for j in range(8))
                   for i in range(8) if i not in source]
    constraints += [tuple(full[offset+i][offset+j] for j in range(8))
                    for i in range(8) if i not in target]
    check = [tuple(row[offset:offset+8]) for row in feedback]
    base = rank(constraints)
    joint = rank(constraints+check)
    return dict(input_size=len(source), output_size=len(target), offset=offset,
                support_space_dimension=8-base,
                feedback_rank_on_support_space=joint-base,
                quiet_tile_dimension=8-joint)


def quiet_tile_branch_witness():
    """A degree-three tile polynomial attains quiet input/output branch weight13."""
    from hd_blockwise_maps import mul, matvec
    from hd_blockwise_shared_maps import matrices

    full, _, feedback, _, _ = matrices()
    source = tuple(mul(mul(x, x ^ 1), x ^ 2) for x in range(8))+(0,)*56
    target = matvec(full, source)
    syndrome = matvec(feedback, source)
    if any(syndrome):
        raise ArithmeticError('A low-degree tile polynomial was not feedback-quiet')
    return dict(input_weight=sum(value != 0 for value in source),
                output_weight=sum(value != 0 for value in target),
                combined_weight=sum(value != 0 for value in source+target),
                feedback_zero=True, witness_scope='one packet only, not an outer codeword')


def base_field_quiet_moment(active=640, precision=160):
    """An all-quiet GF256 first-moment trap; no exit-only/nonquiet assumption."""
    from flint import arb, ctx

    previous = ctx.prec
    ctx.prec = precision
    try:
        rows, low, high = 8192, 3264, 3328
        selected_packets = (205*low+51*high)//64
        equations = 4*(selected_packets-1)
        dimension = 128*active-equations
        if not 1 <= active <= low or dimension < 2:
            raise ValueError('Need a feasible row set and positive quiet dimension')
        log_probability = (205*arb(comb(low, active)).log()
                           +51*arb(comb(high, active)).log()
                           -256*arb(comb(rows, active)).log())/arb(2).log()
        lower = log_probability+8*(dimension-2)
        return dict(active_rows=active, counting_field=256, quiet_dimension_lower=dimension,
                    feedback_equations=equations, selected_packets=selected_packets,
                    output_weight_upper=205*low+51*high,
                    route_event_log2=str(log_probability),
                    quiet_projective_moment_log2_lower=str(lower),
                    moment_exceeds_one=bool(lower > 0),
                    substantial_failure_probability_not_established=True)
    finally:
        ctx.prec = previous


def small_field_descent():
    """Exhaustive GF4-to-GF16 descent in a two-row, nonidentity-feedthrough toy."""
    from hd_blockwise_inner import multiply, dot
    from hd_blockwise_quiet_subcode import nullspace

    root = next(x for x in range(2, 16) if multiply(x, x, 16) ^ x ^ 1 == 0)
    embedding = (0, 1, root, root ^ 1)
    eta = next(x for x in range(16) if x not in embedding)
    components = {embedding[a] ^ multiply(eta, embedding[b], 16): (a, b)
                  for a in range(4) for b in range(4)}
    if len(components) != 16:
        raise ArithmeticError('The extension basis failed')
    for a in range(4):
        for b in range(4):
            if embedding[multiply(a, b, 4)] != multiply(embedding[a], embedding[b], 16):
                raise ArithmeticError('The subfield embedding failed')

    # Two [3,2] evaluation-code rows, followed by coordinate-major placement.
    generator = tuple(tuple(int(j == 2*r) ^ (point if j == 2*r+1 else 0)
                            for j in range(4))
                      for point in range(3) for r in range(2))
    constraints = tuple(tuple(generator[2*i][j] ^ generator[2*i+1][j]
                              for j in range(4)) for i in range(2))
    basis = nullspace(constraints, 4, 4)

    def apply(message, field):
        converted = generator if field == 4 else tuple(tuple(embedding[x] for x in row) for row in generator)
        raw = tuple(dot(row, message, field) for row in converted)
        two = 2 if field == 4 else embedding[2]
        output = tuple(value for i in range(3)
                       for value in (raw[2*i] ^ raw[2*i+1],
                                     raw[2*i] ^ multiply(two, raw[2*i+1], field)))
        return raw, output

    minima = {}
    checked = 0
    for field in (4, 16):
        best = 7
        for scalars in product(range(field), repeat=len(basis)):
            message = [0]*4
            for scalar, vector in zip(scalars, basis):
                for j, value in enumerate(vector):
                    message[j] ^= multiply(scalar, value if field == 4 else embedding[value], field)
            if not any(message):
                continue
            raw, output = apply(message, field)
            if raw[0] ^ raw[1] or raw[2] ^ raw[3]:
                raise ArithmeticError('The word was not quiet before terminal emission')
            best = min(best, sum(value != 0 for value in output))
            if field == 16:
                for component in range(2):
                    base_message = tuple(components[value][component] for value in message)
                    base_raw, base_output = apply(base_message, 4)
                    if base_raw[0] ^ base_raw[1] or base_raw[2] ^ base_raw[3]:
                        raise ArithmeticError('A base component was not quiet')
                    if base_output != tuple(components[value][component] for value in output):
                        raise ArithmeticError('Output did not commute with field descent')
                    if any(base_value and not value for base_value, value in zip(base_output, output)):
                        raise ArithmeticError('A component output exceeded the extension support')
                checked += 1
        minima[field] = best
    return dict(base_field=4, extension_field=16, subfield_embedding=embedding,
                extension_basis=(1, eta), quiet_dimension=len(basis),
                base_minimum=minima[4], extension_minimum=minima[16],
                nonzero_extension_words_checked=checked, scope='quiet words only')


def small_quiet_witnesses():
    """GF4 [4,2] outer and [4,2,3] feedback kernel, with a free terminal state."""
    from hd_blockwise_inner import dot, multiply
    from hd_blockwise_quiet_subcode import nullspace

    rows, regions, dimension, width, state, field = 4, 4, 2, 4, 2, 4
    orders = ((0, 1, 2, 3), (2, 0, 3, 1), (1, 3, 0, 2), (3, 2, 1, 0))
    feedback = ((1, 1, 1, 1), (0, 1, 2, 3))

    def encode_raw(message):
        outer = [tuple(message[2*r] ^ multiply(message[2*r+1], c, field)
                       for c in range(regions)) for r in range(rows)]
        return tuple(outer[r][c] for c in range(regions) for r in orders[c])

    generators = []
    for column in range(rows*dimension):
        message = tuple(int(j == column) for j in range(rows*dimension))
        generators.append(encode_raw(message))
    constraints = [tuple(dot(check, vector[c*width:(c+1)*width], field)
                         for vector in generators)
                   for c in range(regions-1) for check in feedback]
    basis = nullspace(constraints, rows*dimension, field)
    if field**len(basis) > 4096:
        raise ArithmeticError('The toy enumeration exceeded its intended bound')
    active_histogram = {}
    checked = 0
    for coefficients in product(range(field), repeat=len(basis)):
        message = [0]*(rows*dimension)
        for scalar, vector in zip(coefficients, basis):
            for j, value in enumerate(vector):
                message[j] ^= multiply(scalar, value, field)
        if not any(message):
            continue
        selected = {r for r in range(rows) if any(message[2*r:2*r+2])}
        active = len(selected)
        raw = encode_raw(message)
        if active <= state:
            raise ArithmeticError('A sparse quiet word contradicted the MDS support condition')
        total_crowded = 0
        for c in range(regions):
            packet = raw[c*width:(c+1)*width]
            if c < regions-1 and any(dot(check, packet, field) for check in feedback):
                raise ArithmeticError('The nullspace did not satisfy the feedback equations')
            positions = [j for j, r in enumerate(orders[c]) if r in selected]
            total_crowded += crowded_count(positions, rows, width, state+1)
        if total_crowded < 3*active-min(active, width):
            raise ArithmeticError('The crowded-position necessary condition failed')
        active_histogram[active] = active_histogram.get(active, 0)+1
        checked += 1
    return dict(field=field, outer=[regions, dimension], feedback_kernel_minimum=3,
                nullity=len(basis), nonzero_quiet_words_checked=checked,
                active_row_histogram=active_histogram, terminal_free=True)


def main():
    import json
    print(json.dumps(sparse_quiet_receipt()))
    print(json.dumps(small_quiet_witnesses()))
    print(json.dumps(quiet_tile_branch_witness()))
    print(json.dumps(base_field_quiet_moment()))
    print(json.dumps(small_field_descent()))


if __name__ == '__main__':
    main()
