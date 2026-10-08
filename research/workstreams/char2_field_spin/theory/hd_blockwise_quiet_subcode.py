"""Deterministic quiet-state distance cap for identity-feedthrough blockwise maps.

This is an all-setups obstruction, not a probabilistic screen. The outer
needs only injectivity and block length/dimension, not the MDS property.
"""

from fractions import Fraction as F
from hd_blockwise_inner import dot, inverse, multiply


def quiet_cap(message_symbols=1 << 20, outer_length=128, outer_dimension=64,
              width=64, state=8, terminal_free=True):
    if (message_symbols < 1 or not 1 <= outer_dimension <= outer_length
            or message_symbols % outer_dimension or width < 1 or state < 0):
        raise ValueError('Invalid dimensions or nonintegral outer row count')
    rows = message_symbols//outer_dimension
    output = rows*outer_length
    if output % width:
        raise ValueError('Output must contain a whole number of packets')
    packets = output//width
    equations = state*(packets-int(terminal_free))
    active = equations//outer_dimension+1
    available = active <= rows
    return dict(message_symbols=message_symbols, outer=[outer_length, outer_dimension],
                rows=rows, output_symbols=output, width=width, state=state, packets=packets,
                terminal_free=terminal_free, feedback_equations=equations,
                selected_rows=active, guaranteed_nonzero=available,
                quiet_dimension_lower=active*outer_dimension-equations if available else None,
                output_weight_upper=active*outer_length if available else None,
                relative_distance_cap=str(F(active, rows)) if available else None,
                rejects_40_percent=available and F(active, rows) < F(2, 5),
                rejects_45_percent=available and F(active, rows) < F(9, 20),
                deterministic_every_setup=True, uses_CA_zero=False,
                requires_random_routing=False, requires_MDS=False)


def minimum_unexcluded_state(outer_length=128, outer_dimension=64, width=64,
                             distance=F(2, 5), message_symbols=1 << 20):
    if not 0 < distance <= 1:
        raise ValueError('Invalid target distance')
    for state in range(width+1):
        receipt = quiet_cap(message_symbols, outer_length, outer_dimension, width, state)
        if not receipt['guaranteed_nonzero'] or F(receipt['relative_distance_cap']) >= distance:
            return state
    raise ArithmeticError('No finite threshold found in a full-packet state')


def alternative_window_gate(outer_length, width, state):
    """Existing exact balanced-window gate with adjacent windows merged in pairs.

    Each even region contributes a suffix and the next region a prefix.
    Their union is one connected interval, with one state exit per pair.
    A zero lower bound means only that this particular gate does not reject.
    """
    from math import ceil
    from hd_blockwise_geometry import balanced_profile, failure_lower

    profile = balanced_profile(outer_length, width, F(2, 5))
    components = outer_length//2
    epsilon = F(1, 1 << 40)
    bound = failure_lower(profile['mean'], profile['maximum'], components*state)
    minimum = max(0, ceil((profile['mean']-epsilon*profile['maximum'])
                          /(components*(1-epsilon))))
    return dict(outer=[outer_length, outer_length//2], width=width, state=state,
                components=components, exit_constraints=components*state,
                window=profile['window'], high_window=profile['high_window'],
                high_regions=profile['high_regions'], total_weight=profile['total_window'],
                expected_shortened_dimension=float(profile['mean']),
                failure_probability_lower=str(bound),
                rejects_40_percent_40_bits=bound > epsilon,
                minimum_paired_gate_unexcluded_state=minimum,
                zero_lower_bound_is_not_a_certificate=True)


def nullspace(matrix, columns, field=4):
    """Small independent finite-field RREF oracle for the witness regression."""
    work = [list(row) for row in matrix]
    pivots = []
    for column in range(columns):
        pivot = next((r for r in range(len(pivots), len(work)) if work[r][column]), None)
        if pivot is None:
            continue
        top = len(pivots)
        work[top], work[pivot] = work[pivot], work[top]
        scale = inverse(work[top][column], field)
        work[top] = [multiply(scale, x, field) for x in work[top]]
        for r in range(len(work)):
            if r != top and work[r][column]:
                scale = work[r][column]
                work[r] = [a ^ multiply(scale, b, field) for a, b in zip(work[r], work[top])]
        pivots.append(column)
    basis = []
    for column in range(columns):
        if column not in pivots:
            vector = [0]*columns
            vector[column] = 1
            for row, pivot in enumerate(pivots):
                vector[pivot] = work[row][column]
            basis.append(tuple(vector))
    return basis


def finite_field_witness():
    """GF4 [4,2] rows; arbitrary route and time-varying C,A with no CA=0 premise."""
    field, b, k, rows, width, state = 4, 4, 2, 8, 4, 1
    selected = (0, 2, 5, 7)
    packets = b*rows//width
    variables = len(selected)*k
    # The regional permutations are fixed independently of the constraints.
    orders = [tuple((stride*j+c) % rows for j in range(rows))
              for c, stride in enumerate((1, 3, 5, 7))]
    feedback = [tuple(1+(i+2*j) % 3 for j in range(width)) for i in range(packets)]
    emission = [tuple(int(j == i % width) for j in range(width)) for i in range(packets)]
    mixers = [2+i % 2 for i in range(packets)]

    def routed_input(message):
        outer = [[0]*b for _ in range(rows)]
        for index, row in enumerate(selected):
            constant, linear = message[k*index:k*index+k]
            outer[row] = [constant ^ multiply(linear, point, field) for point in range(b)]
        return [outer[row][c] for c in range(b) for row in orders[c]]

    generators = []
    for column in range(variables):
        message = [0]*variables
        message[column] = 1
        generators.append(routed_input(message))
    constraints = [[dot(feedback[i], vector[i*width:(i+1)*width], field)
                    for vector in generators] for i in range(packets-1)]
    basis = nullspace(constraints, variables, field)
    promised = variables-state*(packets-1)
    if len(basis) < promised or promised <= 0:
        raise ArithmeticError('The quiet-subcode dimension bound failed')
    witnesses = []
    for message in basis:
        raw = routed_input(message)
        carry, outputs = 0, []
        emitted_states = []
        for i in range(packets):
            x = tuple(raw[i*width:(i+1)*width])
            emitted_states.append(carry)
            outputs.extend(value ^ multiply(a, carry, field) for value, a in zip(x, emission[i]))
            carry = multiply(mixers[i], carry, field) ^ dot(feedback[i], x, field)
        if any(emitted_states) or outputs != raw or not any(outputs):
            raise ArithmeticError('A quiet message was not a nonzero fixed word')
        if sum(x != 0 for x in outputs) > len(selected)*b:
            raise ArithmeticError('The routed support exceeded the selected rows')
        witnesses.append(dict(message=message, weight=sum(x != 0 for x in outputs),
                              terminal_state=carry, emitted_states=emitted_states))
    ca = [dot(feedback[i], emission[i], field) for i in range(packets)]
    if not any(ca):
        raise ArithmeticError('The test accidentally used only CA=0 maps')
    return dict(field=field, outer=[b, k], rows=rows, width=width, state=state,
                selected_rows=selected, constraints=len(constraints), variables=variables,
                nullity=len(basis), promised_nullity=promised,
                CA_values=ca, uses_time_varying_maps=True, witnesses=witnesses)


def main():
    import json
    for length in (128, 256):
        for state in (4, 8, 12, 13, 15, 16):
            print(json.dumps(quiet_cap(outer_length=length, outer_dimension=length//2, state=state)))
    for length, width, state in ((128, 32, 8), (256, 16, 4)):
        print(json.dumps(quiet_cap(outer_length=length, outer_dimension=length//2,
                                   width=width, state=state)))
        print(json.dumps(alternative_window_gate(length, width, state)))
    print(json.dumps(finite_field_witness()))


if __name__ == '__main__':
    main()
