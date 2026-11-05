#!/usr/bin/env python3
"""Bounded exact checks for a coupled scan over GF(16).

The update is Y = X + D Q, Q_next = rotate(Q) + Y. The field is
GF(2)[t]/(t^4+t+1), theta = t, and each diagonal entry is theta or theta+1.
Rotation moves the last state coordinate to the first position.

These checks establish small-field algebraic identities, not a distance
certificate, a setup-failure bound, or a performance result. They use only
the standard library and do not write files. Run with python -B.
"""

from collections import Counter
from itertools import permutations, product
from random import Random


FIELD_SIZE = 16
MODULUS = 0x13
THETA = 2
COEFFICIENTS = (THETA, THETA ^ 1)


def require(condition, message):
    """Keep checks enabled even if Python is invoked with -O."""
    if not condition:
        raise AssertionError(message)


def polynomial_multiply(a, b):
    result = 0
    while b:
        if b & 1:
            result ^= a
        b >>= 1
        a <<= 1
        if a & FIELD_SIZE:
            a ^= MODULUS
    return result


MUL = tuple(tuple(polynomial_multiply(a, b) for b in range(FIELD_SIZE))
            for a in range(FIELD_SIZE))


def rotate(state):
    return state[-1:] + state[:-1]


def unrotate(state):
    return state[1:] + state[:1]


def add(left, right):
    return tuple(a ^ b for a, b in zip(left, right))


def diagonal(coefficients, state):
    return tuple(MUL[a][q] for a, q in zip(coefficients, state))


def step(inputs, state, coefficients):
    outputs = add(inputs, diagonal(coefficients, state))
    return outputs, add(rotate(state), outputs)


def inverse_step(outputs, next_state, coefficients):
    state = unrotate(add(next_state, outputs))
    return add(outputs, diagonal(coefficients, state)), state


def transition_matrix(coefficients):
    width = len(coefficients)
    matrix = [[0] * width for _ in range(width)]
    for row, coefficient in enumerate(coefficients):
        matrix[row][row] ^= coefficient
        matrix[row][(row - 1) % width] ^= 1
    return matrix


def determinant(matrix):
    """Independent Leibniz expansion; permutation signs vanish in char 2."""
    result = 0
    for columns in permutations(range(len(matrix))):
        term = 1
        for row, column in enumerate(columns):
            term = MUL[term][matrix[row][column]]
        result ^= term
    return result


def matrix_vector(matrix, vector):
    result = []
    for row in matrix:
        value = 0
        for a, b in zip(row, vector):
            value ^= MUL[a][b]
        result.append(value)
    return tuple(result)


def weight(vector):
    return sum(value != 0 for value in vector)


def check_field():
    for a in range(1, FIELD_SIZE):
        require(sum(MUL[a][b] == 1 for b in range(1, FIELD_SIZE)) == 1,
                ("missing or nonunique inverse", a))
    power = 1
    for _ in range(4):
        power = MUL[power][THETA]
    require(power == (THETA ^ 1), "theta does not satisfy the modulus")


def check_local_transitions():
    determinant_cases = 0
    state_cases = 0
    branch_weights = Counter()
    for width in (1, 2, 3):
        zero = (0,) * width
        states = tuple(product(range(FIELD_SIZE), repeat=width))
        for coefficients in product(COEFFICIENTS, repeat=width):
            matrix = transition_matrix(coefficients)
            coefficient_product = 1
            for a in coefficients:
                coefficient_product = MUL[coefficient_product][a]
            det = determinant(matrix)
            require(det == (coefficient_product ^ 1),
                    ("cyclic determinant identity", coefficients, det))
            require(det != 0, ("singular transition below extension degree", coefficients))
            determinant_cases += 1
            for state in states:
                outputs, next_state = step(zero, state, coefficients)
                require(next_state == matrix_vector(matrix, state),
                        ("matrix/update mismatch", coefficients, state))
                require(weight(outputs) == weight(state),
                        ("nonzero diagonal changed support", coefficients, state))
                require(any(next_state) == any(state),
                        ("empty-input state persistence", coefficients, state))

                # These inputs force Y=0. The resulting state must rotate,
                # independently of the coefficient distribution.
                killing_inputs = diagonal(coefficients, state)
                quiet_outputs, quiet_next = step(killing_inputs, state, coefficients)
                require(quiet_outputs == zero and quiet_next == rotate(state),
                        ("zero-output transition identity", coefficients, state))
                require(any(quiet_next) == any(state),
                        ("zero output extinguished nonzero state", coefficients, state))
                if width == 2 and any(state):
                    branch = weight(state) + weight(next_state)
                    require(branch >= 3, ("width-two branch bound", coefficients, state))
                    branch_weights[branch] += 1
                state_cases += 1
    require(determinant_cases == 14, "unexpected determinant case count")
    require(state_cases == 33824, "unexpected state case count")
    require(branch_weights == Counter({3: 240, 4: 780}),
            ("unexpected width-two branch distribution", branch_weights))
    return determinant_cases, state_cases, branch_weights


def check_width_four_boundary():
    # The width<degree argument is unavailable here. Check the actual family
    # rather than presuming it contains a singular matrix at width four.
    singular = 0
    for coefficients in product(COEFFICIENTS, repeat=4):
        coefficient_product = 1
        for a in coefficients:
            coefficient_product = MUL[coefficient_product][a]
        det = determinant(transition_matrix(coefficients))
        require(det == (coefficient_product ^ 1),
                ("width-four determinant identity", coefficients, det))
        singular += det == 0

    coefficients = (THETA, THETA ^ 1, THETA, THETA ^ 1)
    cases = 0
    for state in product(range(FIELD_SIZE), repeat=4):
        outputs, next_state = step(diagonal(coefficients, state), state, coefficients)
        require(outputs == (0, 0, 0, 0) and next_state == rotate(state),
                ("width-four zero-output identity", state))
        require(any(next_state) == any(state),
                ("width-four zero output extinguished nonzero state", state))
        cases += 1
    require(cases == 65536, "unexpected width-four state count")

    # A negative control OUTSIDE {theta,theta+1}: nonzero diagonal entries
    # alone do not imply an invertible empty-input transition.
    bad_coefficients = (1, 1, 1, 1)
    bad_state = (1, 1, 1, 1)
    outputs, next_state = step((0, 0, 0, 0), bad_state, bad_coefficients)
    require(determinant(transition_matrix(bad_coefficients)) == 0,
            "negative control was not singular")
    require(outputs == bad_state and next_state == (0, 0, 0, 0),
            "negative-control state was not extinguished")
    return singular, cases


def check_one_step_inverse():
    vectors = tuple(product(range(FIELD_SIZE), repeat=2))
    cases = 0
    for coefficients in product(COEFFICIENTS, repeat=2):
        for state in vectors:
            for inputs in vectors:
                outputs, next_state = step(inputs, state, coefficients)
                restored_input, restored_state = inverse_step(outputs, next_state, coefficients)
                require(restored_input == inputs and restored_state == state,
                        ("one-step inverse", coefficients, inputs, state))
                require(add(outputs, diagonal(coefficients, state)) == inputs,
                        "known-state input recovery")
                cases += 1
    require(cases == 262144, "unexpected one-step inverse count")
    return cases


def encode(inputs, coefficients):
    state = (0,) * len(inputs[0])
    outputs = []
    for packet, keys in zip(inputs, coefficients):
        packet_output, state = step(packet, state, keys)
        outputs.append(packet_output)
    return tuple(outputs), state


def decode(outputs, coefficients):
    state = (0,) * len(outputs[0])
    inputs = []
    for packet_output, keys in zip(outputs, coefficients):
        inputs.append(add(packet_output, diagonal(keys, state)))
        state = add(rotate(state), packet_output)
    return tuple(inputs), state


def check_full_inverse():
    cases = 0
    for width in (1, 2, 3, 4):
        for length in (1, 2, 7, 19):
            for fixture in range(8):
                rng = Random(100000 * width + 100 * length + fixture)
                coefficients = tuple(tuple(rng.choice(COEFFICIENTS) for _ in range(width))
                                     for _ in range(length))
                inputs = tuple(tuple(rng.randrange(FIELD_SIZE) for _ in range(width))
                               for _ in range(length))
                if fixture == 0:
                    inputs = ((0,) * width,) * length
                outputs, final_state = encode(inputs, coefficients)
                restored, inverse_final_state = decode(outputs, coefficients)
                require(restored == inputs and inverse_final_state == final_state,
                        ("full encoder inverse", width, length, fixture))
                cases += 1
    require(cases == 128, "unexpected full encoder case count")
    return cases


def main():
    check_field()
    determinants, states, branches = check_local_transitions()
    singular_four, states_four = check_width_four_boundary()
    one_step = check_one_step_inverse()
    complete = check_full_inverse()
    print("PASS: GF(16) local algebra; not a distance certificate or timing result")
    print(f"widths_1_to_3_determinants={determinants}; state_cases={states}")
    print(f"width_2_branch_sum_3={branches[3]}; branch_sum_4={branches[4]}")
    print(f"width_4_family_determinants=16; singular={singular_four}; "
          f"zero_output_state_cases={states_four}; outside_family_negative_controls=1")
    print(f"one_step_inverse_cases={one_step}; full_encoder_inverse_cases={complete}")


if __name__ == "__main__":
    main()
