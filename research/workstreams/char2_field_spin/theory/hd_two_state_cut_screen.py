"""Bounded floating-point screen for mixed-interval cuts with more states.

This stationary calculation is NOT a finite-K20 failure certificate and
NOT evidence of distance when it fails to find a cut. The exact two-state
finite obstruction is separately implemented in hd_two_state_cut.py.
"""
import json


def solve(matrix):
    a = [row[:] for row in matrix]
    n = len(a)
    for column in range(n):
        pivot = max(range(column, n), key=lambda i: abs(a[i][column]))
        a[column], a[pivot] = a[pivot], a[column]
        factor = a[column][column]
        if abs(factor) < 1e-14:
            raise ArithmeticError("singular numerical system")
        for j in range(column, n+1):
            a[column][j] /= factor
        for i in range(n):
            if i != column:
                factor = a[i][column]
                if factor:
                    for j in range(column, n+1):
                        a[i][j] -= factor*a[column][j]
    return [row[n] for row in a]


def stationary_statistics(p, numerator, denominator, tails):
    """Actual policy means for maximizing G-(numerator/denominator)*C.

    An interval costs `tails` endpoint dimensions. The score difference
    is a finite Markov chain under independent Bernoulli safe bits. The
    gain/cost derivatives are retained as conditional first moments, with
    the same fixed tie rule as the finite exact interval policy.
    """
    opening = tails*(numerator+denominator)
    values, indices, pending = [], {}, [-opening]
    while pending:
        difference = pending.pop()
        if difference in indices:
            continue
        indices[difference] = len(values)
        values.append(difference)
        clipped = min(0, max(-opening, difference))
        pending.extend((clipped+denominator, clipped-numerator))
    n = len(values)
    transitions = []
    for difference in values:
        clipped = min(0, max(-opening, difference))
        transitions.append(((indices[clipped+denominator], p, 1, 0),
                            (indices[clipped-numerator], 1-p, 0, 1)))
    matrix = [[0.0]*(n+1) for _ in range(n)]
    for i in range(n):
        matrix[i][i] = -1
        for j, probability, _, _ in transitions[i]:
            matrix[j][i] += probability
    matrix[-1] = [1.0]*(n+1)
    stationary = solve(matrix)

    def moment(cost):
        matrix = [[0.0]*(n+1) for _ in range(n)]
        for i, difference in enumerate(values):
            matrix[i][i] = -1
            keep = -opening < difference <= 0
            for j, probability, gain_increment, cost_increment in transitions[i]:
                if keep:
                    matrix[j][i] += probability
                increment = cost_increment if cost else gain_increment
                if difference <= -opening:
                    increment += tails if cost else -tails
                matrix[j][-1] -= probability*stationary[i]*increment
        differences = solve(matrix)
        return sum(value for value, difference in zip(differences, values) if difference > 0)

    gain, cost = moment(False), moment(True)
    return dict(gain=gain, cost=cost, mark=numerator/denominator,
                score=gain-numerator/denominator*cost)


def bounded_screen(tails, density=0.05):
    """Mix two policies on disjoint large blocks to approach G=density/2.

    The displayed support threshold is limiting and numerical. Strict gain
    slack and finite block-boundary losses have not been certified here.
    """
    target = density/2
    points = [dict(gain=0.0, cost=0.0, mark=None)]
    points += [stationary_statistics(density, 1, denominator, tails)
               for denominator in (17, 18, 19)]
    points.append(dict(gain=density, cost=1-density, mark=0.0))
    best = None
    for low in points:
        for high in points:
            if low["gain"] <= target <= high["gain"] and low["gain"] < high["gain"]:
                high_fraction = (target-low["gain"])/(high["gain"]-low["gain"])
                threshold = target+(1-high_fraction)*low["cost"]+high_fraction*high["cost"]
                if best is None or threshold < best["limiting_threshold"]:
                    best = dict(limiting_threshold=threshold,
                                high_policy_fraction=high_fraction,
                                low_policy=low, high_policy=high)
    return dict(tail_dimensions=tails, safe_density=density, **best,
                finite_failure_certificate=False,
                conclusion="numerical stationary mixed-interval screen only")


if __name__ == "__main__":
    print(json.dumps([bounded_screen(tails) for tails in (2, 3, 4, 5, 8)], indent=2))
