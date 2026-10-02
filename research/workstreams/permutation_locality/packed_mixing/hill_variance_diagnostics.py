"""Floating attribution of a regional proposal; never certificate evidence.

The per-part contributions reproduce regional_count.propose's calculation.
Removing a part is an intentionally incomplete counterfactual, not a valid
bound. Endpoint spans describe sensitivity of the saved affine outer dual;
they are not guaranteed improvements obtainable by refining its interval.
"""
from fractions import Fraction as Q
import math


def _logsumexp(values):
    values = list(values)
    if not values:
        return -math.inf
    peak = max(values)
    if peak == -math.inf:
        return peak
    return peak + math.log(sum(math.exp(value-peak) for value in values))


def summarize(parts, outer_logs, inner_logs, common_log, groups):
    """Pure floating bookkeeping, with exact rational interval descriptors."""
    if (not parts or len(parts) != len(outer_logs) or len(parts) != len(inner_logs)
            or type(groups) is not int or groups < 1
            or not math.isfinite(common_log)
            or not all(math.isfinite(x) or x == -math.inf for x in [*outer_logs, *inner_logs])):
        raise ValueError('matching nonempty contributions, finite offset and positive group count required')
    terms = [outer+inner for outer, inner in zip(outer_logs, inner_logs)]
    normalizer = _logsumexp(terms)
    if not math.isfinite(normalizer):
        raise ValueError('at least one finite contribution required for attribution')
    log_two = math.log(2)
    rows = []
    for index, (((lo, hi), dual), outer, inner, term) in enumerate(
            zip(parts, outer_logs, inner_logs, terms)):
        lo, hi = Q(lo), Q(hi)
        if not 0 <= lo < hi <= Q(1, 4) or len(dual) != 3:
            raise ValueError('valid variance interval and three outer dual coefficients required')
        eta, mu, gamma = map(Q, dual)
        remaining = _logsumexp(terms[:index]+terms[index+1:])
        rows.append(dict(
            index=index, interval=[str(lo), str(hi)],
            dual=list(map(str, (eta, mu, gamma))),
            outer_log2=outer/log_two, inner_log2=inner/log_two,
            log2_contribution=(common_log+term)/log_two,
            fraction_of_proposal=math.exp(term-normalizer),
            outer_endpoint_span_bits=float(groups*abs(gamma)*(hi-lo))/log_two,
            without_part_log2=(common_log+remaining)/log_two if math.isfinite(remaining) else None,
            removal_gain_bits=(normalizer-remaining)/log_two if math.isfinite(remaining) else None))
    dominant = max(range(len(rows)), key=lambda i: terms[i])
    return dict(
        diagnostic_only=True, certificate=False,
        log2_proposal=(common_log+normalizer)/log_two,
        common_log2=common_log/log_two,
        dominant_part=dominant,
        descending_parts=sorted(range(len(rows)), key=lambda i: terms[i], reverse=True),
        parts=rows,
        interpretation='Floating attribution only. Removing an interval leaves a gap; '
            'endpoint spans are saved-dual sensitivities, not certified or predicted gains.')


def diagnose(model, cell, witness):
    """Reproduce one regional proposal and expose each variance contribution.

    This performs the same local/placement work as regional_count.propose;
    call from a selected probe, not as a cheap receipt-only parser. Existing
    in-memory proposal scratch is reused under the original checked keys.
    No saved numerical endpoint or proof-status field is consumed.
    """
    import numpy as np
    import regional_count as regional
    import regional_count_probe as floating
    import scalar_cover as sc

    parts, checked = regional.prepare_witness(model, cell, witness)
    lam = Q(witness['parameters'][0])
    scratch = regional.proposal_cache(model)
    local = regional.local_operators(model, witness, _proposal_scratch=scratch)
    size = local[0].nrows()
    arrays = np.array([[[float(matrix[i, j]) for j in range(size)]
                       for i in range(size)] for matrix in local])
    region, region_logs = floating.scaled_placement(arrays, 64)
    scale, count_weights = regional.proposal_count_weights(
        model, cell, witness, parts, checked, scratch)
    _, _, logs = model.family(model.tilt)
    features = np.array(list(map(float, model.features)))
    active = np.array(model.active)
    outer_logs, inner_logs = [], []
    for (interval, dual), log_weights in zip(parts, count_weights):
        matrix, shift = floating.weighted_region(region, log_weights+region_logs)
        inner_logs.append(float(sc.log_power(matrix, sc.REGIONS)+sc.REGIONS*shift))
        eta, mu, gamma = map(float, dual)
        count = sc.G*_logsumexp(logs+eta*features+mu*active+gamma*features*(1-features))
        count -= sc.G*(min(eta*float(x) for x in cell)
                       + min(gamma*float(v) for v in interval))+mu*model.q_min
        outer_logs.append(float(count))
    common = sc.PACKETS*sc.logq(scale)+float(lam)*model.threshold
    result = summarize(parts, outer_logs, inner_logs, float(common), sc.G)
    result.update(cell=list(map(str, cell)), output_tilt=str(lam),
                  threshold=model.threshold, minimum_groups=model.q_min)
    return result
