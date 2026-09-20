"""Decompose retained singleton majorants; this is not a stronger bound.

The zero-state path is a summand of the nonnegative relaxed transfer. Its
size diagnoses that relaxation, not a lower bound on actual code failures.
"""
import json
import math
from pathlib import Path

import numpy as np
import probe_dense_points as probe

check = probe.check


def main():
    folder = check.calibration.HERE
    source = folder/'t256_s18_nested_point_probe.json'
    payload = probe.dedicated.validate_receipt(source)
    record = json.loads((folder/'maps/t256_s18_nested.json').read_text())
    counts = {w:n for w,n in enumerate(check.calibration.smaller_outer.spectrum()) if w and n}
    _,ps,costs = check.typed.category_parameters(counts,128,[sorted(w for w in counts if w!=128),[128]])
    epochs = probe.refinement.Epochs(record)
    rows = []
    for row in payload['results']:
        if row['refined_margin_bits']>=0: continue
        delta = check.Fraction(row['distance_target'])
        cutoff = check.N*delta.numerator//delta.denominator
        witness = row['witness']; z = witness['log_surprisal']; lam = math.exp(z)
        proposal = np.array(witness['proposal']); theta = float(proposal@ps)
        epoch = epochs.at(z)
        matrix = check.typed.dense.epoch_mixture_logs(epoch,256,np.array([theta]))
        full_moment = float(check.fixed.composition.terminal_logs(matrix,check.N//256)[0])
        zero_moment = float(matrix[0,0,0])* (check.N//256)
        point = np.array([row['type_counts']])
        def margin(moment):
            return -float(check.typed.point_logs(point,check.L,128,costs,proposal,moment,cutoff,lam)[0])/math.log(2)
        full = margin(full_moment)
        assert abs(full-row['refined_margin_bits'])<2e-6
        # Conditional degree distribution within the tilted M_00 summand.
        terms = np.array([math.log(math.comb(256,j))+j*math.log(theta)+(256-j)*math.log1p(-theta)+epoch[j,0,0]
                          for j in range(257)])
        normalizer = float(np.logaddexp.reduce(terms))
        assert abs(normalizer-matrix[0,0,0])<1e-10
        probabilities = np.exp(terms-normalizer)
        result = dict(distance_target=row['distance_target'],type_counts=row['type_counts'],
                      full_margin_bits=full,zero_state_path_majorant_margin_bits=margin(zero_moment),
                      full_minus_zero_logmoment_bits=(full_moment-zero_moment)/math.log(2),
                      zero_path_tilted_mean_epoch_weight=float(probabilities@np.arange(257)),
                      zero_path_tilted_weight_four_mass=float(probabilities[4]))
        rows.append(result); print(result,flush=True)
    paths = [Path(__file__),source]
    output = dict(status='DIAGNOSTIC_DECOMPOSITION_OF_RELAXED_BOUND',results=rows,
                  source_sha256={p.relative_to(check.fixed.ROOT).as_posix():check.fixed.sha(p) for p in paths})
    (folder/'t256_s18_nested_path_diagnostic.json').write_text(json.dumps(output,indent=2)+'\n',encoding='utf-8',newline='\n')


if __name__=='__main__': main()
