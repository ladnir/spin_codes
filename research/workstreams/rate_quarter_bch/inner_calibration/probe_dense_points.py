"""Singleton diagnostics distinguish box-width losses from point-bound losses.

These points do not cover the dense domain. Unsuccessful local optimization
is not a lower bound on the best possible witness or on the true failure rate.
"""
import json
import math
from pathlib import Path

import numpy as np
import dedicated_dense as dedicated
import refine_candidates as refinement

check = dedicated.check


def main():
    folder = check.calibration.HERE
    source = folder/'t256_s18_nested_full.json'
    retained = dedicated.validate_receipt(source)
    map_path = folder/'maps/t256_s18_nested.json'
    record = json.loads(map_path.read_text())
    a = {w:n for w,n in enumerate(record['a_counts']) if w and n}
    counts = {w:n for w,n in enumerate(check.calibration.smaller_outer.spectrum()) if w and n}
    tilts = np.arange(-32,5,dtype=float)/4
    model = dedicated.LoggedDense(counts,check.B,256,18,a,record['kernel_counts'],check.L,tilts,
                                 probability_scales=(1.,),bands=[sorted(w for w in counts if w!=128),[128]])
    epochs = refinement.Epochs(record)
    points = [[check.L-q,q,0] for q in (129,256,1024,4096,8192,16384,32768)]
    points += [[0,16384,16384],[8192,8192,16384],[0,0,32768]]
    rows = []
    for row in retained['results']:
        model.cutoff = row['bad_weight']
        refiner = refinement.Refiner(row,record,epochs)
        for point in points:
            box = model.evaluate(np.array(point),np.array(point))
            initial = -box['own_log_bound']/math.log(2)
            improved = refiner.refine(box)
            output = dict(distance_target=row['distance_target'],type_counts=point,
                          grid_margin_bits=initial,refined_margin_bits=-improved['own_log_bound']/math.log(2),
                          witness=improved['witness'])
            rows.append(output)
            print(output,flush=True)
    paths = [Path(__file__),Path(dedicated.__file__),Path(refinement.__file__),Path(refinement.original.__file__),source,map_path]
    payload = dict(status='SINGLETON_DIAGNOSTICS_NOT_DOMAIN_COVERAGE',tag=record['tag'],results=rows,
                   source_sha256={p.relative_to(check.fixed.ROOT).as_posix():check.fixed.sha(p) for p in paths})
    (folder/'t256_s18_nested_point_probe.json').write_text(json.dumps(payload,indent=2)+'\n',encoding='utf-8',newline='\n')


if __name__=='__main__': main()
