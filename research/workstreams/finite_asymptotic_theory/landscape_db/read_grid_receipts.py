"""Read the authenticated finite Q1 grid without mutating its index."""
import csv

import run_complete_q1_grid as grid


def snapshot():
    existing, _, _, models, _ = grid.inventory()
    observations = dict(existing)
    specs = []
    for t in grid.STEPS:
        for s in range(t.bit_length(),21):
            directory = grid.OUT/f't{t}_s{s}'
            if not (directory/'manifest.json').exists():
                continue
            grid.verify(directory)
            with (directory/'q1.csv').open(newline='') as handle:
                for row in csv.DictReader(handle):
                    if grid.key(row) in observations:
                        raise ValueError(f'duplicate grid observation: {grid.key(row)}')
                    observations[grid.key(row)] = row
            specs.append(dict(path=(directory/'q1.csv').relative_to(grid.HERE).as_posix(),
                              base='landscape_db',
                              manifest=(directory/'manifest.json').relative_to(grid.HERE).as_posix(),
                              study=f'activation_aware_complete_grid_t{t}_s{s}_v1',
                              occupation=1,map_tag='nested-pilot',transfer_review_status='activation_aware'))
    planned = list(grid.tuples(models))
    allowed = {grid.key(r) for r in planned if r['native']}
    if not set(observations).issubset(allowed):
        raise ValueError('completed observation outside grid')
    return planned,observations,specs
