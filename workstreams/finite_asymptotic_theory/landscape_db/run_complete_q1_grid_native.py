"""Resume the frozen grid plan with the equivalent native Q1 recurrence.

The original producer remains immutable. Injecting only its coefficient
function preserves maps, witnesses, grid keys and receipt format.
"""
from pathlib import Path
import json

import activation_q1_native as native
import run_complete_q1_grid as grid


def main():
    grid.OUT.mkdir(exist_ok=True)
    engine = grid.OUT / 'native_q1_engine.json'
    encoded = json.dumps(dict(engine='activation_q1_kernel.cpp',
                              binary_sha256=grid.pilot.sha(native.LIBRARY),
                              arithmetic='MSVC /O2 /fp:precise; nearest binary64'), indent=2) + '\n'
    if engine.exists() and engine.read_text() != encoded:
        raise ValueError('native engine differs from frozen grid receipt')
    if not engine.exists():
        engine.write_text(encoded)
    original_inventory = grid.inventory

    def inventory():
        existing, maps, sources, models, fields = original_inventory()
        for path in (Path(__file__), Path(native.__file__),
                     grid.HERE/'activation_q1_kernel.cpp', grid.HERE/'build_activation_q1.ps1', engine):
            sources[path.relative_to(grid.pilot.ROOT).as_posix()] = grid.pilot.sha(path)
        return existing, maps, sources, models, fields

    kernel = native.Q1Kernel()
    grid.inventory = inventory
    grid.q1.coefficient_logs = kernel.coefficients
    grid.main()


if __name__ == '__main__':
    main()
