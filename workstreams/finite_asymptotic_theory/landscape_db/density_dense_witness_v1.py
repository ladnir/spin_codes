"""Dense witness search using character bounds for activation-state density."""
import hashlib
import json
import math
from pathlib import Path

import numpy as np

import close_bch_dense_v2 as cover
import joint_dense_witness_v2 as joint
import syndrome_density_v1 as density

HERE = Path(__file__).resolve().parent


class DensityRefiner(joint.JointRefiner):
    def __init__(self, counts, block, t, s, ac, kernel, length, bands):
        super().__init__(counts, block, t, s, ac, kernel, length, bands)
        self.epochs = density.Epochs(t, s, ac, kernel)
        self.epoch_cache.clear()
        identity = hashlib.sha256(json.dumps([self.size, t, s, ac, kernel,
                   hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                   hashlib.sha256(Path(density.__file__).read_bytes()).hexdigest()], sort_keys=True).encode()).hexdigest()
        directory = HERE/'dense_density_tables_v1'; directory.mkdir(exist_ok=True)
        path = directory/f'{identity}.json'
        if path.exists():
            payload = json.loads(path.read_text())
            if payload['identity'] != identity:
                raise ValueError('density table identity changed')
            self.table = np.array(payload['values'])
        else:
            theta = 1/(1+np.exp(-self.ygrid)); values = []
            for ordinal, x in enumerate(self.xgrid):
                epoch = self.epochs.at(math.exp(x))
                values.append(density.base.terminal_logs(density.base.epoch_mixture(epoch, t, theta), self.size//t)/self.size)
                if ordinal % 50 == 0:
                    print(f'density moment table {ordinal+1}/{len(self.xgrid)}', flush=True)
            self.table = np.array(values)
            path.write_text(json.dumps(dict(identity=identity, values=self.table.tolist()))+'\n')


class DensityCoverRefiner(cover.CoverRefiner, DensityRefiner):
    """Warm starts, joint search, and the density-aware epoch transfer."""
