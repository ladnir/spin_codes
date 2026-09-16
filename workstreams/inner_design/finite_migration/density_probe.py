"""Point-only test of a sharper routing comparison with unchanged IMT maps."""
import argparse
from fractions import Fraction as F
from pathlib import Path

from flint import arb, ctx
import activation_density
import density_comparison
import mixed_dense
import poisson_density_factor as old_density

model = mixed_dense.ladder.model


class Checker(mixed_dense.Checker):
    def __init__(self, exponent):
        super().__init__(exponent)
        self.engine = activation_density.Engine(exponent)

    def bound(self, lo, hi, a, b, witness):
        old = super().bound(lo, hi, a, b, witness)
        nlo, nhi = self.density(a, b)
        if nhi == 1:
            return old
        improved = density_comparison.factor(self.rows, lo, hi, nlo, nhi, self.ps)
        old_factor = arb(old_density.density_factor(self.rows))
        if improved < old_factor:
            return (old + 256 * (improved.log() - old_factor.log())).upper()
        return old


def run(output, verify=False):
    saved = model.base.read(output) if verify else None
    if saved:
        model.authenticate(saved)
        assert saved['status'] == 'IMT_ROUTING_DENSITY_POINTS_NOT_A_COVER'
    elif output.exists():
        raise FileExistsError(output)
    ctx.prec = 512 if verify else 256
    points = []
    for index, (exponent, q) in enumerate(((16, 218), (18, 870))):
        c = Checker(exponent)
        coordinate = F(21, 128)
        witness = saved['points'][index]['witness'] if saved else c.witness(q, q, coordinate, coordinate)
        upper = c.bound(q, q, coordinate, coordinate, witness)
        margin = float(-upper / arb(2).log())
        nlo, nhi = c.density(coordinate, coordinate)
        factor = density_comparison.factor(c.rows, q, q, nlo, nhi, c.ps)
        print('density point', ctx.prec, exponent, q, 'factor', float(factor), 'margin', margin, flush=True)
        bound = model.base.encode(model.exact(upper))
        if saved:
            assert upper <= model.number(model.base.decode(saved['points'][index]['log_bound']))
            assert c.engine.identity() == saved['points'][index]['instance']
        points.append(dict(instance=c.engine.identity(), q=q, coordinate=str(coordinate),
                           witness=witness, log_bound=bound, margin_bits=margin,
                           density_factor=float(factor)))
    if saved:
        model.base.write_new(output.with_name(output.stem + '_replay.json'), dict(
            status='IMT_ROUTING_DENSITY_POINTS_512_REPLAY_PASSED',
            producer_sha256=model.base.sha(output), full_distance_proved=False))
    else:
        model.base.write_new(output, dict(status='IMT_ROUTING_DENSITY_POINTS_NOT_A_COVER',
            points=points, full_distance_proved=False, source_sha256=mixed_dense.ladder.candidate.sources()))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--verify', action='store_true')
    a = p.parse_args()
    run(a.output.resolve(), a.verify)
