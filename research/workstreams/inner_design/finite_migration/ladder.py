"""Extend frozen weight-five IMT transfers to the paper's finite size ladder."""
import argparse
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / 'asymmetric/bch256/weight5'))
import candidate
import screen
from flint import ctx

model = candidate.model
EXPONENTS = (16, 18, 20, 22, 24)


class Engine(candidate.Engine):
    def __init__(self, exponent):
        if exponent not in EXPONENTS:
            raise ValueError('Unsupported message exponent')
        # Construction of the maps does not depend on length. No numerical
        # transfer is evaluated by the frozen initializer.
        super().__init__(20, 'weight5_seed0')
        self.exponent = exponent
        self.length = 1 << (exponent - 7)
        self.output_bits = self.outer_length * self.length
        self.cutoff = self.output_bits // 10
        assert self.length % 128 == 0


class Checker(candidate.Checker):
    def __init__(self, exponent):
        # The initializer builds only length-independent outer shell costs
        # and empty instance-local caches. Replace every length scalar before
        # computing moments; inherited methods then use the new engine.
        super().__init__(20, 'weight5_seed0')
        self.engine = Engine(exponent)
        self.rows, self.cutoff = self.engine.length, self.engine.cutoff
        self.power = exponent - 6
        assert 1 << self.power == self.engine.output_bits // 128
        assert not self.epoch.cache_info().currsize
        assert not self.fixed_moment.cache_info().currsize


def q1(output, exponents, verify=False):
    saved = model.base.read(output) if verify else None
    if saved:
        model.authenticate(saved)
        assert saved['status'] == 'IMT_LADDER_Q1_ONLY'
        exponents = saved['exponents']
    elif output.exists():
        raise FileExistsError('Use a fresh output path')
    ctx.prec = 512 if verify else 256
    results = []
    for index, exponent in enumerate(exponents):
        engine = Engine(exponent)
        old = saved['results'][index] if saved else None
        if old:
            assert old['instance'] == engine.identity()
        result = dict(instance=engine.identity(),
                      q1=screen.q1(engine, old['q1'] if old else None))
        results.append(result)
        print(ctx.prec, 'Q1', exponent, result['q1']['margin_bits'], flush=True)
    if saved:
        assert results == saved['results']
        model.base.write_new(output.with_name(output.stem + '_replay.json'), dict(
            status='IMT_LADDER_Q1_512_BIT_LINEAR_REPLAY_PASSED',
            producer_sha256=model.base.sha(output), full_distance_proved=False))
    else:
        model.base.write_new(output, dict(status='IMT_LADDER_Q1_ONLY',
            exponents=list(exponents), results=results, precision_bits=256,
            full_distance_proved=False, source_sha256=candidate.sources()))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode', choices=('q1',), default='q1')
    parser.add_argument('--m', type=int, nargs='+', choices=EXPONENTS, default=[22, 24])
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    if len(set(args.m)) != len(args.m):
        parser.error('Repeated exponents')
    q1(args.output.resolve(), args.m, args.verify)
