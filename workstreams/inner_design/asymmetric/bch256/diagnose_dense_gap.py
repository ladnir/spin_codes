"""Reproduce the dense point that further geometric subdivision cannot fix."""
from fractions import Fraction as F
from flint import arb, ctx
import dense_search as search
model = search.model


def main():
    output = model.HERE/'DENSE_GAP_POINT.json'
    assert not output.exists(), 'Retained diagnostic is write-once'
    ctx.prec = 256
    checker = search.Checker(20)
    q, v = 3482, F(21,128)
    witness = checker.witness(q,q,v,v)
    log_bound = checker.bound(q,q,v,v,witness)
    # Save an integer exponent with deliberate outward slack for replay.
    bits = model.exact((log_bound/arb(2).log()).upper())
    power = -(-bits.numerator//bits.denominator)
    margin = float(-log_bound/arb(2).log())
    ctx.prec = 512
    checker = search.Checker(20)
    assert checker.bound(q,q,v,v,witness) < power*arb(2).log()
    model.base.write_new(output, dict(status='REPLAYED_DENSE_POINT_BOUND_DOES_NOT_CLOSE',
                                     message_exponent=20, occupation=q, coordinate=str(v),
                                     witness=witness, bound_upper_power=power, margin_bits_diagnostic=margin,
                                     precision_bits=256, replay_precision_bits=512,
                                     is_counterexample_to_code=False, full_distance_proved=False,
                                     source_sha256=model.sources()))
    print('Point margin', margin, 'replayed at 512 bits', flush=True)


if __name__ == '__main__':
    main()
