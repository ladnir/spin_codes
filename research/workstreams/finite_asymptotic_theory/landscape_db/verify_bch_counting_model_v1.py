"""Check fixed-calibration BCH predictions by direct 90-digit count scaling."""
import json
from pathlib import Path

import mpmath as mp

import study_bch_dominance_v1 as study

HERE = Path(__file__).resolve().parent


def main():
    path = HERE/'bch_k_counting_model_v3.json'
    model = json.loads(path.read_text())
    for name, digest in model['source_sha256'].items():
        if study.sha(Path(name)) != digest:
            raise ValueError('Stale model dependency')
    checks = []
    with mp.workdps(90):
        ln2 = mp.log(2)
        for row in model['checks']:
            b, t, s = (row[k] for k in ('block_bits', 'step_bits', 'state_bits'))
            calibration = row['calibration_exponent']; comparison = row['comparison_exponent']
            if calibration != 20 or comparison <= calibration:
                raise ValueError('Expected predictions from K=2^20 at larger K')
            anchor_path = HERE/f'bch_full_reference_b{b}_t{t}_s{s}_e{calibration}.json'
            current_path = HERE/f'bch_full_reference_b{b}_t{t}_s{s}_e{comparison}.json'
            for source in (anchor_path, current_path):
                if model['source_sha256'][str(source)] != study.sha(source):
                    raise ValueError('Prediction reference was not authenticated by the model')
            anchor = json.loads(anchor_path.read_text()); current = json.loads(current_path.read_text())
            l0 = (1 << calibration)//(b//2); length = (1 << comparison)//(b//2)

            def terms(reference):
                return [mp.exp(mp.mpf(c['log_upper'])) for c in reference['components']]

            base = terms(anchor); observed = terms(current)
            predicted = base[0]*mp.mpf(length)/l0 + base[1]*mp.binomial(length, 2)/mp.binomial(l0, 2)
            prediction = -mp.log(predicted)/ln2
            two = -mp.log(observed[0]+observed[1])/ln2
            remainder = mp.log1p(sum(observed[2:])/(observed[0]+observed[1]))/ln2
            full = -mp.log(sum(observed))/ln2
            expected = dict(predicted_two_term_margin_bits=prediction,
                            observed_two_term_margin_bits=two,
                            verified_full_margin_bits=full,
                            count_model_error_bits=prediction-two,
                            q3_and_higher_penalty_bits=remainder,
                            prediction_minus_full_bits=prediction-full)
            errors = {k: float(abs(mp.mpf(row[k])-v)) for k, v in expected.items()}
            if max(errors.values()) > 2e-12:
                raise ArithmeticError('Direct count-scaling comparison failed')
            # Tiny BCH-128 corrections must remain positive and relatively accurate.
            relative = abs(mp.mpf(row['q3_and_higher_penalty_bits'])/remainder-1)
            if remainder <= 0 or relative > mp.mpf('2e-12'):
                raise ArithmeticError('Higher-occupation correction lost precision')
            checks.append(dict(block_bits=b, comparison_exponent=comparison,
                               maximum_absolute_bit_error=max(errors.values()),
                               higher_correction_relative_error=float(relative)))
    if not {(b, e) for b in (64, 128) for e in (22, 24, 26)} <= {
            (row['block_bits'], row['comparison_exponent']) for row in checks}:
        raise ValueError('Expected both blocks at K=2^22, 2^24, and 2^26')
    result = dict(status='VERIFIED_90_DIGIT_COUNT_SCALING', checks=checks,
                  source_sha256={str(p): study.sha(p) for p in (path, Path(__file__))},
                  limitations=['Checks arithmetic of selected numerical bounds.',
                               'Does not prove the model or provide outward certificates.'])
    (HERE/'bch_k_counting_model_verification_v1.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
