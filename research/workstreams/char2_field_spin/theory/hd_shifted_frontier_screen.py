"""Reusable FLOAT witnesses for adaptive shifted-pair frontier probes.

These are isolated support-class screens, not a covering argument. Search
the inexpensive global coefficient bound, then replay the good-pair alias
bound at the selected witness. No benchmark or encoder change is involved.
"""

import argparse
import json

from hd_shifted_general_screen import proposal
from hd_shifted_good_pair_screen import bound


SAVED_WITNESSES = {
    (600, 424, 400, 380): dict(
        mark=.48346682697675003, alpha=.10354771628792783,
        marker_log_odds=[[-3.6670491924654645, -6.711178743827144, -2.304473696143282],
                         [-5.828374644650113, -3.521695518880656, -2.37029381976133]],
        grid64_log2_upper=-392.5101139311737,
        grid128_log2_upper=-557.3664551564143),
    (128, 128, 100, 100): dict(
        mark=.8307402879201827, alpha=.09711434764721995,
        marker_log_odds=[[-6.7061504658983235, -6.70095808586578, -4.519489399636319]]*2,
        grid64_log2_upper=-591.11796782702),
    (256, 256, 200, 200): dict(
        mark=.6918184157995838, alpha=.0979611622417065,
        marker_log_odds=[[-5.63724774919665, -5.637506863945427, -3.4580696611330435]]*2,
        grid64_log2_upper=-801.9637074951461),
    (512, 512, 400, 400): dict(
        mark=.4824969675687487, alpha=.10257997854894543,
        marker_log_odds=[[-4.499238075179427, -4.4982772924200844, -2.309822554333331]]*2,
        grid64_log2_upper=-291.71678752766456),
    (512, 512, 480, 480): dict(
        mark=.4824124480595244, alpha=.09570346725401528,
        marker_log_odds=[[-6.372925394139932, -6.372885843993692, -2.1321024831163555]]*2,
        grid64_log2_upper=-738.5731695398863),
    (512, 512, 500, 400): dict(
        mark=.48988064416643023, alpha=.09108881602998803,
        marker_log_odds=[[-8.530540990701653, -7.9614018694529385, -2.1317049965278967],
                         [-4.562527742805248, -4.5842793897087315, -2.437284040300888]],
        grid64_log2_upper=-546.7438946485927),
    (1024, 1024, 800, 800): dict(
        mark=.23149023644229288, alpha=.11372673628845778,
        marker_log_odds=[[-3.1498249138732786, -3.1498874485983013, -.8957119339551782]]*2,
        grid64_log2_upper=-2028.525510464562),
    (1100, 948, 900, 850): dict(
        mark=.23529614812717517, alpha=.10287299081437842,
        marker_log_odds=[[-3.395931245793809, -5.457087761905962, -.8771729959512455],
                         [-4.534106362553981, -3.152560466961429, -.9600680528813901]],
        grid64_log2_upper=-1004.7103027343401),
    (768, 256, 200, 200): dict(
        mark=.49072233217575395, alpha=.12589686146591356,
        marker_log_odds=[[-2.0073702498707755, -5.2848157672910645, -3.0720522883773462],
                         [-5.2848157672910645, -2.0073702498707755, -3.0720522883773462]],
        grid64_log2_upper=-3212.0713869150204),
    (896, 128, 100, 100): dict(
        mark=.502829781606838, alpha=.1418869392585884,
        marker_log_odds=[[-1.4510232167097832, -6.153312791179192, -3.938186466169325],
                         [-6.153312791179192, -1.4510232167097832, -3.938186466169325]],
        grid64_log2_upper=-3963.5004764255136),
    (512, 512, 256, 256): dict(
        mark=.48260394840441445, alpha=.1183967192010057,
        marker_log_odds=[[-3.1714090310986833, -3.171369506925939, -2.777986789320434]]*2,
        grid64_log2_upper=-3078.295368338877),
    (512, 512, 128, 128): dict(
        mark=.48159150069182743, alpha=.13883584757334735,
        marker_log_odds=[[-2.435354240950426, -2.435320723804513, -3.583244805292005]]*2,
        grid64_log2_upper=-7664.714932078728),
    (512, 512, 450, 128): dict(
        mark=.49511403204576204, alpha=.10292298862679075,
        marker_log_odds=[[-5.771388411089243, -5.61851993349889, -2.1866881921181474],
                         [-2.629069553545417, -2.661481739071303, -3.9775854091981433]],
        grid64_log2_upper=-2357.7188473390124),
}


def screen(profile, *, iterations=240, grid=64, replay=False, progress=None):
    profile = tuple(profile)
    start = SAVED_WITNESSES.get(profile)
    if replay:
        if start is None:
            raise ValueError('No saved witness for this profile')
        witness = start
    else:
        witness = proposal(profile, iterations=iterations, exact_count=True, start=start)
    if progress:
        progress(dict(stage='global-witness', profile=list(profile), witness=witness))
    result = bound(profile, witness['mark'], witness['alpha'], witness['marker_log_odds'],
                   grid=grid, exact_count=True)
    result['global_search_iterations'] = 0 if replay else iterations
    if 'log2_upper_screen' in witness:
        result['global_proposal_log2_upper'] = witness['log2_upper_screen']
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--profile', type=int, nargs=4, default=(600, 424, 400, 380))
    parser.add_argument('--iterations', type=int, default=240)
    parser.add_argument('--grid', type=int, default=64)
    parser.add_argument('--replay', action='store_true')
    args = parser.parse_args()
    print(json.dumps(screen(args.profile, iterations=args.iterations, grid=args.grid,
                            replay=args.replay,
                            progress=lambda value: print(json.dumps(value), flush=True))), flush=True)
