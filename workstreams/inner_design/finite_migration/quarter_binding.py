"""Read-only checks linking the quarter-rate IMT certificate to deployment timings."""
import argparse
import json
from pathlib import Path
import statistics
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / 'asymmetric'))
import verify_certificate as prior
from flint import ctx

ROOT = prior.ROOT
DEPLOYMENT = HERE.parent / 'routing_opt/deployment'
DATA = DEPLOYMENT.parent / 'measurements/deployment'


def replace_once(text, old, new):
    assert text.count(old) == 1, 'Ambiguous or missing reviewed source edit'
    return text.replace(old, new, 1)


def expected_source(original, supported):
    """Replay the reviewed source transformation, without changing any files."""
    code = replace_once(original, '#include "Spin.h"', '#include "Spin.h"\n#include "WorkspaceRouting.h"')
    start = supported.index('Spin::Workspace::Workspace(')
    end = supported.index('\nstd::size_t Spin::Workspace::bytes()', start)
    code = replace_once(code,
        'Spin::Workspace::Workspace(const Spin& code):buckets(code.codeBlocks()),tile(code.tileBlocks()) {}',
        supported[start:end])
    start = supported.index('    if constexpr(requires {Map::tunedWorkspaceRouting;})')
    end = supported.index('    for(std::size_t base=0;base<n;base+=tileSize)', supported.index('\n    }', start) + 6)
    marker = '    for(std::size_t base=0;base<n;base+=tileSize) {'
    code = replace_once(code, marker, supported[start:end] + marker)
    start = supported.index('        if constexpr(workspace_routing::compiled)', supported.index('void Spin::encodeUnchecked'))
    end = supported.index('        if(layout==Layout::Packed24) run<Map128S19,true,true>', start)
    dispatch = supported[start:end].replace('TunedMap<Map128S19>', 'TunedMap<AsymmetricMap>')
    marker = '        if(layout==Layout::Packed24) run<AsymmetricMap,true,true>'
    return replace_once(code, marker, dispatch + marker)


def timing_cell(runs, exponent):
    assert len(runs) == 3
    for run in runs:
        assert run['configuration'] == 't128_s19_asymmetric_greedy3_2_r1'
        assert run['m'] == exponent and run['inplace'] is True
        assert (run['outer_length'], run['outer_dimension']) == (128, 32)
        assert run['trials'] == 101 and run['layout'] == 'packed24'
        assert 0 < run['p10_ms'] <= run['median_ms'] <= run['p90_ms']
    for key in ('output_hash', 'tile_rows', 'workspace_bytes', 'retained_setup_bytes'):
        assert len({run[key] for run in runs}) == 1, key
    medians = [run['median_ms'] for run in runs]
    return dict(medians_ms=medians, median_ms=statistics.median(medians),
                **{key: runs[0][key] for key in ('output_hash', 'tile_rows', 'workspace_bytes', 'retained_setup_bytes')})


def verify():
    ctx.prec = 256
    paths = []

    def read(path):
        paths.append(path)
        return prior.read(path)

    certificate_path = prior.HERE / 'ASYMMETRIC_MARGIN_CERTIFICATE.json'
    certificate = read(certificate_path)
    replay = read(prior.HERE / 'ASYMMETRIC_MARGIN_CERTIFICATE_REPLAY.json')
    assert certificate['status'] == 'OUTWARD_ASYMMETRIC_ALL_OCCUPANCY_CERTIFICATE'
    assert certificate['inner'] == prior.certify.INNER and certificate['precision_bits'] == 256
    assert (certificate['message_bits'], certificate['output_bits']) == (1 << 20, 1 << 22)
    assert replay['status'] == 'HIGHER_PRECISION_REPLAY_PASSED' and replay['precision_bits'] == 512
    assert replay['certificate_sha256'] == prior.sha(certificate_path)
    prior.bind(certificate)
    implementation = read(prior.HERE / 'IMPLEMENTATION.json')
    prior.bind(implementation)
    verification = read(prior.HERE / 'CERTIFICATE_VERIFICATION.json')
    assert verification['status'] == 'ASYMMETRIC_CERTIFICATE_AND_IMPLEMENTATION_BOUND'
    prior.bind(verification)
    engine = prior.certify.Engine(read(prior.HERE.parent / 'NO_CONSTANT_MAP.json'))
    construction = prior.search.g.tv.smaller_outer.construction()
    assert certificate['outer'] == construction
    weights = [w for w, n in enumerate(prior.search.g.tv.smaller_outer.spectrum()) if w and n]
    witnesses = read(prior.search.g.tv.fixed.HERE / 'SMALLER_OUTWARD_WITNESSES.json')
    for key in ('q2_bands', 'adaptive_bands'):
        prior.partition(witnesses[key], weights)
    assert len(certificate['results']) == 2
    for row, delta, target, maximum in zip(certificate['results'], ('33/200', '19/100'), (40, 30), (64, 128)):
        assert row['distance_target'] == delta and row['target_bits'] == target
        assert row['maximum_sparse'] == maximum and len(row['adaptive_terms']) == maximum - 2
        d = prior.search.g.F(delta)
        assert row['bad_weight'] == (1 << 22) * d.numerator // d.denominator
        cover = read(prior.HERE / f'DENSE_greedy3_2_{delta.replace("/", "_")}.json')
        assert cover['candidate'] == 'greedy3_2'
        prior.partition(cover['bands'], weights)
        prior.search.g.check_coverage(cover['selected_boxes'], 32768, maximum + 1)
        assert len(row['dense_terms']) == len(cover['selected_boxes'])
        expected_pairs = list(prior.search.g.tv.fixed.composition.compositions(len(witnesses['q2_bands']), 2))
        assert len(row['q2_terms']) == len(expected_pairs)
        prior.check_dyadic_union([row['q1_upper']] + row['q2_terms'] + row['adaptive_terms'] + row['dense_terms'],
                                row['union_upper'], target)

    selected = next(x for x in implementation['candidates'] if x['name'] == 'asymmetric_greedy3_2_sparse')
    base = HERE.parent / 'generated' / selected['name']
    for name, digest in selected['source_sha256'].items():
        assert prior.sha(base / name) == digest
        paths.append(base / name)
    header, _ = prior.balanced.header(engine.a_columns)
    header += '\nnamespace bare_spin {\nstruct AsymmetricMap : BalancedMap {\n'
    header += 'static constexpr auto feedbackColumns=BalancedMap::columns;\n'
    header += 'static constexpr std::array<std::uint32_t,T> columns{' + ','.join(map(hex, engine.columns)) + '};\n'
    header += 'static constexpr auto groupedColumns=columns;\n'
    header += 'static constexpr std::array<unsigned,S> groupOrder{' + ','.join(map(str, range(19))) + '};\n};\n}\n'
    assert (base / 'AsymmetricMap.h').read_text() == header
    assert (DEPLOYMENT / 'AsymmetricMap.h').read_text() == header
    supported = ROOT / 'workstreams/bare_bch_rm2sub/Spin.cpp'
    assert expected_source((base / 'CandidateSpin.cpp').read_text(), supported.read_text()) == (DEPLOYMENT / 'AsymmetricSpin.cpp').read_text()
    paths.extend([supported, DEPLOYMENT / 'AsymmetricSpin.cpp', DEPLOYMENT / 'AsymmetricMap.h'])

    source_list = DATA / 'sources.sha256'
    for line in source_list.read_text().splitlines():
        digest, name = line.split(maxsplit=1)
        assert prior.sha(ROOT / name) == digest, name
        paths.append(ROOT / name)
    paths.append(source_list)
    log = DATA / 'correctness.log'
    assert '100% tests passed, 0 tests failed out of 12' in log.read_text()
    for name in ('asymmetric_off', 'asymmetric_on', 'asymmetric_fallback', 'asymmetric_portable'):
        assert f'{name}_test' in log.read_text()
    paths.append(log)
    cells = []
    for exponent in (16, 18, 20):
        states = {state: timing_cell([read(DATA / f'asymmetric_{state}-m{exponent}-r{i}.jsonl')
                                    for i in (1, 2, 3)], exponent) for state in ('off', 'on')}
        for key in ('output_hash', 'tile_rows', 'workspace_bytes', 'retained_setup_bytes'):
            assert states['off'][key] == states['on'][key]
        cells.append(dict(message_exponent=exponent, **states,
                          distance_certificate_available=exponent == 20))
    paths.append(Path(__file__).resolve())
    return dict(status='VERIFIED_QUARTER_IMT_DEPLOYMENT_BINDING',
        certificate_scope='BCH [128,32,32], K=2^20, t=128,s=19, greedy3_2, one transvection',
        certified_margins=[prior.certify.old.margin(prior.certify.old.unpack(x['union_upper'])) for x in certificate['results']],
        exact_union_and_coverage_checked=True, exact_map_header_reconstructed=True,
        reviewed_source_transformation_checked=True, retained_correctness_receipts_checked=True,
        new_benchmark_run=False, new_numerical_replay=False, performance=cells,
        timing_policy='In-place, no reset/copy; setup and workspace preparation excluded; median of three process medians, 101 trials each.',
        source_sha256={p.relative_to(ROOT).as_posix(): prior.sha(p) for p in paths})


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path)
    a = p.parse_args()
    result = verify()
    if a.output:
        with a.output.resolve().open('x', encoding='utf-8', newline='\n') as stream:
            json.dump(result, stream, indent=2)
            stream.write('\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'source_sha256'}, indent=2))
