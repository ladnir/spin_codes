"""Check selected finite IMT manuscript claims against authenticated evidence.

This checks transcription and exact unions, not a new interval replay.
The broad IMT plots are Q1 diagnostics; the selected curve has full certificates.
"""
import json
import re
from collections import Counter
import imt_results as evidence
import build_imt_comparison
import build_imt_parameter_figures

ROOT, require = evidence.ROOT, evidence.require
PAPER = ROOT / 'paper'


def read(name):
    return (PAPER / name).read_text()


def words(columns):
    return [sum(((c >> j) & 1) << i for i, c in enumerate(columns)) for j in range(19)]


def check_inner_terminology(source):
    # Frozen evidence URLs retain their identities; their visible labels do not
    # reintroduce the preceding inner into the paper's narrative.
    narrative = re.sub(r'https?://[^\s{}]+', '', source)
    require(not re.search(r'(?:rm2|mr2)[\s_-]*sub', narrative, re.I),
            'Historical inner name in manuscript narrative')


def check():
    data = evidence.load()
    finite, appendix, implementation = map(read, (
        'finite_certificates.tex', 'finite_appendix.tex', 'implementation.tex'))
    rows = re.findall(r'^(16|18|20|22|24) & (.*?) & (.*?) & ([0-9.]+) \\\\', finite, re.M)
    require(len(rows) == 5, 'Expected five selected margins')
    for m_text, l_text, h_text, display in rows:
        m = int(m_text)
        record = data['half'][m]
        require(int(l_text.replace(r'\,', '')) == 2**(m-7), 'Wrong row count')
        require(int(h_text.replace(r'\,', '')) == 2**(m+1)//10, 'Wrong cutoff')
        require(display == f"{record['margin_bits']:.6f}", 'Wrong rounded margin')
        match = re.search(r'\(' + str(m) + r',([0-9.]+)\)', read('figures/imt_certified_curve.tex'))
        require(match and abs(float(match[1])-record['margin_bits']) < 1e-9, 'Wrong curve point')
        maximum = {16: 63, 18: 127, 20: 511, 22: 511, 24: 511}[m]
        expected = f'{m} & {maximum} & ' + ' & '.join(
            f'{v:.6f}' for v in record['component_margin_bits'])
        require(expected in appendix, 'Wrong appendix component margins or cutoff')

    shared = read('structured_imt_appendix.tex')
    maps = re.findall(r'^(\d+) & \\texttt\{([0-9a-f]+)\} & \\texttt\{([0-9a-f]+)\}', shared, re.M)
    require([int(i) for i, _, _ in maps] == list(range(19)), 'Wrong shared map indices')
    inner = data['half'][16]['instance']['inner']
    a, c = words(inner['expansion_columns']), words(inner['feedback_columns'])
    require([int(v, 16) for _, v, _ in maps] == a, 'Wrong shared A map')
    require([int(v, 16) for _, _, v in maps] == c, 'Wrong shared C map')
    spectrum, state = Counter({0: 1}), 0
    for i in range(1, 2**19):
        state ^= a[(i & -i).bit_length()-1]
        spectrum[state.bit_count()] += 1
    require(spectrum == {0: 1, 48: 5166, 56: 110288, 64: 293455, 72: 110128, 80: 5250},
            'Wrong expansion spectrum')
    require(all(column.bit_count() == 5 for column in inner['feedback_columns']),
            'Half feedback is not weight five')

    header = (evidence.BASE / 'routing_opt/deployment/AsymmetricMap.h').read_text()
    arrays = re.findall(r'array<std::uint32_t,T> columns\{([^}]+)\}', header)
    require(len(arrays) == 2, 'Ambiguous quarter map header')
    qa, qc = [[int(v, 16) for v in row.split(',')] for row in arrays]
    require(qa == inner['expansion_columns'], 'Quarter expansion differs')
    require(all(v.bit_count() == 3 for v in qc), 'Quarter feedback is not weight three')
    quarter_words = re.findall(r'^(\d+) & \\texttt\{([0-9a-f]+)\}', appendix, re.M)
    require([int(i) for i, _ in quarter_words] == list(range(19)), 'Wrong quarter map indices')
    require([int(v, 16) for _, v in quarter_words] == words(qc), 'Wrong quarter C map')
    for row, display in zip(data['quarter_proof']['results'], ('0.165', '0.190'), strict=True):
        cutoff = f"{row['bad_weight']:,}".replace(',', r'\,')
        expected = f"{display} & {cutoff} & {row['target_bits']} & {row['margin_bits_diagnostic']:.6f}"
        require(expected in finite, 'Wrong quarter threshold row')
    for key in ('parent_generator_hex', 'subcode_generator_hex'):
        require(data['quarter_proof']['outer'][key][2:] in appendix, 'Wrong quarter polynomial')

    half_times = [cell['summaries']['sparse_pages']['median_ms'] for cell in data['timing']['cells']]
    expected = r'BCH-256, $1/2$ & ' + ' & '.join(f'{v:.3f}' for v in half_times)
    require(expected in implementation, 'Wrong selected timing row')
    require(r'BCH-128, $1/4$ & -- & -- & ' +
            f"{data['quarter_timing']['median_ms']:.3f}" in implementation, 'Wrong quarter timing row')
    for row, key in ((data['timing']['cells'][-1]['summaries']['sparse_pages'], 'process_medians_ms'),
                     (data['quarter_timing'], 'medians_ms')):
        for value in (min(row[key]), max(row[key])):
            require(f'{value:.3f}' in implementation, 'Wrong timing range')
        require(f"{row['retained_setup_bytes']/2**20:.2f}" in implementation, 'Wrong setup memory')
        require('$' + str(row['workspace_bytes']//2**20) + '$' in implementation, 'Wrong workspace memory')
    require('10.110' in implementation, 'Missing current measured transpose time')
    for name in ('abstract.tex', 'introduction.tex', 'implementation.tex'):
        require('11.259' not in read(name), 'Historical headline presented as current')
    require(r'Q_{i+1}=M_iQ_i+C(X_i)' in finite, 'Wrong finite state update')
    require(r'\alpha_i' not in finite and 'four-state' not in appendix, 'Old refresh proof remains')
    require('Q1 alone does not control' in finite and
            'do not certify the different BCH-64/128 diagnostic configurations' in finite,
            'Missing Q1 scope distinction')
    for name in ('[64,32,12]', '[128,64,22]'):
        require(name in finite, 'Missing diagnostic outer identity')
    short = data['half'][16]
    q1_short = short['component_margin_bits'][0]
    require(f'${q1_short:.3f}$ bits' in finite and
            f"${q1_short-short['margin_bits']:.3f}$ bits" in finite,
            'Wrong short-length Q1 margin or higher-occupancy loss')
    require('existing spectrum-model curve' not in finite, 'Unshown historical curve referenced')
    for source in PAPER.glob('*.tex'):
        check_inner_terminology(source.read_text(encoding='utf-8'))
    require('Historical parameter study' not in finite and
            r'\input{figures/parameter_' not in finite, 'Historical figures remain active')
    parameter_check = build_imt_parameter_figures.check(data)
    for name in ('imt_parameter_k_b', 'imt_parameter_s_t', 'imt_parameter_k_s', 'imt_certified_curve'):
        require(r'\input{figures/' + name + '}' in finite, 'Current IMT figure not included')
    require(read('figures/transposed_comparison.tex') == build_imt_comparison.table(data),
            'Stale external comparison table')
    main = read('main.tex')
    require(main.index(r'\input{structured_spin}') < main.index(r'\input{finite_certificates}')
            < main.index(r'\input{scaling_complexity}'), 'Wrong section order')
    return dict(status='SELECTED_FINITE_IMT_INTEGRATION_PASSED', selected_certificates=7,
                matched_timing_cells=4, map_words_checked=57,
                authenticated_files=data['authenticated_files'], full_interval_replay=False,
                parameter_plot_migration_complete=True, parameter_study=parameter_check)


if __name__ == '__main__':
    print(json.dumps(check(), indent=2))
