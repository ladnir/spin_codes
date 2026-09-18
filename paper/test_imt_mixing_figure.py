"""Source-only rejection tests for the mixing explainer's evidence boundary."""
import copy
from pathlib import Path
import unittest

import build_imt_mixing_figure as mixing


def fixture():
    cells = []
    for m in (16,18,20):
        for r in ('1','2','3','4','8','refresh'):
            inner = dict(t=128,s=19,feedback_name='weight5_seed0',
                         transvection_rounds=None if r=='refresh' else int(r),
                         mixing_law='ideal_uniform_nonzero_marginal' if r=='refresh' else 'independent_transvection_product',
                         expansion_columns=[1]*128,feedback_columns=[2]*128)
            instance = dict(message_bits=2**m,output_bits=2**(m+1),cutoff=2**(m+1)//10,
                            outer_length=256,outer_dimension=128,distance='1/10',inner=inner)
            cells.append(dict(exponent=m,rounds=r,instance=instance,q1_upper={'numerator':'1','denominator':str(2**50)},q1_margin_bits=50.0))
    return dict(status='OUTWARD_IMT_MIXING_ROUNDS_Q1',full_distance_proved=False,cells=cells), dict(
        status='IMT_MIXING_ROUNDS_Q1_512_BIT_NATIVE_REPLAY_PASSED',producer_sha256=mixing.GRID_PIN[1],
        cells_checked=18,full_distance_proved=False)


class MixingFigureTests(unittest.TestCase):
    def test_round_plot_is_retired_from_parameter_explainer(self):
        text = Path(__file__).with_name('engineering_appendix.tex').read_text()
        explainer, remainder = text.split(r'\subsection{The certified engineering curve}', 1)
        figures = [r'\input{figures/' + name + '}' for name in
                   ('imt_parameter_k_b',
                    'imt_parameter_s_t', 'imt_parameter_k_s')]
        positions = [explainer.index(figure) for figure in figures]
        self.assertEqual(positions, sorted(positions))
        self.assertNotIn(r'\input{figures/imt_mixing_rounds}', text)
        self.assertLess(explainer.index(r'\emph{full refresh}'), positions[0])

    def test_expected_scope_and_render(self):
        grid,replay = fixture()
        rows = mixing.validate(grid,replay)
        self.assertEqual(len(rows),18)
        text = mixing.render(rows)
        self.assertIn('not full-distance certificates',text)
        self.assertIn('Ideal full refresh',text)

    def test_wrong_replay_rejected(self):
        grid,replay = fixture()
        replay['producer_sha256'] = '0'*64
        with self.assertRaises(ValueError): mixing.validate(grid,replay)

    def test_full_claim_rejected(self):
        grid,replay = fixture()
        grid['full_distance_proved'] = True
        with self.assertRaises(ValueError): mixing.validate(grid,replay)

    def test_changed_maps_rejected(self):
        grid,replay = fixture()
        grid['cells'][1]['instance']['inner']['expansion_columns'][0] = 3
        with self.assertRaises(ValueError): mixing.validate(grid,replay)

    def test_changed_margin_rejected(self):
        grid,replay = fixture()
        grid['cells'][0]['q1_margin_bits'] = 51
        with self.assertRaises(ValueError): mixing.validate(grid,replay)

    def test_missing_cell_rejected(self):
        grid,replay = fixture()
        grid['cells'].pop()
        with self.assertRaises(ValueError): mixing.validate(grid,replay)

    def test_selected_map_binding(self):
        grid,replay = fixture()
        selected = copy.deepcopy(grid['cells'][0]['instance']['inner'])
        mixing.validate(grid,replay,selected)
        selected['feedback_columns'][0] = 3
        with self.assertRaises(ValueError): mixing.validate(grid,replay,selected)


if __name__ == '__main__':
    unittest.main()
