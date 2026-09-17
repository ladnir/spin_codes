"""Reject mismatched length curves, extra rounds, and full-certificate claims."""
import copy
from pathlib import Path
import unittest

import build_imt_length_figure as figure


def fixture():
    maps,cells = {},[]
    for b in (64,128,256):
        for t in ((16,32,64,128) if b==256 else (8,16,32,64)):
            s = min(19 if b==256 else 20,(t.bit_length()-1)*t.bit_length()//2)
            maps[f'b{b}_t{t}'] = dict(t=t,s=s,transvection_rounds=1,spectrum={'1':2**s-1},
                expansion_columns=list(range(1,t+1)),feedback_columns=list(range(1,t+1)))
            for m in figure.exponents(b):
                for refresh in ([False,True] if t==figure.baseline_t(b) else [False]):
                    cells.append(dict(b=b,t=t,s=s,exponent=m,refresh=refresh,sharp=True,
                        q1_margin_bits=40.,outer_rows=2**(m+1)//b,epochs_per_region=2**(m+1)//b//t,
                        dominant_grid_edge=False,full_distance_proved=False))
    return dict(status='BINARY64_IMT_ADAPTIVE_LENGTH_STUDY',cells=cells,maps=maps,full_distance_proved=False),dict(
        status='VERIFIED_BINARY64_IMT_ADAPTIVE_LENGTH_STUDY',producer_sha256=figure.GRID_PIN[1],
        cells_checked=110,log_domain_cells_checked=110,maps_checked=12,largest_difference=0.,full_distance_proved=False)


class LengthFigureTests(unittest.TestCase):
    def test_render_and_selection(self):
        grid,replay = fixture()
        rows = figure.validate(grid,replay)
        self.assertEqual(figure.adaptive(rows,256,16)['t'],128)
        rows[256,32,16,False]['q1_margin_bits'] += .2
        self.assertEqual(figure.adaptive(rows,256,16)['t'],32)
        rows[256,128,16,False]['q1_margin_bits'] += .15
        self.assertEqual(figure.adaptive(rows,256,16)['t'],128)
        text = figure.render(rows)
        self.assertIn('BCH-256',text)
        self.assertIn('one transvection',text)
        self.assertIn('not full certificates',text)

    def test_full_claim_rejected(self):
        grid,replay = fixture()
        grid['cells'][0]['full_distance_proved'] = True
        with self.assertRaises(ValueError): figure.validate(grid,replay)

    def test_extra_round_rejected(self):
        grid,replay = fixture()
        grid['maps']['b256_t128']['transvection_rounds'] = 2
        with self.assertRaises(ValueError): figure.validate(grid,replay)

    def test_missing_cell_rejected(self):
        grid,replay = fixture()
        grid['cells'].pop()
        with self.assertRaises(ValueError): figure.validate(grid,replay)

    def test_baseline_map_binding(self):
        grid,replay = fixture()
        selected = copy.deepcopy(grid['maps']['b256_t128'])
        figure.validate(grid,replay,selected)
        selected['feedback_columns'][0] = 999
        with self.assertRaises(ValueError): figure.validate(grid,replay,selected)

    def test_plot_remains_in_parameter_explainer(self):
        text = Path(__file__).with_name('finite_certificates.tex').read_text()
        prefix = text.split(r'\subsection{The selected construction}')[0]
        self.assertIn(r'\input{figures/imt_parameter_k_b}',prefix)
        self.assertIn(r'p_{\mathrm{cancel}}',prefix)
        self.assertNotIn(r'\input{figures/imt_mixing_rounds}',text)


if __name__ == '__main__': unittest.main()
