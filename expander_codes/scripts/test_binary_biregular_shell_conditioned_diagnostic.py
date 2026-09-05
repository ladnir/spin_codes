#!/usr/bin/env python3

from __future__ import annotations

import unittest

from binary_biregular_activation_diagnostic import exact_logterm
from binary_biregular_shell_conditioned_diagnostic import (
    activation_refined_shell_logterm,
    duration_class_refined_shell_logterm,
    duration_refined_shell_logterm,
    geometric_duration_bands,
    first_wait_refined_shell_logterm,
    shell_conditioned_logterm,
)


class ShellConditionedDiagnosticTests(unittest.TestCase):
    def test_duration_bands_cover_range_once(self) -> None:
        bands = geometric_duration_bands(20, 0.5)
        values = [value for lo, hi in bands for value in range(lo, hi + 1)]
        self.assertEqual(values, list(range(21)))

    def test_shell_conditioned_bound_dominates_exact_term(self) -> None:
        parameters = dict(
            k=15,
            left_degree=6,
            right_degree=3,
            cutoff=4,
            memory=3,
            message_weight=3,
            output_marker=0.8,
        )
        exact = exact_logterm(**parameters)
        bound, _, shells = shell_conditioned_logterm(
            **parameters,
            log_marker_min=-5.0,
            log_marker_max=2.0,
            marker_steps=29,
        )
        self.assertEqual(set(shells), {1, 3})
        self.assertGreaterEqual(bound + 1e-10, exact)

    def test_activation_refinement_dominates_exact_term(self) -> None:
        parameters = dict(
            k=15,
            left_degree=6,
            right_degree=3,
            cutoff=4,
            memory=3,
            message_weight=3,
            output_marker=0.8,
        )
        exact = exact_logterm(**parameters)
        bound, selected, _, stats = activation_refined_shell_logterm(
            **parameters,
            log_marker_min=-5.0,
            log_marker_max=2.0,
            marker_steps=15,
            refined_shells=1,
            activation_band_width=2,
            log_activation_limit=5.0,
            activation_steps=15,
            coupled_activation_input=False,
        )
        self.assertEqual(selected, [3])
        self.assertGreaterEqual(stats["improved_entries"], 0)
        self.assertGreaterEqual(bound + 1e-10, exact)

    def test_duration_refinement_dominates_exact_term(self) -> None:
        parameters = dict(
            k=15,
            left_degree=6,
            right_degree=3,
            cutoff=4,
            memory=3,
            message_weight=3,
            output_marker=0.8,
        )
        exact = exact_logterm(**parameters)
        bound, selected, bands, _, stats = duration_refined_shell_logterm(
            **parameters,
            log_marker_min=-5.0,
            log_marker_max=2.0,
            marker_steps=15,
            refined_shells=1,
            duration_relative_width=0.5,
            log_duration_limit=5.0,
            duration_steps=15,
        )
        self.assertEqual(selected, [3])
        self.assertEqual(bands[0], (0, 0))
        self.assertGreaterEqual(stats["improved_entries"], 0)
        self.assertGreaterEqual(bound + 1e-10, exact)

    def test_duration_class_refinement_dominates_exact_term(self) -> None:
        parameters = dict(
            k=15,
            left_degree=6,
            right_degree=3,
            cutoff=4,
            memory=3,
            message_weight=3,
            output_marker=0.8,
        )
        exact = exact_logterm(**parameters)
        bound, selected, _, stats = duration_class_refined_shell_logterm(
            **parameters,
            log_marker_min=-5.0,
            log_marker_max=2.0,
            marker_steps=29,
            refined_shells=1,
        )
        self.assertEqual(selected, [3])
        self.assertGreaterEqual(stats["improved_entries"], 0)
        self.assertGreaterEqual(bound + 1e-10, exact)

    def test_first_wait_refinement_dominates_exact_term(self) -> None:
        parameters = dict(
            k=15,
            left_degree=6,
            right_degree=3,
            cutoff=4,
            memory=3,
            message_weight=3,
            output_marker=0.8,
        )
        exact = exact_logterm(**parameters)
        bound, selected, cut, _, stats = first_wait_refined_shell_logterm(
            **parameters,
            log_marker_min=-5.0,
            log_marker_max=2.0,
            marker_steps=29,
            refined_shells=1,
            cut_fraction=0.5,
        )
        self.assertEqual(selected, [3])
        self.assertEqual(cut, 2)
        self.assertGreaterEqual(stats["improved_entries"], 0)
        self.assertGreaterEqual(bound + 1e-10, exact)


if __name__ == "__main__":
    unittest.main()
