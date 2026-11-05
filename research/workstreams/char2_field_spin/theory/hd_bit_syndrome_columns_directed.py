"""Positive scaled column contraction, with every labelled row retained.

Inputs are exported upper bounds in hd_positive_upward.Scaled form, not
logs or raw FLOAT estimates. The five marked edges are00/atom0,00/atomB,
01,10,11. For a labelled row, gap matrices occur before, between, and
after its marked positions. The terminal coefficient depends on the full
edge-count tuple and that row's coefficient class.

The DP preserves matrix chronology. Only commuting scalar group factors
are represented by the five edge counts. All arithmetic is delegated to
the upward positive primitives: no unrounded multiply, addition or sum
occurs here. This module builds no provider or model-specific weight and
does not claim a distance certificate.
"""
import numpy as np

import hd_positive_upward as up
from hd_post_refresh_quiet_rows import transition_geometry


ENDPOINTS = ((0, 0), (0, 0), (0, 1), (1, 0), (1, 1))


def count_geometry(marked):
    """Frozen lexicographic five-count order, also allowing the empty toy."""
    if type(marked) is not int or not 0 <= marked <= 8:
        raise ValueError('Explicit marked-position count0..8 required')
    return transition_geometry(marked) if marked else (((0, 0, 0, 0, 0),), ())


def counted_columns(terminal, class_indices, edge_masses, marked,
                    gap_products, gap_codes, *, batch_size=128,
                    check_time=lambda: None):
    """Return (Scaled[row,2,2], metadata) for an ordered count DP.

    terminal[class,count] already contains the scalar group power, but no
    marked-edge mass. gap_products indexes arbitrary chronological gap
    words under the same positive measure; gap_codes[row,position] selects
    them. Their construction and theorem are the caller's responsibility.
    """
    counts, layers = count_geometry(marked)
    if not all(isinstance(value, up.Scaled) for value in (terminal, edge_masses, gap_products)):
        raise ValueError('Terminal, edge and gap values must be Scaled upper arrays')
    classes, codes = np.asarray(class_indices), np.asarray(gap_codes)
    if (terminal.mantissas.ndim != 2 or terminal.mantissas.shape[1] != len(counts)
            or not len(terminal.mantissas)
            or classes.ndim != 1 or not len(classes) or classes.dtype.kind not in 'iu'
            or np.min(classes) < 0 or np.max(classes) >= len(terminal.mantissas)
            or edge_masses.mantissas.shape != (5,)
            or gap_products.mantissas.ndim != 3 or gap_products.mantissas.shape[1:] != (2, 2)
            or not len(gap_products.mantissas)
            or codes.shape != (len(classes), marked+1) or codes.dtype.kind not in 'iu'
            or np.min(codes) < 0 or np.max(codes) >= len(gap_products.mantissas)
            or type(batch_size) is not int or not 1 <= batch_size <= 512):
        raise ValueError('Complete class/count tables, five edges and ordered gap geometry required')
    batches = []
    for start in range(0, len(classes), batch_size):
        check_time()
        stop = min(len(classes), start+batch_size)
        block = codes[start:stop]
        current = gap_products[block[:, 0]][None, ...]
        for position, (next_size, targets) in enumerate(layers):
            check_time()
            following = up.zeros((next_size, stop-start, 2, 2))
            after = gap_products[block[:, position+1]]
            for label, (incoming, outgoing) in enumerate(ENDPOINTS):
                if edge_masses.mantissas[label] == 0:
                    continue
                # The first matrix retains its initial state and the
                # last matrix retains its terminal state. The chosen
                # rank-one marked edge connects incoming to outgoing.
                term = up.multiply(current[..., incoming, None], after[None, :, None, outgoing, :])
                term = up.multiply(term, edge_masses[label])
                following = up.injective_add(following, targets[label], term)
            current = following
        # Gather and transpose change only indexing, never the exact scale.
        scalar = up.Scaled(terminal.mantissas[classes[start:stop]].T, terminal.exponent)
        weighted = up.multiply(current, scalar[:, :, None, None])
        batches.append(up.sum_axis0(weighted))
    result = up.concatenate(batches, axis=0)
    return result, dict(rowset_count=len(classes), final_atom_count_classes=len(counts),
        coefficient_classes=len(terminal.mantissas), dp_passes=1, marked_positions=marked,
        batch_size=batch_size, batches=len(batches), gap_word_slots=len(gap_products.mantissas),
        maximum_dp_layer_bytes=len(counts)*min(batch_size, len(classes))*4*8,
        edge_order=list(ENDPOINTS), count_order=counts,
        arithmetic='positive-scaled-binary64-with-primitive-upward-rounding',
        labelled_rows_preserved=True, orbit_multiplicity_factor=False,
        terminal_includes_edge_mass=False, chronological_matrices=True,
        regional_provider_built=False, whole_code_certificate=False)
