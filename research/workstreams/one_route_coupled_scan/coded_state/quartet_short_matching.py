r"""Capacity-two matching through at most one temporal arc, on one fixed window.

Each left source(t,i) has unit capacity and can match to its own input,
the four inputs of its next quartet, or the four inputs of its quartet
two updates ahead. Non-H inputs have capacity2; H inputs have capacity0.
All neighbors lie inside the window. This uses exactly the H masks and
partitions sampled by quartet_fixed_row_flow.py for the same seed.

The two-update match embeds in the full flow as state(t,i)->state(t+1,i)
followed by a next-quartet input arc. That unit temporal arc is used only
by source(t,i); all other matching routes use no temporal arc. Thus every
integral matching is a feasible full flow. Conversely, an integral maximum
full flow minimizing temporal use has no carry-through vertex with an
unused own source: replacing the incoming path prefix by that source would
reduce temporal use. At every remaining carry-through vertex, pair the
incoming carry to an input exit and its own source to the outgoing carry.
Every resulting path has at most one temporal arc. A carried path ending
at the next own input is already a next-quartet neighbor. Therefore this
matching equals the full flow for unit own-source/temporal capacities and
integral input capacities. Its deficit/2 is exactly the full physical tau.

For a left set S, define N_t=S_t union Gamma_t(S_(t-1) union S_(t-2)),
with missing earlier sets empty. The exact Hall cut is
 n*T-|S|+2|N(S)\H|.
There is no death term in this identity. The returned Hall set is a cut
witness for the restricted matching, not an actual dimension obstruction.

Reuse the compact unit-source Dinic implementation. Check the canonical
Hall neighborhood and each matched route independently. CLI runs one
32/64-update width512 window, with optional same-instance full-flow replay.
No jobs, random sampling, or other computation is performed on import.

Parent verification:55296 exhaustive Hall sets,81 exact minima,89 route
embedding checks,and12 endpoint checks passed. On the seed0 width512,T64
window, matching and full flow both equal32301: deficit467,tau233.5.
The matching uses10724 own,15701 next,and5876 two-ahead routes; the parent
reported8.846 seconds. These are structural diagnostics, not probability
or construction-level certificates.
"""

import argparse
from itertools import combinations
from operator import index
import os
import shutil
import threading
import time


SM_WIDTH, SM_REGION, SM_GLOBAL_UPDATES, SM_GLOBAL_THRESHOLD = 512, 16, 4096, 52430


def _sm_dependency(name):
    if name in globals():
        return globals()[name]
    import quartet_fixed_row_flow as module
    return getattr(module, name)


def short_matching_hall_identity(source_masks, h_masks, partitions, width):
    """Independent Hall neighborhood, including the second earlier layer."""
    width, updates = _sm_dependency('_fr_validate_instance')(h_masks, partitions, width)
    if len(source_masks) != updates:
        raise ValueError('Hall source sequence has the wrong length')
    full = (1 << width)-1
    volume = outside = total_neighbors = 0
    neighbor_masks = []
    for t, mask in enumerate(source_masks):
        if not isinstance(mask, int) or mask < 0 or mask & ~full:
            raise ValueError('Hall source mask is outside the coordinate set')
        old = (source_masks[t-1] if t else 0) | (source_masks[t-2] if t >= 2 else 0)
        neighbors = mask
        if old:
            for quartet in partitions[t]:
                quartet_mask = sum(1 << coordinate for coordinate in quartet)
                if old & quartet_mask:
                    neighbors |= quartet_mask
        volume += mask.bit_count()
        outside += (neighbors & ~h_masks[t]).bit_count()
        total_neighbors += neighbors.bit_count()
        neighbor_masks.append(neighbors)
    return dict(hall_source_size=volume, hall_noneligible_neighbors=outside,
                hall_neighbor_size=total_neighbors, neighbor_masks=tuple(neighbor_masks),
                cut_capacity=width*updates-volume+2*outside,
                integer_deficit_witness=volume-2*outside)


def _sm_validate_matching(graph, h_masks, partitions, width, slots, starts, expected):
    """Check integral assignments and their conflict-free full-flow embedding."""
    loads = [0]*slots
    temporal_sources = set()
    matched = own = next_quartet = two_ahead = 0
    to, following, capacity = graph.to, graph.following, graph.capacity
    for source in range(slots):
        edge, assigned = graph.head[source], 0
        source_time, source_coordinate = divmod(source, width)
        while edge != -1:
            target = to[edge]
            if not (edge & 1) and slots <= target < 2*slots:
                amount = capacity[edge ^ 1]
                if amount:
                    if amount != 1 or assigned:
                        raise ArithmeticError('A left vertex has a non-unit or repeated assignment')
                    assigned = 1
                    input_slot = target-slots
                    input_time, input_coordinate = divmod(input_slot, width)
                    if (h_masks[input_time] >> input_coordinate) & 1:
                        raise ArithmeticError('A matching entered an H input')
                    delay = input_time-source_time
                    if delay == 0:
                        if input_coordinate != source_coordinate:
                            raise ArithmeticError('Own-input match changed its coordinate')
                        own += 1
                    elif delay in (1, 2):
                        quartet = partitions[input_time][starts[input_time][source_coordinate]//4]
                        if input_coordinate not in quartet:
                            raise ArithmeticError('A delayed match is outside its target quartet')
                        if delay == 1:
                            next_quartet += 1
                        else:
                            # This original source alone can use this temporal arc.
                            if source in temporal_sources:
                                raise ArithmeticError('A unit temporal arc was used twice')
                            temporal_sources.add(source)
                            two_ahead += 1
                    else:
                        raise ArithmeticError('Matching route exceeds the permitted window delay')
                    loads[input_slot] += 1
                    if loads[input_slot] > 2:
                        raise ArithmeticError('A non-H input exceeds capacity2')
                    matched += 1
            edge = following[edge]
    if matched != expected or own+next_quartet+two_ahead != matched:
        raise ArithmeticError('Assignment count disagrees with matching flow')
    return dict(matched_own_input=own, matched_next_quartet=next_quartet,
                matched_two_updates_ahead=two_ahead,
                occupied_noneligible_inputs=sum(load > 0 for load in loads),
                full_flow_embedding_verified=True, maximum_temporal_arc_load=int(bool(two_ahead)))


def short_matching_solve(h_masks, partitions, width=512):
    """Solve the matching, equal to full flow, and return its exact Hall cut."""
    width, updates = _sm_dependency('_fr_validate_instance')(h_masks, partitions, width)
    slots = width*updates
    if slots > 32768:
        raise ValueError('This diagnostic is bounded to at most32768 left vertices')
    graph = _sm_dependency('_FRDinic')(2*slots+2)
    source, sink, infinity = 2*slots, 2*slots+1, slots+1
    flattened, starts = [], []
    for partition in partitions:
        flat = tuple(coordinate for quartet in partition for coordinate in quartet)
        index_by_coordinate = [0]*width
        for start in range(0, width, 4):
            index_by_coordinate[flat[start]] = start
            index_by_coordinate[flat[start+1]] = start
            index_by_coordinate[flat[start+2]] = start
            index_by_coordinate[flat[start+3]] = start
        flattened.append(flat)
        starts.append(index_by_coordinate)
    build_started = time.monotonic()
    add = graph.add_edge
    for t in range(updates):
        current, own_input_base = t*width, slots+t*width
        H = h_masks[t]
        has_next, has_two = t+1 < updates, t+2 < updates
        if has_next:
            next_input_base = own_input_base+width
            next_partition, next_starts = flattened[t+1], starts[t+1]
        if has_two:
            two_input_base = own_input_base+2*width
            two_partition, two_starts = flattened[t+2], starts[t+2]
        for coordinate in range(width):
            left = current+coordinate
            add(source, left, 1)
            add(left, own_input_base+coordinate, infinity)
            if not (H >> coordinate) & 1:
                add(own_input_base+coordinate, sink, 2)
            if has_next:
                start = next_starts[coordinate]
                add(left, next_input_base+next_partition[start], infinity)
                add(left, next_input_base+next_partition[start+1], infinity)
                add(left, next_input_base+next_partition[start+2], infinity)
                add(left, next_input_base+next_partition[start+3], infinity)
            if has_two:
                start = two_starts[coordinate]
                add(left, two_input_base+two_partition[start], infinity)
                add(left, two_input_base+two_partition[start+1], infinity)
                add(left, two_input_base+two_partition[start+2], infinity)
                add(left, two_input_base+two_partition[start+3], infinity)
    build_seconds = time.monotonic()-build_started
    started = time.monotonic()
    flow = graph.max_flow(source, sink)
    flow_seconds = time.monotonic()-started
    embedding = _sm_validate_matching(graph, h_masks, partitions, width, slots, starts, flow)
    reached = graph.reachable(source)
    if reached[sink]:
        raise ArithmeticError('Residual sink remains reachable')
    masks, reached_inputs = [], []
    for t in range(updates):
        mask = input_mask = 0
        base = t*width
        for coordinate in range(width):
            if reached[base+coordinate]:
                mask |= 1 << coordinate
            if reached[slots+base+coordinate]:
                input_mask |= 1 << coordinate
        masks.append(mask)
        reached_inputs.append(input_mask)
    hall = short_matching_hall_identity(masks, h_masks, partitions, width)
    cut = graph.cut_capacity(reached)
    if cut != flow or hall['cut_capacity'] != flow or not 0 <= flow <= slots:
        raise ArithmeticError('Independent Hall identity disagrees with maximum matching')
    for t, neighbors in enumerate(hall['neighbor_masks']):
        if neighbors & ~reached_inputs[t] or (reached_inputs[t] & ~neighbors & ~h_masks[t]):
            raise ArithmeticError('Residual right set disagrees with the canonical neighborhood')
    deficit = slots-flow
    if deficit != hall['integer_deficit_witness']:
        raise ArithmeticError('Hall witness does not attain the matching deficit')
    return dict(width=width, updates=updates, matching_flow=flow, hall_cut_capacity=cut,
                integer_deficit=deficit, physical_tau=deficit/2,
                physical_tau_upper=deficit/2,  # Retained compatibility alias; equality now proved.
                hall_source_size=hall['hall_source_size'],
                hall_noneligible_neighbors=hall['hall_noneligible_neighbors'],
                hall_neighbor_size=hall['hall_neighbor_size'], hall_source_masks=tuple(masks),
                graph_nodes=graph.nodes, original_edges=len(graph.to)//2,
                dinic_phases=graph.phases, unit_augmentations=graph.augmentations,
                build_seconds=build_seconds, flow_seconds=flow_seconds,
                independent_hall_identity_verified=True, maximum_match_delay=2,
                equals_full_flow_by_path_normal_form=True,
                neighbors_outside_window=False, cut_is_not_an_actual_support_obstruction=True,
                full_code_certificate=False, **embedding)


def short_matching_probe(updates=32, seed=0, compare_full=False):
    """Same sampled window as the full-flow helper; no independent resampling."""
    updates, seed = index(updates), index(seed)
    if updates not in (32, 64):
        raise ValueError('Production diagnostics permit only32 or64 updates')
    H, partitions = _sm_dependency('fixed_row_flow_sample')(SM_WIDTH, updates, SM_REGION, 4096, seed)
    quotas = tuple(sum(H[t].bit_count() for t in range(begin, begin+SM_REGION))
                   for begin in range(0, updates, SM_REGION))
    if any(quota != 4096 for quota in quotas):
        raise ArithmeticError('Regional H quota is not exact')
    result = short_matching_solve(H, partitions, SM_WIDTH)
    greedy = _sm_dependency('fixed_row_flow_greedy')(H, partitions, SM_WIDTH)
    if result['integer_deficit'] > greedy['greedy_integer_loss']:
        raise ArithmeticError('Optimal matching is worse than its feasible greedy subfamily')
    result.update(greedy)
    if compare_full:
        full = _sm_dependency('fixed_row_flow_solve')(H, partitions, SM_WIDTH)
        if result['matching_flow'] != full['max_flow']:
            raise ArithmeticError('Short matching disagrees with the full-flow equality theorem')
        result.update(full_flow_value=full['max_flow'], full_flow_integer_deficit=full['integer_deficit'],
                      full_flow_physical_tau=full['physical_tau'], same_instance_full_flow_compared=True)
    result.update(seed=seed, eligible_per_region=4096, sampled_regional_quotas=quotas,
                  global_physical_threshold=SM_GLOBAL_THRESHOLD,
                  linear_window_threshold=SM_GLOBAL_THRESHOLD*updates/SM_GLOBAL_UPDATES,
                  exceeds_linear_window_share=(result['integer_deficit']*SM_GLOBAL_UPDATES
                                               > 2*SM_GLOBAL_THRESHOLD*updates),
                  upper_exceeds_linear_window_share=(result['integer_deficit']*SM_GLOBAL_UPDATES
                                                    > 2*SM_GLOBAL_THRESHOLD*updates),
                  comparison_is_not_a_global_certificate=True,
                  hall_source_masks=tuple(hex(mask) for mask in result['hall_source_masks']))
    return result


def short_matching_checks():
    """Exhaust Hall sets and full-flow support cuts independently on tiny layers."""
    quartet = ((0, 1, 2, 3),)
    balanced = [sum(1 << coordinate for coordinate in chosen) for chosen in combinations(range(8), 4)]
    cases = [(2, mask) for mask in [0, 255]+balanced]
    cases += [(3, mask) for mask in (0, 4095, 15, 240, 3840, 1365, 2730, 255, 4080)]
    histories = minima = embedding_checks = terminal_checks = 0
    for updates, eligible in cases:
        H = tuple((eligible >> (4*t)) & 15 for t in range(updates))
        partitions = (quartet,)*updates
        slots, hall_min, full_min = 4*updates, 4*updates, 4*updates
        for packed in range(1 << slots):
            masks = tuple((packed >> (4*t)) & 15 for t in range(updates))
            volume = sum(mask.bit_count() for mask in masks)
            neighbors_outside = actual_outside = deaths = 0
            for t, mask in enumerate(masks):
                previous = masks[t-1] if t else 0
                before_previous = masks[t-2] if t >= 2 else 0
                # Width4 has one quartet, so these direct identities use no graph lookup.
                hall_neighbors = mask | (15 if previous or before_previous else 0)
                actual_footprint = mask | (15 if previous else 0)
                neighbors_outside += (hall_neighbors & ~H[t]).bit_count()
                actual_outside += (actual_footprint & ~H[t]).bit_count()
                deaths += (previous & ~mask).bit_count()
            hall_cut = slots-volume+2*neighbors_outside
            full_cut = slots-volume+deaths+2*actual_outside
            hall_min, full_min = min(hall_min, hall_cut), min(full_min, full_cut)
            assert short_matching_hall_identity(masks, H, partitions, 4)['cut_capacity'] == hall_cut
            histories += 1
        matching = short_matching_solve(H, partitions, 4)
        full = _sm_dependency('fixed_row_flow_solve')(H, partitions, 4)
        greedy = _sm_dependency('fixed_row_flow_greedy')(H, partitions, 4)
        assert matching['matching_flow'] == hall_min
        assert full['max_flow'] == full_min
        assert greedy['greedy_feasible_flow'] <= hall_min == full_min
        assert matching['full_flow_embedding_verified']
        minima += 1
        embedding_checks += 1
    for H, expected in (((0,), 0), ((15,), 4), ((0, 15), 4), ((15, 0), 0)):
        matching = short_matching_solve(H, (quartet,)*len(H), 4)
        assert matching['integer_deficit'] == expected
        terminal_checks += 1
    # Nontrivial varying partitions test the two-ahead path embedding and inequalities.
    for seed in range(8):
        H, partitions = _sm_dependency('fixed_row_flow_sample')(8, 4, 2, 8, seed)
        matching = short_matching_solve(H, partitions, 8)
        full = _sm_dependency('fixed_row_flow_solve')(H, partitions, 8)
        greedy = _sm_dependency('fixed_row_flow_greedy')(H, partitions, 8)
        assert greedy['greedy_feasible_flow'] <= matching['matching_flow'] == full['max_flow']
        assert matching['full_flow_embedding_verified']
        altered = (((0, 2, 4, 6), (1, 3, 5, 7)),)+partitions[1:]
        assert short_matching_solve(H, altered, 8)['matching_flow'] == matching['matching_flow']
        embedding_checks += 1
        terminal_checks += 1
    print(dict(exact_hall_source_histories=histories, exact_exhaustive_minima=minima,
               matching_full_flow_embedding_checks=embedding_checks,
               exact_terminal_checks=terminal_checks, capacity_two=True,
               no_numerical_probability_claim=True, full_code_certificate=False), flush=True)


def short_matching_main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--checks', action='store_true')
    mode.add_argument('--probe', action='store_true')
    parser.add_argument('--updates', type=int, choices=(32, 64), default=32)
    parser.add_argument('--seed', type=int, default=0)
    parser.add_argument('--compare-full', action='store_true')
    args = parser.parse_args()
    import resource
    if shutil.disk_usage('/tmp').free < 20*(1 << 30):
        raise RuntimeError('20GiB /tmp free-space guard')
    resource.setrlimit(resource.RLIMIT_AS, (2*(1 << 30), 2*(1 << 30)))
    timer = threading.Timer(175, lambda: os._exit(124))
    timer.daemon = True
    timer.start()
    started = time.monotonic()
    try:
        if args.checks:
            short_matching_checks()
        else:
            result = short_matching_probe(args.updates, args.seed, args.compare_full)
            print(dict(seconds=round(time.monotonic()-started, 3), **result), flush=True)
    finally:
        timer.cancel()


if __name__ == '__main__':
    short_matching_main()
