r"""One bounded fixed-H flow diagnostic; no field amplitudes or encoder timing.

For a window of T updates and n coordinates, use state and input vertices
at each(t,i). Add source->state capacity1; state->same-coordinate next
state capacity1; state->own input and each input of its NEXT quartet
capacity INF>n*T; non-H input->sink capacity2. No edge crosses either
window boundary. In particular there is no final-state death or next hit.

For the source-side state masks S, canonical V_t=S_t union Gamma_t(Sprev)
gives cut capacity n*T-L+D+2|V\H|. Thus the integer flow deficit is
max_S(L-D-2|V\H|), and the physical tau is half this deficit. The returned
cut is checked independently against its canonical support/footprint.
The same instance also has a backward greedy feasible flow. Its loss is
an upper bound on the exact deficit; it never uses the temporal arcs.

Each16-update region samples exactly4096 H positions for n512. A local
random.Random(seed) samples these subsets and fresh quartet partitions;
the first partition is sampled but unused. The model draws no outside
neighbors. This is one structural window diagnostic, not a probability
certificate, global obstruction witness, or encoder benchmark.

Production CLI permits exactly one32/64-update width512 window. Exact
checks exhaust all256 state histories of width4,T2 for every balanced
H mask, plus all-H/all-non-H controls. No sampling or job runs on import.
"""

import argparse
from array import array
from collections import deque
from itertools import combinations
from operator import index
import os
import random
import shutil
import threading
import time


FR_WIDTH, FR_REGION, FR_GLOBAL_UPDATES = 512, 16, 4096
FR_GLOBAL_THRESHOLD = 52430


class _FRDinic:
    """Compact residual arrays; every original source edge has unit capacity."""

    __slots__ = ('head', 'to', 'following', 'capacity', 'nodes', 'phases', 'augmentations')

    def __init__(self, nodes):
        self.nodes = nodes
        self.head = [-1]*nodes
        self.to, self.following, self.capacity = array('i'), array('i'), array('i')
        self.phases = self.augmentations = 0

    def add_edge(self, source, target, capacity):
        edge = len(self.to)
        self.to.append(target)
        self.following.append(self.head[source])
        self.capacity.append(capacity)
        self.head[source] = edge
        self.to.append(source)
        self.following.append(self.head[target])
        self.capacity.append(0)
        self.head[target] = edge+1

    def _levels(self, source, sink):
        head, to, following, capacity = self.head, self.to, self.following, self.capacity
        level = [-1]*self.nodes
        level[source] = 0
        queue = deque((source,))
        while queue:
            vertex = queue.popleft()
            next_level = level[vertex]+1
            edge = head[vertex]
            while edge != -1:
                target = to[edge]
                if capacity[edge] > 0 and level[target] < 0:
                    level[target] = next_level
                    queue.append(target)
                edge = following[edge]
        return level if level[sink] >= 0 else None

    def max_flow(self, source, sink):
        """Iterative blocking paths; integer unit source arcs make every push1."""
        to, following, capacity = self.to, self.following, self.capacity
        edge = self.head[source]
        while edge != -1:
            if capacity[edge] != 1:
                raise ValueError('This bounded Dinic implementation requires unit source arcs')
            edge = following[edge]
        flow = 0
        while True:
            level = self._levels(source, sink)
            if level is None:
                return flow
            self.phases += 1
            current = self.head.copy()
            vertices, path = [source], []
            while vertices:
                vertex = vertices[-1]
                if vertex == sink:
                    # A positive integer path starts with an unsaturated unit arc.
                    for chosen in path:
                        capacity[chosen] -= 1
                        capacity[chosen ^ 1] += 1
                    flow += 1
                    self.augmentations += 1
                    del vertices[1:]
                    path.clear()
                    continue
                edge = current[vertex]
                while edge != -1 and (capacity[edge] <= 0 or level[to[edge]] != level[vertex]+1):
                    edge = following[edge]
                current[vertex] = edge
                if edge == -1:
                    # No new forward residual capacity appears in this level graph.
                    level[vertex] = -1
                    vertices.pop()
                    if path:
                        failed = path.pop()
                        current[vertices[-1]] = following[failed]
                else:
                    path.append(edge)
                    vertices.append(to[edge])

    def reachable(self, source):
        reached = bytearray(self.nodes)
        reached[source] = 1
        queue = deque((source,))
        head, to, following, capacity = self.head, self.to, self.following, self.capacity
        while queue:
            vertex = queue.popleft()
            edge = head[vertex]
            while edge != -1:
                target = to[edge]
                if capacity[edge] > 0 and not reached[target]:
                    reached[target] = 1
                    queue.append(target)
                edge = following[edge]
        return reached

    def cut_capacity(self, reached):
        """Recover each original forward capacity from its residual pair."""
        total = 0
        for vertex in range(self.nodes):
            if not reached[vertex]:
                continue
            edge = self.head[vertex]
            while edge != -1:
                if not edge & 1 and not reached[self.to[edge]]:
                    total += self.capacity[edge]+self.capacity[edge ^ 1]
                edge = self.following[edge]
        return total


def fixed_row_flow_sample(width=512, updates=32, region_updates=16,
                          eligible_per_region=None, seed=0):
    """Sample one reproducible instance; H labels first, then all partitions."""
    width, updates, region_updates, seed = map(index, (width, updates, region_updates, seed))
    if width < 4 or width % 4 or updates < 1 or region_updates < 1 or updates % region_updates:
        raise ValueError('Require quartet width and a positive integral number of regions')
    regional_slots = width*region_updates
    eligible = regional_slots//2 if eligible_per_region is None else index(eligible_per_region)
    if not 0 <= eligible <= regional_slots:
        raise ValueError('Regional H quota is outside the slot range')
    rng = random.Random(seed)
    h_masks = [0]*updates
    for begin in range(0, updates, region_updates):
        for slot in rng.sample(range(regional_slots), eligible):
            offset, coordinate = divmod(slot, width)
            h_masks[begin+offset] |= 1 << coordinate
    partitions = []
    for _ in range(updates):
        permutation = list(range(width))
        rng.shuffle(permutation)
        partitions.append(tuple(tuple(permutation[start:start+4]) for start in range(0, width, 4)))
    return tuple(h_masks), tuple(partitions)


def _fr_validate_instance(h_masks, partitions, width):
    width = index(width)
    if width < 4 or width % 4 or not h_masks or len(h_masks) != len(partitions):
        raise ValueError('Require nonempty matching H/partition sequences and quartet width')
    full = (1 << width)-1
    for mask, partition in zip(h_masks, partitions):
        if not isinstance(mask, int) or mask < 0 or mask & ~full:
            raise ValueError('H mask is outside the coordinate set')
        if len(partition) != width//4 or any(len(quartet) != 4 for quartet in partition):
            raise ValueError('Each partition must consist of quartets')
        if sorted(coordinate for quartet in partition for coordinate in quartet) != list(range(width)):
            raise ValueError('Quartet coordinates must partition the whole layer')
    return width, len(h_masks)


def fixed_row_flow_cut_identity(state_masks, h_masks, partitions, width):
    """Recompute canonical footprints directly from masks, independently of graph arcs."""
    width, updates = _fr_validate_instance(h_masks, partitions, width)
    if len(state_masks) != updates:
        raise ValueError('State history length mismatch')
    full, previous = (1 << width)-1, 0
    volume = deaths = outside = footprint = 0
    footprint_masks = []
    for t, state in enumerate(state_masks):
        if not isinstance(state, int) or state < 0 or state & ~full:
            raise ValueError('State mask is outside the coordinate set')
        live = 0
        if t:
            for quartet in partitions[t]:
                mask = sum(1 << coordinate for coordinate in quartet)
                if previous & mask:
                    live |= mask
        covered = live | state
        volume += state.bit_count()
        deaths += (previous & ~state).bit_count()
        outside += (covered & ~h_masks[t]).bit_count()
        footprint += covered.bit_count()
        footprint_masks.append(covered)
        previous = state
    # previous starts empty; there is no iteration after the final state.
    capacity = width*updates-volume+deaths+2*outside
    return dict(volume=volume, deaths=deaths, covered_outside=outside, footprint=footprint,
                cut_capacity=capacity, integer_score=volume-deaths-2*outside,
                footprint_masks=tuple(footprint_masks))


def fixed_row_flow_greedy(h_masks, partitions, width=512):
    """Feasible backward flow using own/future input arcs, never temporal arcs."""
    width, updates = _fr_validate_instance(h_masks, partitions, width)
    last = h_masks[-1]
    remaining = [0 if (last >> i) & 1 else 1 for i in range(width)]
    losses = [0]*updates
    losses[-1] = last.bit_count()
    for t in range(updates-2, -1, -1):
        current = [0]*width
        H = h_masks[t]
        for quartet in partitions[t+1]:
            capacity = sum(remaining[i] for i in quartet)
            eligible = sum((H >> i) & 1 for i in quartet)
            losses[t] += max(0, eligible-capacity)
            noneligible_future = min(4-eligible, max(0, capacity-eligible))
            for coordinate in quartet:
                if (H >> coordinate) & 1:
                    current[coordinate] = 0
                elif noneligible_future:
                    current[coordinate] = 2
                    noneligible_future -= 1
                else:
                    current[coordinate] = 1
        remaining = current
    loss = sum(losses)
    return dict(greedy_integer_loss=loss, greedy_physical_upper=loss/2,
                greedy_feasible_flow=width*updates-loss,
                greedy_losses_by_update=tuple(losses), greedy_used_temporal_arcs=False)


def fixed_row_flow_solve(h_masks, partitions, width=512):
    """Solve one supplied instance and return its source-side state masks."""
    width, updates = _fr_validate_instance(h_masks, partitions, width)
    slots = width*updates
    if slots > 32768:
        raise ValueError('This diagnostic is bounded to at most32768 state vertices')
    source, sink, infinity = 2*slots, 2*slots+1, slots+1
    graph = _FRDinic(2*slots+2)
    quartet_start, flattened = [], []
    for partition in partitions:
        flat = tuple(coordinate for quartet in partition for coordinate in quartet)
        starts = [0]*width
        for start in range(0, width, 4):
            starts[flat[start]] = start
            starts[flat[start+1]] = start
            starts[flat[start+2]] = start
            starts[flat[start+3]] = start
        flattened.append(flat)
        quartet_start.append(starts)
    build_started = time.monotonic()
    add_edge = graph.add_edge
    for t in range(updates):
        current_base = t*width
        has_next = t+1 < updates
        if has_next:
            next_base = current_base+width
            next_input_base = slots+next_base
            next_partition, starts = flattened[t+1], quartet_start[t+1]
        H = h_masks[t]
        for coordinate in range(width):
            state, own_input = current_base+coordinate, slots+current_base+coordinate
            add_edge(source, state, 1)
            add_edge(state, own_input, infinity)
            if not (H >> coordinate) & 1:
                add_edge(own_input, sink, 2)
            if has_next:
                add_edge(state, next_base+coordinate, 1)
                start = starts[coordinate]
                add_edge(state, next_input_base+next_partition[start], infinity)
                add_edge(state, next_input_base+next_partition[start+1], infinity)
                add_edge(state, next_input_base+next_partition[start+2], infinity)
                add_edge(state, next_input_base+next_partition[start+3], infinity)
    build_seconds = time.monotonic()-build_started
    flow_started = time.monotonic()
    flow = graph.max_flow(source, sink)
    flow_seconds = time.monotonic()-flow_started
    reached = graph.reachable(source)
    if reached[sink]:
        raise ArithmeticError('Residual sink remains reachable')
    state_masks, input_masks = [], []
    for t in range(updates):
        state_mask = input_mask = 0
        base = t*width
        for coordinate in range(width):
            if reached[base+coordinate]:
                state_mask |= 1 << coordinate
            if reached[slots+base+coordinate]:
                input_mask |= 1 << coordinate
        state_masks.append(state_mask)
        input_masks.append(input_mask)
    identity = fixed_row_flow_cut_identity(state_masks, h_masks, partitions, width)
    graph_cut = graph.cut_capacity(reached)
    if flow != graph_cut or flow != identity['cut_capacity'] or not 0 <= flow <= slots:
        raise ArithmeticError('Independent cut identity disagrees with maximum flow')
    for t, footprint in enumerate(identity['footprint_masks']):
        if footprint & ~input_masks[t] or (input_masks[t] & ~footprint & ~h_masks[t]):
            raise ArithmeticError('Source-side inputs disagree with the canonical footprint')
    deficit = slots-flow
    if deficit != identity['integer_score']:
        raise ArithmeticError('Integer deficit disagrees with the support score')
    return dict(width=width, updates=updates, max_flow=flow, cut_capacity=graph_cut,
                integer_deficit=deficit, physical_tau=deficit/2,
                volume=identity['volume'], deaths=identity['deaths'],
                footprint=identity['footprint'], covered_outside=identity['covered_outside'],
                state_masks=tuple(state_masks), graph_nodes=graph.nodes,
                original_edges=len(graph.to)//2, dinic_phases=graph.phases,
                unit_augmentations=graph.augmentations,
                build_seconds=build_seconds, flow_seconds=flow_seconds,
                independent_cut_identity_verified=True, terminal_reset=False,
                neighbors_outside_window=False, global_obstruction_claim=False,
                full_code_certificate=False)


def fixed_row_flow_probe(updates=32, seed=0):
    """One width512 window, with exact4096-slot H quotas in each region."""
    updates, seed = index(updates), index(seed)
    if updates not in (32, 64):
        raise ValueError('The bounded production probe permits only32 or64 updates')
    h_masks, partitions = fixed_row_flow_sample(FR_WIDTH, updates, FR_REGION, 4096, seed)
    quotas = tuple(sum(h_masks[t].bit_count() for t in range(begin, begin+FR_REGION))
                   for begin in range(0, updates, FR_REGION))
    if any(quota != 4096 for quota in quotas):
        raise ArithmeticError('Sampled H quota is not exact')
    result = fixed_row_flow_solve(h_masks, partitions, FR_WIDTH)
    greedy = fixed_row_flow_greedy(h_masks, partitions, FR_WIDTH)
    if result['integer_deficit'] > greedy['greedy_integer_loss']:
        raise ArithmeticError('Greedy feasible loss is below the exact deficit')
    result.update(greedy)
    result.update(seed=seed, eligible_per_region=4096, sampled_regional_quotas=quotas,
                  global_physical_threshold=FR_GLOBAL_THRESHOLD,
                  linear_window_threshold=FR_GLOBAL_THRESHOLD*updates/FR_GLOBAL_UPDATES,
                  exceeds_linear_window_share=(result['integer_deficit']*FR_GLOBAL_UPDATES
                                               > 2*FR_GLOBAL_THRESHOLD*updates),
                  comparison_is_not_a_global_certificate=True,
                  state_masks=tuple(hex(mask) for mask in result['state_masks']))
    return result


def fixed_row_flow_checks():
    """Exhaust tiny support sets independently of the residual algorithm."""
    partition = ((0, 1, 2, 3),)
    partitions = (partition, partition)
    balanced = [sum(1 << coordinate for coordinate in chosen) for chosen in combinations(range(8), 4)]
    h_values = [0, 255]+balanced
    histories = exact_minima = terminal_checks = quota_checks = greedy_checks = 0
    for eligible in h_values:
        h_masks = (eligible & 15, eligible >> 4)
        expected = 8
        for bits in range(256):
            states = (bits & 15, bits >> 4)
            # Direct two-layer identity, with no graph or shared footprint helper.
            first, last = states
            first_v, last_v = first, last | (15 if first else 0)
            volume = first.bit_count()+last.bit_count()
            death = (first & ~last).bit_count()
            outside = (first_v & ~h_masks[0]).bit_count()+(last_v & ~h_masks[1]).bit_count()
            capacity = 8-volume+death+2*outside
            expected = min(expected, capacity)
            independent = fixed_row_flow_cut_identity(states, h_masks, partitions, 4)
            assert independent['cut_capacity'] == capacity
            histories += 1
        solved = fixed_row_flow_solve(h_masks, partitions, 4)
        assert solved['max_flow'] == expected
        assert solved['integer_deficit'] == 8-expected
        greedy = fixed_row_flow_greedy(h_masks, partitions, 4)
        assert solved['integer_deficit'] <= greedy['greedy_integer_loss']
        greedy_checks += 1
        exact_minima += 1
    for h_masks, expected_deficit in (((0,), 0), ((15,), 4), ((0, 15), 4), ((15, 0), 0), ((15, 15), 8)):
        solved = fixed_row_flow_solve(h_masks, (partition,)*len(h_masks), 4)
        assert solved['integer_deficit'] == expected_deficit
        terminal_checks += 1
    # Changing the unused first partition cannot change the graph or cut.
    H, sampled = fixed_row_flow_sample(8, 4, 2, 8, 17)
    alternative = (((0, 2, 4, 6), (1, 3, 5, 7)),)+sampled[1:]
    left, right = fixed_row_flow_solve(H, sampled, 8), fixed_row_flow_solve(H, alternative, 8)
    assert left['max_flow'] == right['max_flow']
    assert left['integer_deficit'] <= fixed_row_flow_greedy(H, sampled, 8)['greedy_integer_loss']
    greedy_checks += 1
    terminal_checks += 1
    for width, updates, region, seed in ((4, 2, 2, 0), (8, 4, 2, 1), (8, 4, 2, 19)):
        H, partitions = fixed_row_flow_sample(width, updates, region, seed=seed)
        H2, partitions2 = fixed_row_flow_sample(width, updates, region, seed=seed)
        assert (H, partitions) == (H2, partitions2)
        for begin in range(0, updates, region):
            assert sum(H[t].bit_count() for t in range(begin, begin+region)) == width*region//2
            quota_checks += 1
        _fr_validate_instance(H, partitions, width)
        solved = fixed_row_flow_solve(H, partitions, width)
        assert solved['integer_deficit'] <= fixed_row_flow_greedy(H, partitions, width)['greedy_integer_loss']
        greedy_checks += 1
    assert FR_GLOBAL_THRESHOLD*32/FR_GLOBAL_UPDATES == 409.609375
    print(dict(exact_state_histories=histories, exact_minimum_cuts=exact_minima,
               exact_terminal_checks=terminal_checks, exact_quota_checks=quota_checks,
               exact_greedy_upper_checks=greedy_checks,
               independent_cut_identity=True, no_numerical_probability_claim=True,
               full_code_certificate=False), flush=True)


def fixed_row_flow_main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--checks', action='store_true')
    mode.add_argument('--probe', action='store_true')
    parser.add_argument('--updates', type=int, choices=(32, 64), default=32)
    parser.add_argument('--seed', type=int, default=0)
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
            fixed_row_flow_checks()
        else:
            result = fixed_row_flow_probe(args.updates, args.seed)
            print(dict(seconds=round(time.monotonic()-started, 3), **result), flush=True)
    finally:
        timer.cancel()


if __name__ == '__main__':
    fixed_row_flow_main()
