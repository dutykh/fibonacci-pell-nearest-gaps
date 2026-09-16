# Authors: Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science
# and Technology, Abu Dhabi, UAE) and Prof. Laurent Vuillon (Univ. Savoie Mont
# Blanc, CNRS, LAMA, Chambery, France).
# Copyright (C) 2026 Dr. Denys Dutykh and Prof. Laurent Vuillon.
# Distributed under the GNU Lesser General Public License, version 2.1;
# see the LICENSE file of this package.
"""Complete endpoint potential, primitive terminals, and exact equality language."""
from collections import deque
from hashlib import sha256
import json
from math import gcd
from time import monotonic
import elementary as h

B = 89
CAP = 24000

def canon(s):
    return s if h.det(s) > 0 else (s[2], s[3], s[0], s[1])


def fingerprint(obj):
    return sha256(json.dumps(obj, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def reach(initial, adjacency):
    found = set(initial)
    queue = deque(initial)
    while queue:
        for target in adjacency[queue.popleft()]:
            if target not in found:
                found.add(target)
                queue.append(target)
    return found


def block(word):
    result = h.ONE
    for letter in word:
        result = h.mul(result, h.DOUBLE[letter])
    return result


def language_audit(start, live_edges, terminal):
    tokens = ((1, 2, 1, 1, 1), (2, 1, 1, 1, 1))
    prefixes = {word[:i] for word in tokens for i in range(len(word))}

    def advance(state, letter):
        if state == "initial":
            return () if letter == 1 else None
        if state is None:
            return None
        extended = state + (letter,)
        if extended in tokens:
            return ()
        return extended if extended in prefixes else None

    initial = (start, "initial")
    found, queue = {initial}, deque((initial,))
    while queue:
        matrix_state, regex_state = queue.popleft()
        assert (matrix_state == terminal) == (regex_state == ())
        for letter in (1, 2):
            target = (live_edges.get((matrix_state, letter)), advance(regex_state, letter))
            if target not in found:
                assert len(found) < 32
                found.add(target)
                queue.append(target)
    return len(found)


def check():
    started = monotonic()
    reflections = h.matrices(8, 5, B)
    original_seeds, direct, visits, depth = h.prefix_partition(reflections, B, started)
    assert len(original_seeds) == 11 and not direct
    seeds = [h.Seed(s.branch, s.optional, s.prefix, canon(s.state), s.weight)
             for s in original_seeds]
    states, lookup, edges = [], {}, []

    def add(s):
        if s not in lookup:
            if len(states) >= CAP:
                raise RuntimeError("INCONCLUSIVE: quotient-state cap exceeded")
            lookup[s] = len(states)
            states.append(s)
        return lookup[s]

    for seed in seeds:
        add(seed.state)
    cursor = 0
    while cursor < len(states):
        if cursor % 128 == 0:
            h.check_clock(started, "independent axis89 closure")
        state = states[cursor]
        h.audit_reduced_state(state, B)
        assert h.content(state) == 1 and h.det(state) == B * B
        row = []
        for letter in (1, 2):
            edge = h.step(state, letter)
            # Independently verify equivariance of the row-exchange quotient.
            swapped = (state[2], state[3], state[0], state[1])
            other = h.step(swapped, letter)
            assert canon(other.target) == canon(edge.target)
            assert other.emitted == edge.emitted
            row.append((add(canon(edge.target)), edge.weight, edge.emitted))
        edges.append(row)
        cursor += 1
    assert len(states) == 11420 and len(edges) == 11420
    n = len(states)
    infinity = 10**100
    distance = [infinity] * n
    for seed in seeds:
        i = lookup[seed.state]
        distance[i] = min(distance[i], seed.weight)
    order = sorted(range(n), key=states.__getitem__)
    sweeps = 0
    while True:
        h.check_clock(started, "the ordered edge sweeps")
        sweeps += 1
        if sweeps > n:
            raise RuntimeError("INCONCLUSIVE: relaxation did not stabilize")
        changed = False
        for i in order:
            if distance[i] == infinity:
                continue
            for j, weight, _emitted in edges[i]:
                candidate = distance[i] + weight
                if candidate < distance[j]:
                    distance[j] = candidate
                    changed = True
        if not changed:
            break
    assert all(value < infinity for value in distance)
    assert all(seed.weight >= distance[lookup[seed.state]] for seed in seeds)
    assert all(distance[i] + weight >= distance[j]
               for i, row in enumerate(edges) for j, weight, _ in row)
    accepting = {i: h.terminal_weight(s) for i, s in enumerate(states)
                 if gcd(*h.terminal_vector(s)) == B}
    assert len(accepting) == 116
    assert min(distance[i] + final for i, final in accepting.items()) == 0
    assert all(distance[i] + final >= 0 for i, final in accepting.items())

    tight = {(i, letter): j for i, row in enumerate(edges)
             for letter, (j, weight, _) in enumerate(row, 1)
             if distance[i] + weight == distance[j]}
    adjacency, reverse = [[] for _ in states], [[] for _ in states]
    for (i, _), j in tight.items():
        adjacency[i].append(j)
        reverse[j].append(i)
    tight_seeds = [s for s in seeds if s.weight == distance[lookup[s.state]]]
    terminals = {i for i, final in accepting.items() if distance[i] + final == 0}
    live = reach([lookup[s.state] for s in tight_seeds], adjacency) & reach(terminals, reverse)
    live_seeds = [s for s in tight_seeds if lookup[s.state] in live]
    expected = ((119, 39, 41, 80), (115, 37, 47, 84),
                (89, 21, 0, 89), (89, 0, 68, 89),
                (131, 21, 47, 68), (121, 37, 41, 78))
    edge_indices = ((0, 1, 1), (1, 1, 2), (1, 2, 3), (2, 2, 4),
                    (3, 1, 4), (4, 1, 5), (5, 1, 0))
    assert live == {lookup[s] for s in expected}
    assert terminals & live == {lookup[expected[1]]}
    assert live_seeds == [h.Seed(0, 1, (), expected[0], -1)]
    live_edges = {(states[i], letter): states[j] for (i, letter), j in tight.items()
                  if i in live and j in live}
    assert live_edges == {(expected[i], letter): expected[j]
                          for i, letter, j in edge_indices}
    product_size = language_audit(expected[0], live_edges, expected[1])

    # Two token intertwiners and one terminal identity prove the target language
    # for arbitrarily many tokens. No bounded-word sampling is used.
    q0 = h.word_matrix((1, 1, 1))
    K = (47, 84, 68, -47)
    H = (89, 68, 0, -89)
    token0, token1 = (1, 2, 1, 1, 1), (2, 1, 1, 1, 1)
    assert h.mul(reflections[0], q0) == h.mul(q0, K)
    assert h.mul(K, h.J[2]) == h.mul(h.J[2], H)
    assert h.mul(K, block(token0)) == h.mul(block(token1), K)
    assert h.mul(K, block(token1)) == h.mul(block(token0), K)
    assert h.terminal_vector(q0) == (8, 5)
    assert h.mul(H, H) == (B * B, 0, 0, B * B)

    graph_record = [[list(states[i]), letter, list(states[j]), weight, emitted]
                    for i in order for letter, (j, weight, emitted) in enumerate(edges[i], 1)]
    potential_record = [[list(states[i]), distance[i]] for i in order]
    corridor = [[list(s), distance[lookup[s]],
                 h.terminal_weight(s) if lookup[s] in accepting else None]
                for s in expected]
    result = {"status": "PASS", "B": B, "states": n, "edges": 2 * n,
              "seeds": len(seeds), "prefix_visits": visits, "prefix_depth": depth,
              "direct": len(direct), "primitive_terminals": len(accepting),
              "minimum": 0, "edge_sweeps": sweeps,
              "graph_sha256": fingerprint(graph_record),
              "potential_sha256": fingerprint(potential_record),
              "live_states": corridor,
              "live_edges": [[list(expected[i]), letter, list(expected[j]),
                              edges[lookup[expected[i]]][letter - 1][1]]
                             for i, letter, j in edge_indices],
              "language_product_states": product_size,
              "language": "q = 111 chi((12111|21111)*) 2",
              "all_seeds": [[s.branch, s.optional, list(s.prefix), list(s.state), s.weight]
                            for s in seeds]}
    return result, (states, lookup, edges, seeds, distance, accepting)
