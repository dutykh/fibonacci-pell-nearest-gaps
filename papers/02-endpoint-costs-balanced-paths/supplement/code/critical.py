# Authors: Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science
# and Technology, Abu Dhabi, UAE) and Prof. Laurent Vuillon (Univ. Savoie Mont
# Blanc, CNRS, LAMA, Chambery, France).
# Copyright (C) 2026 Dr. Denys Dutykh and Prof. Laurent Vuillon.
# Distributed under the GNU Lesser General Public License, version 2.1;
# see the LICENSE file of this package.
"""Critical zero-slack components and an exact longest-path bound."""
from collections import deque
from fractions import Fraction
from time import monotonic
import elementary as h

B = 89


def check(endpoint_record, graph):
    started = monotonic()
    states, index, edges, seeds, dist, accepting = graph
    n = len(states)
    order = sorted(range(n), key=states.__getitem__)
    graph_hash = endpoint_record["graph_sha256"]
    potential_hash = endpoint_record["potential_sha256"]
    assert graph_hash == "58cace292af6432c957548afa9da02e76d1cfa1e446d53212ef2067504b48d02"
    assert potential_hash == "e341181736aafbcb3faef3a922d4b7672b7ac1082a0cfbf433ed27717462170d"
    zero, reverse = [[] for _ in states], [[] for _ in states]
    for i, row in enumerate(edges):
        for letter, (j, w, _) in enumerate(row, 1):
            slack = dist[i] + w - dist[j]
            assert slack >= 0
            if slack == 0:
                zero[i].append((j, letter))
                reverse[j].append(i)

    # Genuine postorder: each vertex is appended only when its own edge
    # iterator is exhausted. Merely reversing discovery order is incorrect.
    seen, finished = set(), []
    for start in order:
        if start in seen:
            continue
        seen.add(start)
        stack = [(start, iter(zero[start]))]
        while stack:
            u, iterator = stack[-1]
            successor = next(iterator, None)
            if successor is None:
                stack.pop()
                finished.append(u)
            elif successor[0] not in seen:
                v = successor[0]
                seen.add(v)
                stack.append((v, iter(zero[v])))
    assert len(finished) == n
    component, groups = [-1] * n, []
    for start in reversed(finished):
        if component[start] >= 0:
            continue
        number = len(groups)
        component[start] = number
        group, stack = [], [start]
        while stack:
            u = stack.pop()
            group.append(u)
            for v in reverse[u]:
                if component[v] < 0:
                    component[v] = number
                    stack.append(v)
        groups.append(sorted(group, key=states.__getitem__))

    dag = [set() for _ in groups]
    indegree = [0] * len(groups)
    for u, row in enumerate(zero):
        for v, _ in row:
            a, b = component[u], component[v]
            if a != b and b not in dag[a]:
                dag[a].add(b)
                indegree[b] += 1
    queue = deque(i for i, val in enumerate(indegree) if val == 0)
    topo = []
    while queue:
        a = queue.popleft()
        topo.append(a)
        for b in dag[a]:
            indegree[b] -= 1
            if indegree[b] == 0:
                queue.append(b)
    assert len(topo) == len(groups), "SCC condensation must be acyclic"

    critical, node_weight = [], [Fraction(0) for _ in groups]
    for number, group in enumerate(groups):
        if len(group) == 1 and not any(v == group[0] for v, _ in zero[group[0]]):
            continue
        root = group[0]
        potentials = {root: (0, 0)}
        queue = deque((root,))
        while queue:
            u = queue.popleft()
            length, twos = potentials[u]
            for v, letter in zero[u]:
                if component[v] == number and v not in potentials:
                    potentials[v] = (length + 1, twos + (letter == 2))
                    queue.append(v)
        assert set(potentials) == set(group)
        constraints = []
        for u in group:
            for v, letter in zero[u]:
                if component[v] == number:
                    du = potentials[u][0] + 1 - potentials[v][0]
                    dm = potentials[u][1] + (letter == 2) - potentials[v][1]
                    constraints.append((du, dm))
        theta = next(Fraction(dm, du) for du, dm in constraints if du)
        assert 0 <= theta <= 1
        assert all(Fraction(dm) == theta * du for du, dm in constraints)
        f = {u: Fraction(m) - theta * length for u, (length, m) in potentials.items()}
        assert all(Fraction(letter == 2) - theta == f[v] - f[u]
                   for u in group for v, letter in zero[u] if component[v] == number)
        oscillation = max(f.values()) - min(f.values())
        node_weight[number] = oscillation + 1
        critical.append({"size": len(group), "theta": str(theta), "oscillation": str(oscillation),
                         "states": [[list(states[u]), str(f[u])] for u in group],
                         "edges": [[list(states[u]), letter, list(states[v])]
                                   for u in group for v, letter in zero[u]
                                   if component[v] == number]})
    assert sorted(record["size"] for record in critical) == [1, 1, 1, 1, 3, 3, 6]
    assert {record["theta"] for record in critical} == {"0", "1/5", "1/3", "1"}
    longest = node_weight.copy()
    predecessor = [None] * len(groups)
    for a in topo:
        for b in dag[a]:
            candidate = longest[a] + 1 + node_weight[b]
            if candidate > longest[b]:
                longest[b], predecessor[b] = candidate, a
    end = max(range(len(groups)), key=longest.__getitem__)
    L = longest[end]
    witness_components = []
    while end is not None:
        witness_components.append(end)
        end = predecessor[end]
    witness_components.reverse()
    assert sum((node_weight[i] for i in witness_components), Fraction()) + len(witness_components) - 1 == L
    assert all(longest[b] >= longest[a] + 1 + node_weight[b]
               for a, row in enumerate(dag) for b in row)
    h.check_clock(started, "critical check")
    result = {"status": "PASS", "B": B, "states": n, "graph_sha256": graph_hash,
              "potential_sha256": potential_hash, "zero_edges": sum(map(len, zero)),
              "components": len(groups), "cyclic_components": critical,
              "L": str(L), "K_graph": str(L + 1), "maximum_seed_prefix": 2,
              "K_full_word": str(L + 3),
              "longest_path_component_sizes": [len(groups[i]) for i in witness_components],
              "longest_path_component_weights": [str(node_weight[i]) for i in witness_components]}
    return result
