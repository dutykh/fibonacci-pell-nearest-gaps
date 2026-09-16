#!/usr/bin/env python3
# Authors: Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science
# and Technology, Abu Dhabi, UAE) and Prof. Laurent Vuillon (Univ. Savoie Mont
# Blanc, CNRS, LAMA, Chambery, France).
# Copyright (C) 2026 Dr. Denys Dutykh and Prof. Laurent Vuillon.
# Distributed under the GNU Lesser General Public License, version 2.1;
# see the LICENSE file of this package.
"""Independently reconstruct the fixed norm-89 endpoint and critical graph.

The only local mathematical import is code/unit_reducer.py. No function of the
first computation is used: this one performs every row subtraction singly, and
reaches the potential, the strong components and the condensation bound by a
work-list relaxation, an iterative strong-component search and sink-first
remaining budgets.
"""
from __future__ import annotations

from collections import deque
from dataclasses import replace
from fractions import Fraction
from hashlib import sha256
import json
from math import gcd
from pathlib import Path
import resource
import signal
import sys

WRITE_REFERENCE = sys.argv[1:] == ["--write-reference"]
if not __debug__ or (sys.argv[1:] and not WRITE_REFERENCE):
    raise SystemExit("Run without optimization; the only accepted option is --write-reference.")

HERE = Path(__file__).resolve().parent
STATE_CAP = 24000
RELAXATION_CAP = 1000000
B = 89


def limits():
    for kind, requested in ((resource.RLIMIT_AS, 128 * 1024 * 1024),
                            (resource.RLIMIT_CPU, 30)):
        soft, hard = resource.getrlimit(kind)
        ceiling = requested
        if soft != resource.RLIM_INFINITY:
            ceiling = min(ceiling, soft)
        if hard != resource.RLIM_INFINITY:
            ceiling = min(ceiling, hard)
        resource.setrlimit(kind, (ceiling, ceiling))

    def expired(_signum, _frame):
        raise RuntimeError("INCONCLUSIVE: the second computation reached its 30-second wall cap")

    signal.signal(signal.SIGALRM, expired)
    signal.alarm(30)


def read_data(name):
    path = HERE / "data" / name
    if not path.is_file():
        raise SystemExit(f"Missing data file: {path.name}")
    return json.loads(path.read_text(encoding="utf-8"))


def swap(s):
    return s[2], s[3], s[0], s[1]


def canon(s):
    # This lexicographic convention differs from the first computation's sign.
    return min(s, swap(s))


def positive(s):
    return s if s[0] * s[3] > s[1] * s[2] else swap(s)


def fingerprint(obj):
    return sha256(json.dumps(obj, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def reachable(initial, adjacency):
    found = set(initial)
    queue = deque(found)
    while queue:
        for v in adjacency[queue.popleft()]:
            if v not in found:
                found.add(v)
                queue.append(v)
    return found


def endpoint(g, reference, table):
    reflections = ((39, 80, 80, -39), (80, 39, 39, -80))
    seeds, direct, visits, depth = [], [], 0, 0
    for i, r in enumerate(reflections):
        assert g.determinant(r) == -B * B and gcd(*r) == 1
        ss, dd, vv, hh = g.split_prefixes(r, i)
        seeds.extend(replace(s, state=canon(s.state)) for s in ss)
        direct.extend(dd)
        visits += vv
        depth = max(depth, hh)
    assert visits <= 64
    assert (len(seeds), len(direct), visits, depth) == (11, 0, 16, 2)
    seed_rows = [[s.reflection, s.optional or None, list(s.prefix),
                  list(positive(s.state)), s.weight] for s in seeds]
    assert seed_rows == reference["all_seeds"]

    states, index, edges = [], {}, []

    def add(s):
        if s not in index:
            if len(states) >= STATE_CAP:
                raise RuntimeError("INCONCLUSIVE: independent quotient-state cap reached")
            index[s] = len(states)
            states.append(s)
        return index[s]

    for seed in seeds:
        add(seed.state)
    cursor = 0
    while cursor < len(states):
        s = states[cursor]
        assert min(s) >= 0 and gcd(*s) == 1
        assert abs(g.determinant(s)) == B * B
        assert not (s[0] >= s[2] and s[1] >= s[3])
        assert not (s[2] >= s[0] and s[3] >= s[1])
        row = []
        for digit in (1, 2):
            target, weight = g.consume(s, digit)
            other, other_weight = g.consume(swap(s), digit)
            assert canon(target) == canon(other) and weight == other_weight
            row.append((add(canon(target)), weight))
        edges.append(row)
        cursor += 1
    n = len(states)
    assert n == 11420 and 2 * n == 22840

    # Work-list relaxation is independent of the first computation's ordered sweeps.
    infinity = 10 ** 100
    distance = [infinity] * n
    queue, queued = deque(), set()
    for seed in seeds:
        i = index[seed.state]
        if seed.weight < distance[i]:
            distance[i] = seed.weight
            if i not in queued:
                queue.append(i)
                queued.add(i)
    pops = 0
    while queue:
        i = queue.popleft()
        queued.remove(i)
        pops += 1
        if pops > RELAXATION_CAP:
            raise RuntimeError("INCONCLUSIVE: independent shortest-path work cap reached")
        for j, weight in edges[i]:
            if distance[i] + weight < distance[j]:
                distance[j] = distance[i] + weight
                if j not in queued:
                    queue.append(j)
                    queued.add(j)
    assert all(value < infinity for value in distance)
    seed_slacks = [s.weight - distance[index[s.state]] for s in seeds]
    edge_slacks = [distance[i] + w - distance[j]
                   for i, row in enumerate(edges) for j, w in row]
    assert min(seed_slacks) == 0 and min(edge_slacks) == 0
    terminals = {i: g.finish(s) for i, s in enumerate(states)
                 if gcd(*g.append_terminal(s)) == B}
    terminal_slacks = [distance[i] + w for i, w in terminals.items()]
    assert len(terminals) == 116 and min(terminal_slacks) == 0

    # Check each actual integer, not only its fingerprint. Every successor was
    # generated above, so this also checks closure on exactly the listed set.
    order = sorted(range(n), key=lambda i: positive(states[i]))
    full_table = [[list(positive(states[i])), distance[i]] for i in order]
    assert table == full_table, "full state/potential table mismatch"
    graph_record = [[list(positive(states[i])), digit, list(positive(states[j])), w, w + 2 * digit]
                    for i in order for digit, (j, w) in enumerate(edges[i], 1)]
    graph_hash, potential_hash = fingerprint(graph_record), fingerprint(full_table)
    assert graph_hash == reference["graph_sha256"]
    assert potential_hash == reference["potential_sha256"]

    tight = [[(digit, j) for digit, (j, w) in enumerate(row, 1)
              if distance[i] + w == distance[j]] for i, row in enumerate(edges)]
    adjacency = [[v for _, v in row] for row in tight]
    reverse = [[] for _ in states]
    for i, row in enumerate(adjacency):
        for j in row:
            reverse[j].append(i)
    endpoints = {i for i, w in terminals.items() if distance[i] + w == 0}
    tight_seeds = [s for s in seeds if s.weight == distance[index[s.state]]]
    live = (reachable((index[s.state] for s in tight_seeds), adjacency)
            & reachable(endpoints, reverse))
    live_seeds = [s for s in tight_seeds if index[s.state] in live]
    assert len(live) == 6 and len(live_seeds) == 1
    start = index[live_seeds[0].state]
    assert (live_seeds[0].reflection, live_seeds[0].optional,
            live_seeds[0].prefix, live_seeds[0].weight) == (0, 1, (), -1)
    ends = endpoints & live
    assert len(ends) == 1
    terminal = next(iter(ends))
    live_rows = [[list(positive(states[i])), distance[i], terminals.get(i)] for i in live]
    assert sorted(live_rows) == sorted(reference["live_states"])
    live_edges = [[list(positive(states[i])), digit, list(positive(states[j])),
                   edges[i][digit - 1][1]]
                  for i in live for digit, j in tight[i] if j in live]
    assert sorted(live_edges) == sorted(reference["live_edges"])
    assert len(live_edges) == 7

    # Cutting at the unique accepting vertex gives a finite first-return DAG.
    # This proves the entire regular language without a bounded word sample.
    live_adj = {i: [(d, j) for d, j in tight[i] if j in live] for i in live}
    assert live_adj[start] == [(1, terminal)]
    first_returns = []
    todo = [(terminal, (), frozenset())]
    while todo:
        i, word, seen = todo.pop()
        assert i not in seen, "nonterminal cycle in the first-return graph"
        for d, j in live_adj[i]:
            extended = word + (d,)
            if j == terminal:
                first_returns.append(extended)
            else:
                todo.append((j, extended, seen | {i}))
    tokens = ((1, 2, 1, 1, 1), (2, 1, 1, 1, 1))
    assert sorted(first_returns) == list(tokens)

    features = {start: -2}  # Q-5P for optional 1 and empty paired prefix.
    queue = deque((start,))
    while queue:
        i = queue.popleft()
        for digit, j in live_adj[i]:
            value = features[i] + 2 - 10 * (digit == 2)
            if j in features:
                assert features[j] == value
            else:
                features[j] = value
                queue.append(j)
    assert features[terminal] == 0

    def block(word):
        matrix = g.IDENTITY
        for digit in word:
            matrix = g.multiply(matrix, g.BLOCKS[digit])
        return matrix

    q0 = g.multiply(g.multiply(g.P1, g.P1), g.P1)
    kernel = (47, 84, 68, -47)
    final = (89, 68, 0, -89)
    assert g.multiply(reflections[0], q0) == g.multiply(q0, kernel)
    assert g.multiply(kernel, g.P2) == g.multiply(g.P2, final)
    for first, second in (tokens, tokens[::-1]):
        assert g.multiply(kernel, block(first)) == g.multiply(block(second), kernel)
    assert g.append_terminal(q0) == (8, 5)
    assert g.multiply(final, final) == (B * B, 0, 0, B * B)
    for key, value in {"B": B, "states": n, "edges": 2 * n,
                       "seeds": len(seeds), "direct": 0,
                       "primitive_terminals": 116, "minimum": 0,
                       "prefix_visits": visits, "prefix_depth": depth}.items():
        assert reference[key] == value
    result = {"states": n, "edges": 2 * n, "seeds": len(seeds),
              "direct": 0, "primitive_terminals": len(terminals), "minimum": 0,
              "relaxation_steps": pops, "live_equality_states": len(live),
              "first_returns": ["".join(map(str, w)) for w in sorted(first_returns)],
              "terminal_Q_minus_5P": features[terminal],
              "graph_sha256": graph_hash, "potential_sha256": potential_hash}
    return result, states, tight, depth


def critical(states, adjacency, reference):
    n = len(states)
    number, low, component = [-1] * n, [-1] * n, [-1] * n
    active, active_set, groups = [], set(), []
    counter = 0
    for root in range(n):
        if number[root] >= 0:
            continue
        number[root] = low[root] = counter
        counter += 1
        active.append(root)
        active_set.add(root)
        stack = [(root, iter(adjacency[root]), None)]
        while stack:
            u, iterator, parent = stack[-1]
            nxt = next(iterator, None)
            if nxt is not None:
                _, v = nxt
                if number[v] < 0:
                    number[v] = low[v] = counter
                    counter += 1
                    active.append(v)
                    active_set.add(v)
                    stack.append((v, iter(adjacency[v]), u))
                elif v in active_set:
                    low[u] = min(low[u], number[v])
                continue
            stack.pop()
            if parent is not None:
                low[parent] = min(low[parent], low[u])
            if low[u] == number[u]:
                group = []
                while True:
                    v = active.pop()
                    active_set.remove(v)
                    component[v] = len(groups)
                    group.append(v)
                    if v == u:
                        break
                groups.append(group)
    assert counter == n and not active and not active_set
    assert len(groups) == 11411
    dag = [set() for _ in groups]
    weights = [Fraction(0) for _ in groups]
    critical_rows = []
    for c, group in enumerate(groups):
        for u in group:
            for _, v in adjacency[u]:
                if component[v] != c:
                    dag[c].add(component[v])
        if len(group) == 1 and not any(v == group[0] for _, v in adjacency[group[0]]):
            continue
        values = {group[0]: (0, 0)}
        todo, relations = [group[0]], []
        while todo:
            u = todo.pop()
            length, twos = values[u]
            for digit, v in adjacency[u]:
                if component[v] != c:
                    continue
                proposal = length + 1, twos + (digit == 2)
                if v not in values:
                    values[v] = proposal
                    todo.append(v)
                else:
                    relations.append((proposal[0] - values[v][0], proposal[1] - values[v][1]))
        assert set(values) == set(group)
        theta = next(Fraction(m, length) for length, m in relations if length)
        assert all(theta * length == m for length, m in relations)
        potential = {u: Fraction(m) - theta * length for u, (length, m) in values.items()}
        assert all(potential[v] - potential[u] == (digit == 2) - theta
                   for u in group for digit, v in adjacency[u] if component[v] == c)
        oscillation = max(potential.values()) - min(potential.values())
        weights[c] = oscillation + 1
        offset = potential[min(group, key=lambda u: positive(states[u]))]
        rows = sorted([[list(positive(states[u])), str(potential[u] - offset)] for u in group])
        edges = sorted([[list(positive(states[u])), digit, list(positive(states[v]))]
                        for u in group for digit, v in adjacency[u] if component[v] == c])
        critical_rows.append({"size": len(group), "theta": str(theta),
                              "oscillation": str(oscillation), "states": rows, "edges": edges})
    critical_rows.sort(key=lambda row: row["states"])
    expected_rows = sorted(reference["cyclic_components"], key=lambda row: row["states"])
    assert critical_rows == expected_rows, "critical states, coboundaries or edges differ"

    # The strong-component search emits sinks first. The budget recurrence is
    # the reverse of the first computation's forward longest-path calculation.
    assert all(v < u for u, row in enumerate(dag) for v in row)
    budget = []
    for u, row in enumerate(dag):
        budget.append(weights[u] + max([Fraction(0), *(1 + budget[v] for v in row)]))
    maximum = max(budget)
    assert maximum == 58 and str(maximum) == reference["L"]
    assert str(maximum + 1) == reference["K_graph"] == "59"
    zero_edges = sum(map(len, adjacency))
    assert zero_edges == reference["zero_edges"] == 12684
    assert reference["components"] == len(groups)
    assert sorted((row["size"], row["theta"], row["oscillation"]) for row in critical_rows) == sorted(
        [(6, "1/5", "1"), (3, "1/3", "2/3"), (3, "1/3", "2/3"),
         (1, "0", "0"), (1, "0", "0"), (1, "1", "0"), (1, "1", "0")])
    return {"zero_edges": zero_edges, "components": len(groups),
            "critical_components": [[r["size"], r["theta"], r["oscillation"]]
                                    for r in critical_rows],
            "L": str(maximum), "K_graph": 59}


def main():
    limits()
    sys.path.insert(0, str(HERE / "code"))
    import unit_reducer as g

    reference = read_data("reference-output.json")
    table = read_data("axis89-potential.json")
    assert reference["format_version"] == 1 and reference["status"] == "PASS"
    endpoint_record, states, tight, depth = endpoint(g, reference["endpoint"], table)
    critical_record = critical(states, tight, reference["critical"])
    assert depth == reference["critical"]["maximum_seed_prefix"] == 2
    assert reference["balanced_count_bound"] == "min_theta |m-theta*n| < 59*(Delta+1)+2"
    assert reference["genuine_count_bound"] == "min_theta |P-theta*Q| <= 118*Delta+121"
    result = {"status": "PASS", "B": B, "method": "unit row subtractions, work-list relaxation, iterative strong components, sink-first budgets",
              "endpoint": endpoint_record, "critical": critical_record,
              "prefix_allowance": depth, "full_potential_entries_compared": len(table)}
    recorded = HERE / "data" / "independent-reference.json"
    rendered = json.dumps(result, sort_keys=True, indent=2) + "\n"
    if WRITE_REFERENCE:
        recorded.write_text(rendered, encoding="utf-8")
    else:
        if not recorded.is_file():
            raise SystemExit("Missing data/independent-reference.json; the package is incomplete.")
        if json.loads(recorded.read_text(encoding="utf-8")) != json.loads(rendered):
            raise SystemExit("The values of the second computation differ from the recorded ones.")
    print(json.dumps(result, sort_keys=True))
    signal.alarm(0)


if __name__ == "__main__":
    main()
