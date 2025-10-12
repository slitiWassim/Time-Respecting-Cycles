import os 
import sys
import time

import pandas as pd
import numpy as np 

from raphtory import algorithms as rp
from raphtory import Graph


from collections import defaultdict
from bisect import bisect_right




def time_respecting_cycles(rolling_g, max_length=None, max_cycles=None, max_duration=None):
    """

    Enumerate all time-respecting simple cycles in a temporal graph with duration less 'max_duration'.

    Parameters
    ----------
    rolling_g : Raphtory Graph 
        Raphtory Temporal Knowledge Graph
    max_length : int, optional
        Maximum number of nodes in a cycle (to prune long cycles).
    max_cycles : int, optional
        Stop after this many cycles have been found (Prevent RAM overload).
    max_duration in (ms) : int , optional
        Only return cycles with duration less than max_duration
            
        
    Yields
    ------
    (cycle_nodes, cycle_timestamps)
        Lists of node names and edges timestamps in temporal chronological order.( Time Respecting Cycles)
    """

    cycle_count = [0]
    adjacency = defaultdict(dict)

    # ---- Build adjacency ----
    for node in rolling_g.nodes:
        for edge in node.out_edges:
            times = sorted(edge.history().tolist())
            if times:
                adjacency[node.name][edge.dst.name] = times

    def _out_neighbors(node_id, prev_time):
        for nbr, times in adjacency[node_id].items():
            if prev_time is None:
                for t in times:
                    yield nbr, t
            else:
                idx = bisect_right(times, prev_time)
                for t in times[idx:]:
                    yield nbr, t

    # ---- Tarjan’s SCC ----
    def strongly_connected_components():
        index, lowlink = {}, {}
        stack, on_stack = [], set()
        sccs, idx = [], 0

        def dfs(v):
            nonlocal idx
            index[v] = lowlink[v] = idx
            idx += 1
            stack.append(v)
            on_stack.add(v)
            for w in adjacency.get(v, {}).keys():
                if w not in index:
                    dfs(w)
                    lowlink[v] = min(lowlink[v], lowlink[w])
                elif w in on_stack:
                    lowlink[v] = min(lowlink[v], index[w])
            if lowlink[v] == index[v]:
                comp = set()
                while True:
                    w = stack.pop()
                    on_stack.remove(w)
                    comp.add(w)
                    if w == v:
                        break
                sccs.append(comp)

        for v in adjacency.keys():
            if v not in index:
                dfs(v)
        return sccs

    # ---- Johnson-like cycle search ----
    def _johnson_cycle_search(start):
        path, times_path = [start], []
        blocked = set()
        B = defaultdict(set)

        def unblock(u):
            if u in blocked:
                blocked.remove(u)
                for v in list(B[u]):
                    B[u].remove(v)
                    unblock(v)

        def backtrack(v, last_time):
            closed = False
            blocked.add(v)

            for w, t_next in _out_neighbors(v, last_time):
                # Duration pruning
                if times_path and max_duration and t_next - times_path[0] > max_duration:
                    continue

                if w == start:
                    duration = t_next - times_path[0] if times_path else 0
                    if (
                        (max_length is None or len(path) <= max_length)
                        and (max_duration is None or duration <= max_duration)
                    ):
                        cycle_count[0] += 1
                        yield (path[:] + [start], times_path[:] + [t_next])
                        closed = True
                        if max_cycles and cycle_count[0] >= max_cycles:
                            return
                elif w not in path and (max_length is None or len(path) < max_length):
                    path.append(w)
                    times_path.append(t_next)
                    for result in backtrack(w, t_next):
                        yield result
                        closed = True
                        if max_cycles and cycle_count[0] >= max_cycles:
                            return
                    path.pop()
                    times_path.pop()

            if closed:
                unblock(v)
            else:
                for w, _ in _out_neighbors(v, None):
                    B[w].add(v)

        yield from backtrack(start, None)

    # ---- Run SCCs ----
    sccs = strongly_connected_components()
    print(f"Found {len(sccs)} SCCs")
    components = [c for c in sccs if len(c) > 1]

    for comp in components:
        for start in comp:
            for cycle in _johnson_cycle_search(start):
                yield cycle
                if max_cycles and cycle_count[0] >= max_cycles:
                    return





