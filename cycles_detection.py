from collections import defaultdict
import numpy as np 
from itertools import chain
from raphtory import algorithms as rp
from raphtory import Graph


# ---- Main function ----
def time_respecting_cycles(
    rolling_g,
    max_length=None,
    max_cycles=None,
    max_duration=None,
    visited_memo=True
):
    """
    Enumerate **time-respecting cycles** in a Raphtory temporal directed graph.

    This algorithm identifies all elementary cycles where the sequence of edge timestamps
    is **strictly increasing** — ensuring that the traversal follows the arrow of time.
    It combines a structural depth-first cycle search (based on Johnson's algorithm [1]_)
    with efficient temporal validation using stored edge timestamp histories.

    The procedure builds preliminary (min, max) timestamp intervals for each edge to quickly
    prune impossible transitions.

    Parameters
    ----------
    rolling_g : Raphtory Graph
        A Raphtory temporal knowledge graph. 

    max_length : int, optional
        Maximum number of distinct nodes allowed in a cycle.
        Limits cycle depth to control computational growth.
        Default is ``None`` (no structural limit).

    max_cycles : int, optional
        Stop after this many **valid** cycles have been found.
        Prevents memory overload in dense graphs. Default is ``None``.

    max_duration : int, optional
        Maximum duration (in milliseconds) between the earliest and latest
        timestamps within a cycle. Cycles exceeding this temporal span are skipped.
        Default is ``None`` (no duration filter).

    visited_memo : bool, optional, default=True
        Enables memoization of previously visited node-interval states during
        backtracking. This dramatically reduces redundant recursion and speeds
        up traversal on large graphs.

    Yields
    ------
    (cycle_nodes, cycle_timestamps) : 
        - **cycle_nodes** : list of node names forming a cycle (start node repeated at end).
        - **cycle_timestamps** : list of timestamps (strictly increasing) associated with the edges.

        Each yielded pair corresponds to a time-respecting cycle

    Notes
    -----
    The algorithm proceeds in several conceptual stages:


    1. **Interval-based temporal adjacency (min, max) construction:**  
       For every edge, build `(t_min, t_max)` intervals 

    2. **SCC decomposition:**  
       Identify strongly connected components (SCCs) using Raphtory's built-in method.
       Following the implementation of NetworkX for detecting cycles in directed graphs,
       Johnson's algorithm is enhanced with well-known preprocessing techniques by restricting
       the search to strongly connected components. https://github.com/networkx/networkx/blob/main/networkx/algorithms/cycles.py

    3. **Johnson-style search:**  
       Perform a modified Johnson's backtracking algorithm within each SCC to
       enumerate structural cycles, pruned by temporal interval compatibility.
      
    4. **Cycle validation:**  
       Each structural cycle is checked via ``validate_cycle()`` to ensure that
       an increasing sequence of timestamps exists along the edges.

    5. **Yielding results:**  
       Only temporally valid cycles are yielded. Duration, length, and count
       constraints are enforced throughout.


    References
    ----------
    .. [1] Johnson, D. B. (1975). *Finding all the elementary circuits of a directed graph.*
           SIAM Journal on Computing, 4(1), 77-84. https://doi.org/10.1137/0204007

    Examples
    --------
    >>> cycles = time_respecting_cycles(rolling_g, max_length=4, max_cycles=500)
    >>> for nodes, times in cycles:
    ...     print(nodes, times)
    ('A', 'B', 'C', 'A') [102, 134, 180]
    """

    # ---- Strongly connected components ----
    components = [list(value.name) for _, value in rp.strongly_connected_components(rolling_g).groups() if len(value) > 1]
    
    ## only keep nodes that could be in a cycle
    nodes = list(chain.from_iterable(components))
    rolling_g = rolling_g.subgraph(nodes)

    cycle_count = [0]
    adjacency = defaultdict(dict)

    # ---- Build interval-based temporal adjacency (min, max) ----
    for node in rolling_g.nodes:
        for edge in node.out_edges:
            times = edge.history().tolist()
            if times is None or len(times) == 0:
                continue

            t_min, t_max = int(min(times)), int(max(times))
            adjacency[node.name][edge.dst.name] = {
                "interval": (t_min, t_max)
            }


    

    # ---- Main cycle search ----
    for comp in components:
        for start in comp:
            for raw_cycle_nodes,_ in johnson_cycle_search(start,
                                                          adjacency,
                                                          cycle_count,
                                                          max_length,
                                                          max_cycles,
                                                          max_duration,
                                                          visited_memo):
                
                edge_cycle = list(zip(raw_cycle_nodes[:-1], raw_cycle_nodes[1:]))
                valid_cycles = validate_cycle(rolling_g,edge_cycle)
                if valid_cycles:
                    cycle_count[0] += len(valid_cycles)
                    for vc in valid_cycles:
                        yield vc
                    if max_cycles and cycle_count[0] >= max_cycles:
                        return



# ---- Johnson cycle search ----
#   The depth-first traversal search that adapts Johnson's classical 
#   cycle enumeration to temporal graphs, pruning infeasible paths using 
#   timestamp interval compatibility checks.

def johnson_cycle_search(
    start,
    adjacency,
    cycle_count,
    max_length,
    max_cycles,
    max_duration,
    visited_memo
):
    path = [start]
    times_path = []
    blocked = set()
    B = defaultdict(set)
    visited_states = set()

    def unblock(u):
        if u in blocked:
            blocked.remove(u)
            for v in list(B[u]):
                B[u].remove(v)
                unblock(v)

    def backtrack(v, prev_interval):
        closed = False
        blocked.add(v)

        key = (v, None if prev_interval is None else (prev_interval[0], prev_interval[1], len(path)))
        if visited_memo and key in visited_states:
            blocked.remove(v)
            return
        if visited_memo:
            visited_states.add(key)

        for w, next_interval, _ in out_neighbors(v, prev_interval, adjacency):
            if max_duration and times_path:
                duration = next_interval[1] - times_path[0][0]
                if duration > max_duration:
                    continue

            if w == start:
                duration = next_interval[1] - times_path[0][0] if times_path else 0
                if ((max_length is None or len(path) <= max_length)
                    and (max_duration is None or duration <= max_duration)):
                    yield (path[:] + [start], times_path[:] + [next_interval])
                    closed = True
                    if max_cycles and cycle_count[0] >= max_cycles:
                        return
            elif w not in path and (max_length is None or len(path) < max_length):
                path.append(w)
                times_path.append(next_interval)
                for result in backtrack(w, next_interval):
                    yield result
                    closed = True
                    if max_cycles and cycle_count[0] >= max_cycles:
                        return
                path.pop()
                times_path.pop()

        if closed:
            unblock(v)
        else:
            for w_key in adjacency[v].keys():
                B[w_key].add(v)

    yield from backtrack(start, None)



# ---- Helper:  Fast temporal feasibility check for edge intervals. ----
def is_interval_compatible(prev_interval, next_interval):
    return next_interval[1] > prev_interval[0]


# ---- Generates feasible outgoing edges. ----
def out_neighbors(node_id, prev_interval, adjacency):
    for nbr, data in adjacency[node_id].items():
        next_interval = data["interval"]

        if prev_interval is None:
            yield nbr, next_interval, False
            continue

        if is_interval_compatible(prev_interval, next_interval):
            yield nbr, next_interval, False
            continue


# ---- validate_cycle : Checks whether an edge cycle satisfies time-respecting constraints. ----
def validate_cycle(rolling_g,edge_cycle):
    time_lists = []
    for src, dst in edge_cycle:
        edge = rolling_g.edge(src=src, dst=dst)
        if edge is None:
            return []
        times = edge.history().tolist()
        if times is None or len(times) == 0:
            return []
        time_lists.append(times)

    node_path = [edge_cycle[0][0]] + [dst for _, dst in edge_cycle]
    results = []

    # ---- Simple fallback ----
    def dfs(i, prev_t, combo):
        if i == len(time_lists):
            results.append((tuple(node_path), combo[:]))
            return
        for t in time_lists[i]:
            if t > prev_t:
                combo.append(t)
                dfs(i + 1, t, combo)
                combo.pop()

    dfs(0, float("-inf"), [])
    return results




