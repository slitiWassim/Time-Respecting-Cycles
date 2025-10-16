


## How the algorithm works 

### 1) **Interval-based temporal adjacency (min, max) construction**

For every edge in the graph, the algorithm builds an **interval representation** $(t_{min}, t_{max})$ summarizing the earliest and latest observed timestamps. This provides a **lightweight temporal abstraction** used for fast preliminary pruning before fine-grained timestamp validation.

**Why Building an interval-based temporal adjacency (min, max) ?** 
significantly improves performance by avoiding repeated per-edge timestamp scans.  Instead of recomputing the minimum and maximum times every time an edge is accessed, these values are pre-computed once and stored.

### 2) Strongly Connected Components decomposition

Following the [NetworkX implementation of Johnson’s algorithm](https://github.com/networkx/networkx/blob/main/networkx/algorithms/cycles.py) ,  the search is restricted to Strongly Connected Components .

Strongly Connected Components (SCCs) are identified using **Raphtory** built-in method.  

### 3) Johnson cycle search

Running the algorithm by examining cycles from every temporal edge individually would cause an exponential increase in complexity due to the vast number of temporal edges in temporal graphs. To improve efficiency, we instead focus on identifying **potentially time-ordered cycles**. This is achieved by locating structural cycles where consecutive edges $e_t$​ and $e_{t+1}$​ satisfy the temporal consistency condition  $min⁡(e_t.times)<max⁡(e_{t+1}.times)$.

This filtering step significantly reduces the search space while preserving cycles that are likely to be temporally valid.

Within each strongly connected component (SCC), a **modified Johnson backtracking algorithm** is then employed to enumerate structural cycles. Each candidate cycle undergoes **early pruning** through an interval compatibility check (is_interval_compatible),to ensure temporal feasibility between consecutive edges.

### 4) **Cycle validation**

Each candidate (structural) cycle is passed through `validate_cycle()`, which performs a fine grained temporal validation  by checking for a strictly **increasing sequence of timestamps** across edges.   Only cycles satisfying full temporal consistency are accepted.

Instead of computing the full **Cartesian product** of all timestamp combinations (which can explode combinatorially), it performs **incremental DFS pruning**:

- At each step, it picks the next timestamp that is strictly larger than the last chosen one.    
- It stops early when no valid next timestamp exists.

