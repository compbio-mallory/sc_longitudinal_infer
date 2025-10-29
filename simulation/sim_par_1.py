#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse, random, copy, sys
from collections import defaultdict, OrderedDict, deque
import numpy as np

# ------------------------------- Data containers -------------------------------

class Edge:
    def __init__(self, p, c):
        self.p = p
        self.c = c

class Node:
    def __init__(self, e, p, c):
        self.e = e
        self.p = p
        self.c = c  # "NA" or "id1;id2;..."

# ------------------------------- Tree I/O -------------------------------

def read_tree(tree_f):
    """
    Read gen_tree.py output with columns:
      TP  ID  PID  EL  Perc  Leaf  Mut
    Returns:
      parent:    {child -> parent}
      children:  {parent -> [child,...]}
      edge_len:  {(pid,cid) -> float}
      nodes:     set(node IDs)
      root_id:   id with parent -1
      node_info: {id -> dict(tp, pid, el, perc, leaf, mutflag)}
    """
    parent = {}
    children = defaultdict(list)
    edge_len = {}
    nodes = set()
    node_info = {}
    root_id = None

    with open(tree_f) as f:
        _ = f.readline()
        for line in f:
            if not line.strip(): continue
            tp_s, cid_s, pid_s, el_s, perc_s, leaf_s, mut_s = line.strip().split("\t")
            cid = int(cid_s); pid = int(pid_s)
            parent[cid] = pid
            nodes.add(cid)
            node_info[cid] = {
                "tp": float(tp_s),
                "pid": pid,
                "el": float(el_s),
                "perc": float(perc_s),
                "leaf": int(leaf_s),
                "mutflag": (mut_s.lower()=="true" if isinstance(mut_s,str) else bool(mut_s))
            }
            if pid != -1:
                children[pid].append(cid)
                edge_len[(pid, cid)] = float(el_s)
            else:
                root_id = cid

    return parent, children, edge_len, nodes, root_id, node_info

def compute_ancestors(parent, nodes):
    """ ancestors[v] = set of ancestors (excluding v) """
    ancestors = {}
    for v in nodes:
        s = set()
        u = parent.get(v, -1)
        while u != -1 and u in parent:
            s.add(u)
            u = parent[u]
        ancestors[v] = s
    return ancestors

def incomparable(a, b, ancestors):
    """ True if a and b are not in ancestor/descendant relation """
    if a == b: return False
    if a in ancestors.get(b,set()): return False
    if b in ancestors.get(a,set()): return False
    return True

# ------------------------------- Small utils -------------------------------

def init(n, m):
    return [[0]*m for _ in range(n)]

def print_matrix(M, n, matrix_file):
    with open(matrix_file, "w") as f:
        for i in range(n):
            f.write("\t".join(str(x) for x in M[i]) + "\n")

def count_total_value(M, n, m, value):
    return sum(1 for i in range(n) for j in range(m) if M[i][j] == value)

def print_dict(out_dict, out_file):
    with open(out_file, 'w') as f:
        for key,val in out_dict.items():
            f.write("{}\t{}\n".format(key, val))

# ------------------------------- Cells (legacy) -------------------------------

def sampleCells(noOfTimepoints, rng_py):
    # Legacy feel
    candidates = [100,300,600,1000]
    assigned = []
    for _ in range(1, int(noOfTimepoints)):
        assigned.append(rng_py.choice(candidates))
    return assigned

def distribute_SNVcells(cell_n, tree_f, out_f):
    """
    For each timepoint (integer >=2), assign cells to nodes at that TP
    proportionally to 'Perc'. Legacy behavior preserved.
    """
    with open(tree_f,"r") as file:
        _ = file.readline()
        line = file.readline().rstrip('\n')
        tp_nodes = OrderedDict()
        while(line != ""):
            line_a = line.split('\t')
            timepoint = float(line_a[0])
            if timepoint == 1.0 or not timepoint.is_integer():
                line = file.readline().rstrip('\n')
                continue
            nodeID = int(line_a[1])
            perc = float(line_a[4])
            tp_nodes.setdefault(timepoint, []).append(f"{nodeID};{perc}")
            line = file.readline().rstrip('\n')

    SNVcell_dict = {}
    prev_cell = -1
    i = 0
    for tp, items in tp_nodes.items():
        cell_num = int(cell_n[i]); i += 1
        total = 0
        last_node = None
        for s in items:
            nodeID_s, perc_s = s.split(";")
            nodeID = int(nodeID_s); perc = float(perc_s)
            assign = int(round(cell_num * perc))
            if total + assign > cell_num:
                assign = cell_num - total
            if assign == 0:
                prev_node = nodeID - 1
                key_prev = f"{int(tp)}_{prev_node}"
                prev_cells = SNVcell_dict[key_prev].split(";")
                SNVcell_dict[key_prev] = ";".join(prev_cells[:-1])
                SNVcell_dict[f"{int(tp)}_{nodeID}"] = prev_cells[-1]
                assign = 1
                total += 1
                last_node = nodeID
                continue
            if assign > 0:
                cell_ids = str(prev_cell + 1)
                for j in range(prev_cell + 2, prev_cell + assign + 1):
                    cell_ids += ";" + str(j)
                prev_cell += assign
                SNVcell_dict[f"{int(tp)}_{nodeID}"] = cell_ids
                last_node = nodeID
                total += assign
        if total < cell_num and last_node is not None:
            diff = cell_num - total
            key_last = f"{int(tp)}_{last_node}"
            ids = SNVcell_dict[key_last].split(";")
            last_id = int(ids[-1])
            ext = ";".join(str(k) for k in range(last_id+1, last_id+diff+1))
            SNVcell_dict[key_last] = SNVcell_dict[key_last] + (";" + ext if ext else "")
            prev_cell = last_id + diff

    print_dict(SNVcell_dict, out_f)
    return SNVcell_dict, prev_cell + 1

# ------------------------------- Baseline mutations -------------------------------

def draw_base_mutations(children, edge_len, mc, rng_np, node_info):
    """
    For each edge pid->cid:
      if Mut==True on child cid:
         draw K ~ Poisson( mc * EL(pid,cid) ) and assign those new IDs to child.
    Returns:
      node_new: {node -> [new ids introduced ON incoming edge]}
      next_id : next fresh mutation id
    """
    node_new = defaultdict(list)
    next_id = 0
    for pid, kids in children.items():
        for cid in kids:
            if not node_info[cid]["mutflag"]:
                continue
            lam = float(mc) * float(edge_len.get((pid, cid), 1.0))
            k = rng_np.poisson(lam) if lam > 0 else 0
            if k > 0:
                ids = list(range(next_id, next_id + k))
                node_new[cid].extend(ids)
                next_id += k
    return node_new, next_id

# ------------------------------- Parallel: helpers -------------------------------

def pick_incomparable_extras(all_nodes, ancestors, already, need_k, rng_np):
    """
    Pick up to need_k nodes that are incomparable to *every* node
    already chosen (including previously picked extras).
    Greedy but works well for our sizes.
    """
    if need_k <= 0:
        return []
    cand = list(all_nodes)
    rng_np.shuffle(cand)
    chosen = []                     # extras picked so far
    for v in cand:
        # compare against base origins + already-chosen extras
        if all(incomparable(v, u, ancestors) for u in (already + chosen)):
            chosen.append(v)
            if len(chosen) == need_k:
                break
    return chosen

def inject_new_mutations(n_new, p_extra, non_root_nodes, ancestors, rng_np, node_new, next_id):
    """
      - Create n_new NEW mutation IDs.
      - For each new mut:
          * choose 1 random origin edge (child node)
          * then add up to p_extra additional origins that are pairwise incomparable
            to ALL previously chosen origins for this mutation.
      - Update node_new for every origin.
      - Return parallel_origins ONLY for muts that have >=2 origins.
    """
    parallel_origins = defaultdict(set)

    all_nodes = list(non_root_nodes)
    if not all_nodes or n_new <= 0:
        return parallel_origins, next_id, []

    new_ids = []
    for _ in range(n_new):
        base = rng_np.choice(all_nodes)
        origins = [base]
        extras = pick_incomparable_extras(all_nodes, ancestors, origins, p_extra, rng_np)
        origins.extend(extras)

        mid = next_id; next_id += 1
        new_ids.append(mid)
        for node in origins:
            node_new[node].append(mid)
        if len(origins) >= 2:
            parallel_origins[mid] = set(origins)

    return parallel_origins, next_id, new_ids

# ------------------------------- Back: helpers -------------------------------

def compute_descendants(children, nodes):
    """ descendants[v] = set of descendants (excluding v) """
    descendants = {v:set() for v in nodes}
    # BFS from each node
    for v in nodes:
        q = deque(children.get(v, []))
        while q:
            u = q.popleft()
            descendants[v].add(u)
            for w in children.get(u, []):
                q.append(w)
    return descendants

def pick_losses_in_subtree(origin_child, back_k, descendants, ancestors, rng_np):
    """
    Choose up to back_k child-nodes (edges) under origin_child such that all selected
    losses are pairwise incomparable (i.e., no ancestor/descendant among them).
    """
    subs = list(descendants.get(origin_child, set()))
    if back_k <= 0 or not subs:
        return []
    rng_np.shuffle(subs)
    chosen = []
    for v in subs:
        if all(incomparable(v, u, ancestors) for u in chosen):
            chosen.append(v)
            if len(chosen) == back_k:
                break
    return chosen

def apply_back_losses(node_all, children, back_origin_losses, parent):
    """
    Remove mutation mid from all nodes in the subtrees rooted at each loss child.
    back_origin_losses: { mid -> (origin_child, [loss_child,...]) }
    """
    # map child -> list of all descendants including itself
    def subtree_nodes(start):
        out = []
        q = deque([start])
        while q:
            u = q.popleft()
            out.append(u)
            for v in children.get(u, []):
                q.append(v)
        return out

    for mid, (origin_child, losses) in back_origin_losses.items():
        for loss_child in losses:
            for v in subtree_nodes(loss_child):
                if mid in node_all.get(v, set()):
                    node_all[v].discard(mid)

# ------------------------------- Propagation -------------------------------

def propagate_to_descendants(root_id, children, node_new):
    """
    node_all[v] = all mutations present at node v after inheritance.
    node_new[v] are the new ones introduced on edge(parent->v).
    """
    node_all = defaultdict(set)
    q = deque([root_id])
    while q:
        u = q.popleft()
        for v in children.get(u, []):
            s = set(node_all[u])
            s.update(node_new.get(v, []))
            node_all[v] = s
            q.append(v)
    return node_all

# ------------------------------- Write outputs -------------------------------

def write_mut(prefix, node_all):
    out_f = prefix + ".mut.csv"
    with open(out_f, "w") as f:
        for nid in sorted(node_all.keys()):
            if node_all[nid]:
                f.write(f"{nid}\t{';'.join(str(x) for x in sorted(node_all[nid]))}\n")
    return out_f

def write_parallel(prefix, parallel_origins):
    out_f = prefix + ".parallel.csv"
    with open(out_f, "w") as f:
        for mid in sorted(parallel_origins.keys()):
            nodes = sorted(parallel_origins[mid])
            if nodes:
                f.write(f"{mid}\t{';'.join(str(n) for n in nodes)}\n")
    return out_f

def write_back(prefix, back_set, back_origin_losses):
    # list of back mutation IDs (one per line)
    out_ids = prefix + ".back.csv"
    with open(out_ids, "w") as f:
        for mid in sorted(back_set):
            f.write(f"{mid}\n")
    # mid \t origin_child \t loss1;loss2;...
    out_map = prefix + ".back_origin.csv"
    with open(out_map, "w") as f:
        for mid in sorted(back_origin_losses.keys()):
            origin_child, losses = back_origin_losses[mid]
            loss_str = ";".join(str(x) for x in sorted(losses)) if losses else ""
            f.write(f"{mid}\t{origin_child}\t{loss_str}\n")
    return out_ids, out_map

# ------------------------------- G / D matrix (legacy) -------------------------------

def add_missing(G, n, m, missingP, rng_py):
    total = n * m
    k = int(total * missingP)
    idxs = rng_py.sample(range(total), k)
    for i in idxs:
        r = i // m
        c = i % m
        G[r][c] = 3
    return G

def add_FPFNs(D, n, m, alpha, beta, rng_py):
    zero_positions = [(i,j) for i in range(n) for j in range(m) if D[i][j] == 0]
    one_positions  = [(i,j) for i in range(n) for j in range(m) if D[i][j] == 1]
    k_fp = int(len(zero_positions) * alpha)
    k_fn = int(len(one_positions)  * beta)
    for (i,j) in rng_py.sample(zero_positions, k_fp):
        D[i][j] = 1
    for (i,j) in rng_py.sample(one_positions,  k_fn):
        D[i][j] = 0
    return D

def mutation_matrix(node_all, SNVcell_dict, tree_f, n, m, missingP, alpha, beta,
                    G_matrix_file, D_matrix_file, D_miss_file, rng_py):
    # Build per-cell mutation lists
    SNVcell_mut = {}
    for tpnode, cells in SNVcell_dict.items():
        nodeID = int(tpnode.split("_")[1])
        muts_here = sorted(list(node_all.get(nodeID, set())))
        for cell in cells.split(";"):
            if cell == "": continue
            SNVcell_mut[cell] = muts_here

    G = init(n, m)
    for cellID_s, muts in SNVcell_mut.items():
        r = int(cellID_s)
        for mid in muts:
            if mid < m:
                G[r][mid] = 1

    print_matrix(G, n, G_matrix_file)
    total_zeros = count_total_value(G, n, m, 0)
    total_ones  = count_total_value(G, n, m, 1)
    print("total zeros for G: {}; total ones for G: {}".format(total_zeros, total_ones))

    D_miss = add_missing(copy.deepcopy(G), n, m, missingP, rng_py)
    print_matrix(D_miss, n, D_miss_file)
    total_zeros = count_total_value(D_miss, n, m, 0)
    total_ones  = count_total_value(D_miss, n, m, 1)
    print("total zeros after missing: {}; total ones after missing: {}".format(total_zeros, total_ones))

    D = add_FPFNs(copy.deepcopy(D_miss), n, m, alpha, beta, rng_py)
    total_zeros = count_total_value(D, n, m, 0)
    total_ones  = count_total_value(D, n, m, 1)
    print("total zeros after FPFN: {}; total ones after FPFN: {}".format(total_zeros, total_ones))

    print_matrix(D, n, D_matrix_file)
    return G

# ------------------------------- Main -------------------------------

def main():
    if len(sys.argv) <= 1:
        print("""
Generate simulated data with controlled NEW parallel & back mutations.

Usage: python sim_par_1.py -f tree.csv -P prefix [options]

Parallel knobs:
  -n/--new-muts INT           number of NEW mutations to add [3]
  -p/--parallel-edges INT     EXTRA parallel origins per new mut (total = 1+p) [0]
  --par-pct FLOAT             instead of -n, use percentage of baseline muts (rounded)

Back knobs (operate on NEW mutations only):
  --back-n INT                number of NEW muts to back-mutate [0]
  --back-k INT                max pairwise-incomparable loss edges per back mut [0]
  --back-pct FLOAT            instead of --back-n, use percentage of baseline muts

these are kept as it is :
  -a --alpha, -b --beta, -m --missing_rate, -t --timepoints, -mc --mut_const
  -G --Gmatrix, -cellstp --cellstp, -FPFN --fpfn
  -L (deprecated alias for -p)
  --seed INT
""")
        sys.exit(0)

    parser = argparse.ArgumentParser(description='Generate simulated data with controlled NEW parallel & back mutations.')
    parser.add_argument('-a', '--alpha', type=float, default=0.01)
    parser.add_argument('-b', '--beta', type=float, default=0.2)
    parser.add_argument('-m', '--missing_rate', type=float, default=0.2)
    parser.add_argument('-t','--timepoints', type=int, default=3)
    parser.add_argument('-FPFN','--fpfn', default=False)
    parser.add_argument('-mc', '--mut_const', type=float, default=1.7)
    parser.add_argument('-f', '--tree_file', required=True)
    parser.add_argument('-P', '--prefix', required=True)
    parser.add_argument('-G', '--Gmatrix', default="NA")
    parser.add_argument('-cellstp', '--cellstp', default="NA")
    # parallel knobs
    parser.add_argument('-n', '--new-muts', type=int, default=3, dest='new_muts')
    parser.add_argument('-p', '--parallel-edges', type=int, default=0, dest='parallel_edges')
    parser.add_argument('--par-pct', type=float, default=None, dest='par_pct')
    # back knobs
    parser.add_argument('--back-n', type=int, default=0, dest='back_n')
    parser.add_argument('--back-k', type=int, default=0, dest='back_k')
    parser.add_argument('--back-pct', type=float, default=None, dest='back_pct')
    # seed
    parser.add_argument('--seed', type=int, default=1)

    args = parser.parse_args()


    alpha = float(args.alpha)
    beta  = float(args.beta)
    missing = float(args.missing_rate)
    mc = float(args.mut_const)
    tree_f = args.tree_file
    prefix = args.prefix
    n_new = int(args.new_muts)
    p_extra = int(args.parallel_edges)
    back_n = int(args.back_n)
    back_k = int(args.back_k)
    par_pct = args.par_pct
    back_pct = args.back_pct

    # RNGs compatible with old NumPy installs
    rng_np = np.random.RandomState(args.seed)
    class _PyRNG:
        def __init__(self, seed):
            random.seed(seed)
        def sample(self, population, k):
            return random.sample(population, k)
        def choice(self, arr):
            return random.choice(arr)
    rng_py = _PyRNG(args.seed)

    print("For file --> ", tree_f)

    # Read tree
    parent, children, edge_len, nodes, root_id, node_info = read_tree(tree_f)
    non_root_nodes = set(n for n in nodes if parent.get(n, -1) != -1)

    # 1) Baseline new mutations per Mut flag
    node_new, next_id = draw_base_mutations(children, edge_len, mc, rng_np, node_info)

    # If pct knobs provided, map to counts using baseline total
    baseline_total = (max((max(s) if s else -1) for s in node_new.values()) + 1) if node_new else 0
    if par_pct is not None:
        n_new = max(0, int(round(baseline_total * float(par_pct) / 100.0)))
    if back_pct is not None:
        back_n = max(0, int(round(baseline_total * float(back_pct) / 100.0)))

    # 2) Inject NEW mutations with p-extra parallel edges each
    ancestors = compute_ancestors(parent, nodes)
    parallel_origins, next_id, new_ids = inject_new_mutations(
        n_new=n_new,
        p_extra=p_extra,
        non_root_nodes=non_root_nodes,
        ancestors=ancestors,
        rng_np=rng_np,
        node_new=node_new,
        next_id=next_id
    )

    # 3) Propagate (pre-back) to get full sets at nodes
    node_all = propagate_to_descendants(root_id, children, node_new)

    # 4) Back-mutations **only on NEW IDs** (subset of new_ids)
    back_set = set()
    back_origin_losses = {}  # mid -> (origin_child, [loss_child,...])

    if back_n > 0 and new_ids:
        # pick up to back_n from new_ids (deterministic-ish via RNG)
        ids_pool = list(new_ids)
        rng_np.shuffle(ids_pool)
        pick = ids_pool[:min(back_n, len(ids_pool))]

        # For each chosen mid, decide an origin (where it was introduced) and pick losses
        # We map mid -> set(origin nodes) from node_new (origins are the nodes that contain mid but NOT its parent)
        # construct origin lookup
        origins_by_mid = defaultdict(list)
        for v, ids in node_new.items():
            for mid in ids:
                origins_by_mid[mid].append(v)

        descendants = compute_descendants(children, nodes)

        for mid in pick:
            # choose one origin_child among origins (there is always >=1)
            olist = origins_by_mid.get(mid, [])
            if not olist:
                continue
            origin_child = rng_np.choice(olist)

            # choose up to back_k incomparable losses inside the subtree of origin_child
            losses = pick_losses_in_subtree(origin_child, back_k, descendants, ancestors, rng_np)

            back_set.add(mid)
            back_origin_losses[mid] = (origin_child, losses)

        # Apply removals
        apply_back_losses(node_all, children, back_origin_losses, parent)

    # 5) Write mut/parallel/back (ground truth)
    mut_f = write_mut(prefix, node_all)
    par_f = write_parallel(prefix, parallel_origins)
    back_ids_f, back_map_f = write_back(prefix, back_set, back_origin_losses)
    total_muts = len(set(m for s in node_all.values() for m in s))
    print("Total unique mutations (after parallel+back):", total_muts)
    print("Wrote:", mut_f)
    print("Wrote:", par_f)
    print("Wrote:", back_ids_f)
    print("Wrote:", back_map_f)

    # 6) Cells + G/D (legacy)
    cell_n = sampleCells(args.timepoints, rng_py)
    print(" Assigned cells ", cell_n)
    SNVcell_f = prefix + ".SNVcell.csv"
    SNVcell_dict, n_cells = distribute_SNVcells(cell_n, tree_f, SNVcell_f)

    m_mut = (max((max(s) if s else -1) for s in node_all.values()) + 1) if node_all else 0
    print(" Final total no of cells ", n_cells)
    print(" Total cells ", n_cells, " total mutations ", m_mut)

    Gf   = prefix + ".G.csv"
    Df   = prefix + ".D.csv"
    Miss = prefix + ".miss.csv"
    _ = mutation_matrix(node_all, SNVcell_dict, tree_f, n_cells, m_mut, missing, alpha, beta, Gf, Df, Miss, rng_py)

    print("="*108)

if __name__ == "__main__":
    main()
