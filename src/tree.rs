//! The taxonomy the EM estimates over
//!
//! The tree is a [`phylo`] [`SimpleRootedTree`]. Each node carries, in phylo's slots:
//! its *taxa* = the reported label (a taxid, the header of an index sequence, or `root`),
//! its *edge weight* = the conditional proportion `q_{c|v}` of the node given its parent
//! (the root keeps 1), its *zeta* = the expected read count `N_v` of its subtree, and its
//! *annotation* = the rank (`sequence` for an index record, `root` for the synthetic root).
//! Node ids are dense and never raw taxids: the arena is indexed by id.
//!
//! Before the EM: [`reference_tree`] builds the lineages of every reference taxid from an
//! NCBI `nodes.dmp` ([`TaxNodes`]) and hangs **every index sequence as its own leaf** under
//! its taxid, so taxids are internal vertices and the leaves are the sequences. Without a
//! nodes file the taxids hang under the root; without a taxid map the sequences do. After
//! alignment, [`contract_to_observed`] keeps the leaves that received a read together with
//! every vertex on their root paths and deletes the rest. A [`Partition`] says which leaf
//! owns each likelihood column (sequence): one sequence per leaf, until a collapse pools
//! the sequences beneath a vertex into it (the pooled leaf's likelihood is their mean).
//!
//! During the EM: [`init_conditionals`], [`leaf_marginals`] / [`contig_pi`],
//! [`subtree_counts`], [`prune_and_collapse`], [`m_step`] and [`penalty`] implement one
//! iteration; [`assign_read`] is the final descent to the deepest vertex whose posterior
//! meets the confidence level.

use crate::group::Grouping;
use crate::utils::EMProb;
use crate::SeqId;
use anyhow::{bail, Context, Result};
use phylo::node::NodeID;
use phylo::prelude::*;
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};

pub type TaxTree = SimpleRootedTree<String, EMProb, EMProb>;
pub type TaxNode = Node<String, EMProb, EMProb>;

/// Column marker for a sequence whose leaf was pruned.
pub const NO_NODE: NodeID = usize::MAX;
/// Id of the root (never deleted, never renumbered).
pub(crate) const ROOT: NodeID = 0;
const ROOT_RANK: &str = "root";
/// Rank annotation of an index sequence at a tip.
pub const SEQ_RANK: &str = "sequence";
const NO_RANK: &str = "-";

/// Parent and rank of every taxid in an NCBI `nodes.dmp` (or a `child parent [rank]` file).
#[derive(Debug, Default)]
pub struct TaxNodes {
    parent: HashMap<u32, u32>,
    rank: HashMap<u32, String>,
}

impl TaxNodes {
    pub fn from_file(path: &str) -> Result<Self> {
        let text = std::fs::read_to_string(path)
            .with_context(|| format!("failed to read taxonomy nodes '{}'", path))?;
        Self::parse(&text).with_context(|| format!("taxonomy nodes '{}'", path))
    }

    /// One node per line: `taxid | parent | rank | …` (`nodes.dmp`) or whitespace-separated
    /// `taxid parent [rank]`. Blank lines and lines starting with `#` are skipped. The root
    /// is the node whose parent is itself.
    pub fn parse(text: &str) -> Result<Self> {
        let mut parent = HashMap::new();
        let mut rank = HashMap::new();
        for (lineno, line) in text.lines().enumerate() {
            let trimmed = line.trim();
            if trimmed.is_empty() || trimmed.starts_with('#') {
                continue;
            }
            let fields: Vec<&str> = if trimmed.contains('|') {
                trimmed.split('|').map(str::trim).collect()
            } else {
                trimmed.split_whitespace().collect()
            };
            let (Some(t), Some(p)) = (fields.first(), fields.get(1)) else {
                bail!(
                    "line {}: expected `taxid parent [rank]`, got {:?}",
                    lineno + 1,
                    line
                );
            };
            let parse = |s: &str| -> Result<u32> {
                match s.parse::<u32>() {
                    Ok(v) if v > 0 => Ok(v),
                    _ => bail!(
                        "line {}: taxid '{}' is not a positive integer",
                        lineno + 1,
                        s
                    ),
                }
            };
            let (t, p) = (parse(t)?, parse(p)?);
            if let Some(&prev) = parent.get(&t) {
                if prev != p {
                    bail!(
                        "line {}: taxid {} has parents {} and {}",
                        lineno + 1,
                        t,
                        prev,
                        p
                    );
                }
            }
            parent.insert(t, p);
            let r = fields.get(2).copied().filter(|r| !r.is_empty());
            rank.insert(t, r.unwrap_or(NO_RANK).to_string());
        }
        Ok(Self { parent, rank })
    }

    pub fn len(&self) -> usize {
        self.parent.len()
    }

    pub fn is_empty(&self) -> bool {
        self.parent.is_empty()
    }

    pub fn contains(&self, taxid: u32) -> bool {
        self.parent.contains_key(&taxid)
    }

    /// The parent of a taxid; `None` at the taxonomy root or for an unknown taxid.
    pub fn parent(&self, taxid: u32) -> Option<u32> {
        self.parent.get(&taxid).copied().filter(|&p| p != taxid)
    }

    pub fn rank(&self, taxid: u32) -> &str {
        self.rank.get(&taxid).map_or(NO_RANK, String::as_str)
    }
}

pub(crate) fn new_node(id: NodeID, label: &str, rank: &str) -> TaxNode {
    let mut node = TaxNode::new(id);
    node.set_taxa(Some(label.to_string()));
    node.set_annotation(Some(rank.into()));
    node.set_weight(None);
    node.set_zeta(Some(0.0));
    node
}

pub(crate) fn new_tree() -> TaxTree {
    let mut tree = TaxTree::new(ROOT);
    tree.set_node(new_node(ROOT, "root", ROOT_RANK));
    tree
}

/// True for the synthetic root above the taxonomy (a read stopping there is unclassified).
pub fn is_synthetic_root(tree: &TaxTree, id: NodeID) -> bool {
    rank(tree, id) == ROOT_RANK
}

/// The label a node is reported under.
pub fn label(tree: &TaxTree, id: NodeID) -> &str {
    tree.get_node(id)
        .and_then(|n| n.get_taxa())
        .map_or("", String::as_str)
}

pub fn rank(tree: &TaxTree, id: NodeID) -> &str {
    tree.get_node(id)
        .and_then(|n| n.get_annotation())
        .unwrap_or(NO_RANK)
}

/// True for a leaf that is an index sequence (as opposed to a collapsed taxon).
pub fn is_sequence(tree: &TaxTree, id: NodeID) -> bool {
    rank(tree, id) == SEQ_RANK
}

fn weight(tree: &TaxTree, id: NodeID) -> EMProb {
    tree.get_node(id)
        .and_then(|n| n.get_weight())
        .unwrap_or(0.0)
}

/// `q_{v|parent}`, the share of its parent's mass a node receives (the root: 1).
pub fn conditional(tree: &TaxTree, id: NodeID) -> EMProb {
    if id == tree.get_root_id() {
        1.0
    } else {
        weight(tree, id)
    }
}

fn set_weight(tree: &mut TaxTree, id: NodeID, w: EMProb) {
    if let Some(n) = tree.get_node_mut(id) {
        n.set_weight(Some(w));
    }
}

/// `N_v`, the expected reads in the subtree of `id`.
pub fn count(tree: &TaxTree, id: NodeID) -> EMProb {
    tree.get_node(id).and_then(|n| n.get_zeta()).unwrap_or(0.0)
}

fn set_count(tree: &mut TaxTree, id: NodeID, n: EMProb) {
    if let Some(node) = tree.get_node_mut(id) {
        node.set_zeta(Some(n));
    }
}

pub(crate) fn children(tree: &TaxTree, id: NodeID) -> Vec<NodeID> {
    tree.get_node(id)
        .map(|n| n.get_children().to_vec())
        .unwrap_or_default()
}

/// Node ids in pre-order (parents before children).
pub fn preorder(tree: &TaxTree) -> Vec<NodeID> {
    tree.preord_ids(tree.get_root_id())
        .expect("the root is a node of its tree")
        .collect()
}

/// Node ids in post-order (children before parents), from `start`.
fn postorder(tree: &TaxTree, start: NodeID) -> Vec<NodeID> {
    tree.postord_ids(start)
        .expect("start is a node of the tree")
        .collect()
}

/// Length of a per-node array indexed by id.
pub(crate) fn arena_len(tree: &TaxTree) -> usize {
    tree.get_node_ids().max().map_or(0, |m| m + 1)
}

/// The reference taxonomy before alignment: every index sequence at a tip under its taxid.
pub struct RefTree {
    pub tree: TaxTree,
    /// Tip of each index sequence, by [`SeqId`].
    pub leaf_of_seq: Vec<NodeID>,
    /// Taxids of the map that the nodes file does not know (hung under the root).
    pub n_unknown: usize,
    /// Sequences without a taxid (hung under the root).
    pub n_unmapped: usize,
}

impl RefTree {
    pub fn n_nodes(&self) -> usize {
        self.tree.num_nodes()
    }

    /// Tips = index sequences.
    pub fn n_tips(&self) -> usize {
        self.leaf_of_seq.len()
    }
}

/// Build the tree over the index: the lineage of every taxid of `grouping` from `nodes`
/// (every taxid directly under the root when there is no nodes file), then one leaf per
/// index sequence under its taxid. Taxids absent from `nodes` hang under the root (one
/// warning each); sequences the map did not name are leaves under the root.
pub fn reference_tree(grouping: &Grouping, nodes: Option<&TaxNodes>) -> Result<RefTree> {
    let mut tree = new_tree();
    let mut next_id: NodeID = ROOT + 1;
    let mut by_taxid: HashMap<u32, NodeID> = HashMap::new();
    let mut unknown: HashSet<u32> = HashSet::new();

    let mut node_of_group = vec![ROOT; grouping.n_groups()];
    for (col, slot) in node_of_group.iter_mut().enumerate() {
        if let Some(taxid) = grouping.taxid_of(SeqId::new(col as u32)) {
            *slot = place_taxid(
                &mut tree,
                &mut next_id,
                &mut by_taxid,
                &mut unknown,
                nodes,
                taxid,
            )?;
        }
    }

    let mut leaf_of_seq = Vec::with_capacity(grouping.n_sequences());
    let mut n_unmapped = 0;
    for s in 0..grouping.n_sequences() {
        let seq = SeqId::new(s as u32);
        let parent = node_of_group[grouping.group(seq).index()];
        if parent == ROOT {
            n_unmapped += 1;
        }
        let id = next_id;
        next_id += 1;
        tree.add_child(parent, new_node(id, grouping.seq_name(seq), SEQ_RANK));
        leaf_of_seq.push(id);
    }

    Ok(RefTree {
        tree,
        leaf_of_seq,
        n_unknown: unknown.len(),
        n_unmapped,
    })
}

/// The node of `taxid`, inserting its missing ancestors top-down.
fn place_taxid(
    tree: &mut TaxTree,
    next_id: &mut NodeID,
    by_taxid: &mut HashMap<u32, NodeID>,
    unknown: &mut HashSet<u32>,
    nodes: Option<&TaxNodes>,
    taxid: u32,
) -> Result<NodeID> {
    if let Some(&id) = by_taxid.get(&taxid) {
        return Ok(id);
    }
    let known = match nodes {
        Some(n) if n.contains(taxid) => Some(n),
        Some(_) => {
            if unknown.insert(taxid) {
                eprintln!(
                    "Warning: taxid {} is not in the taxonomy nodes; placed directly under the root",
                    taxid
                );
            }
            None
        }
        None => None,
    };
    let Some(nodes) = known else {
        let id = *next_id;
        *next_id += 1;
        tree.add_child(ROOT, new_node(id, &taxid.to_string(), NO_RANK));
        by_taxid.insert(taxid, id);
        return Ok(id);
    };
    // Climb until a taxid that already has a node, or the taxonomy root.
    let mut lineage = vec![taxid];
    let mut attach = ROOT;
    let mut cur = taxid;
    while let Some(p) = nodes.parent(cur) {
        if let Some(&id) = by_taxid.get(&p) {
            attach = id;
            break;
        }
        if lineage.contains(&p) {
            bail!("taxonomy nodes contain a cycle through taxid {}", p);
        }
        lineage.push(p);
        cur = p;
    }
    for &t in lineage.iter().rev() {
        let id = *next_id;
        *next_id += 1;
        tree.add_child(attach, new_node(id, &t.to_string(), nodes.rank(t)));
        by_taxid.insert(t, id);
        attach = id;
    }
    Ok(attach)
}

/// Which leaf owns each likelihood column (sequence), and how many sequences each leaf holds
/// (`|P_ℓ|`: 1 for a sequence, the pooled set for a collapsed vertex).
#[derive(Clone, Debug)]
pub struct Partition {
    leaf_of_col: Vec<NodeID>,
    size: HashMap<NodeID, usize>,
}

impl Partition {
    pub fn new(leaf_of_col: Vec<NodeID>) -> Self {
        let mut p = Self {
            leaf_of_col,
            size: HashMap::new(),
        };
        p.recount();
        p
    }

    fn recount(&mut self) {
        self.size.clear();
        for &leaf in &self.leaf_of_col {
            if leaf != NO_NODE {
                *self.size.entry(leaf).or_insert(0) += 1;
            }
        }
    }

    pub fn n_cols(&self) -> usize {
        self.leaf_of_col.len()
    }

    /// The leaf of a column, or [`NO_NODE`] once pruned.
    pub fn leaf_of(&self, col: usize) -> NodeID {
        self.leaf_of_col[col]
    }

    /// `|P_ℓ|`: sequences (with reads) a leaf holds.
    pub fn size_of(&self, leaf: NodeID) -> usize {
        self.size.get(&leaf).copied().unwrap_or(0)
    }
}

/// Restrict the reference tree to the leaves with at least one read and the vertices on
/// their root paths (writeup §1.3; unary vertices are kept). `cols` are the likelihood
/// columns (sequences) that kept an alignment, in column order.
pub fn contract_to_observed(ref_tree: &RefTree, cols: &[SeqId]) -> Result<(TaxTree, Partition)> {
    if cols.is_empty() {
        bail!("no read aligned to any reference");
    }
    let leaf_of_col: Vec<NodeID> = cols
        .iter()
        .map(|s| ref_tree.leaf_of_seq[s.index()])
        .collect();
    let mut keep: HashSet<NodeID> = HashSet::new();
    for &leaf in &leaf_of_col {
        if keep.contains(&leaf) {
            continue;
        }
        let path = ref_tree
            .tree
            .node_to_root_ids(leaf)
            .map_err(|e| anyhow::anyhow!("taxonomy contraction failed: {}", e))?;
        for id in path {
            if !keep.insert(id) {
                break;
            }
        }
    }
    let mut tree = ref_tree.tree.clone();
    for v in postorder(&tree, ROOT) {
        if !keep.contains(&v) {
            tree.remove_node(v);
        }
    }
    set_weight(&mut tree, ROOT, 1.0);
    Ok((tree, Partition::new(leaf_of_col)))
}

/// Set every conditional from counts: `q_{c|v} = N_c / N_v` (0 when `N_v = 0`).
fn conditionals_from_counts(tree: &mut TaxTree) {
    let root = tree.get_root_id();
    set_weight(tree, root, 1.0);
    for v in preorder(tree) {
        let ch = children(tree, v);
        if ch.is_empty() {
            continue;
        }
        let total: EMProb = ch.iter().map(|&c| count(tree, c)).sum();
        for c in ch {
            let q = if total > 0.0 {
                count(tree, c) / total
            } else {
                0.0
            };
            set_weight(tree, c, q);
        }
    }
}

/// Initial conditionals from per-column initial proportions (treated as counts).
pub fn init_conditionals(tree: &mut TaxTree, part: &Partition, pi_contig: &[EMProb]) {
    subtree_counts(tree, part, pi_contig);
    conditionals_from_counts(tree);
}

/// Uniform initial conditionals: `q_{c|v} = 1 / |ch(v)|` at every internal vertex.
pub fn init_uniform(tree: &mut TaxTree) {
    let root = tree.get_root_id();
    set_weight(tree, root, 1.0);
    for v in preorder(tree) {
        let ch = children(tree, v);
        let q = 1.0 / ch.len() as EMProb;
        for c in ch {
            set_weight(tree, c, q);
        }
    }
}

/// `π_v` for every node id (product of conditionals along the root path; pruned ids 0).
pub fn leaf_marginals(tree: &TaxTree) -> Vec<EMProb> {
    let mut pi = vec![0.0; arena_len(tree)];
    let root = tree.get_root_id();
    for v in preorder(tree) {
        pi[v] = match tree.get_node_parent_id(v) {
            Some(p) if v != root => pi[p] * weight(tree, v),
            _ => 1.0,
        };
    }
    pi
}

/// Per-column mixing proportion `π̃_s = π_ℓ(s) / |P_ℓ(s)|` from per-node marginals
/// (0 for a pruned column).
pub fn contig_pi(part: &Partition, pi_node: &[EMProb]) -> Vec<EMProb> {
    (0..part.n_cols())
        .map(|col| {
            let leaf = part.leaf_of(col);
            if leaf == NO_NODE {
                0.0
            } else {
                pi_node[leaf] / part.size_of(leaf) as EMProb
            }
        })
        .collect()
}

/// Store `N_v` for every node from per-column expected counts.
pub fn subtree_counts(tree: &mut TaxTree, part: &Partition, e_contig: &[EMProb]) {
    let mut n = vec![0.0; arena_len(tree)];
    for (col, &e) in e_contig.iter().enumerate() {
        let leaf = part.leaf_of(col);
        if leaf != NO_NODE {
            n[leaf] += e;
        }
    }
    let root = tree.get_root_id();
    for v in postorder(tree, root) {
        if v != root {
            if let Some(p) = tree.get_node_parent_id(v) {
                n[p] += n[v];
            }
        }
    }
    for v in tree.get_node_ids().collect::<Vec<_>>() {
        set_count(tree, v, n[v]);
    }
}

/// What one pruning pass did.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct PruneReport {
    /// Vertices deleted below a kept vertex (with their subtrees).
    pub pruned: usize,
    /// Vertices that became leaves over the sequences beneath them.
    pub collapsed: usize,
    /// The root itself collapsed: every child fell below `ρ`.
    pub root_collapsed: bool,
}

impl PruneReport {
    pub fn changed(&self) -> bool {
        self.pruned + self.collapsed > 0
    }
}

/// The local test of the writeup's pruning step, decided for every vertex at once (counts
/// never grow down the tree, so this equals the top-down traversal): a vertex with
/// `N_v < ρ` is deleted; a vertex with `N_v ≥ ρ` whose children all have `N_c < ρ` is
/// collapsed into a leaf holding every sequence beneath it. The root is never deleted.
pub fn prune_and_collapse(tree: &mut TaxTree, part: &mut Partition, rho: EMProb) -> PruneReport {
    let root = tree.get_root_id();
    let mut deleted: HashSet<NodeID> = HashSet::new();
    let mut collapse: HashSet<NodeID> = HashSet::new();
    for v in preorder(tree) {
        if v != root && count(tree, v) < rho {
            deleted.insert(v);
            continue;
        }
        let ch = children(tree, v);
        if !ch.is_empty() && ch.iter().all(|&c| count(tree, c) < rho) {
            collapse.insert(v);
        }
    }
    if deleted.is_empty() {
        return PruneReport::default();
    }
    // Nearest surviving ancestor of a deleted vertex.
    let survivor = |tree: &TaxTree, v: NodeID| -> Option<NodeID> {
        let mut a = tree.get_node_parent_id(v);
        while let Some(p) = a {
            if !deleted.contains(&p) {
                return Some(p);
            }
            a = tree.get_node_parent_id(p);
        }
        None
    };
    let pruned = deleted
        .iter()
        .filter(|&&v| survivor(tree, v).is_none_or(|p| !collapse.contains(&p)))
        .count();
    for col in 0..part.leaf_of_col.len() {
        let leaf = part.leaf_of_col[col];
        if leaf == NO_NODE || !deleted.contains(&leaf) {
            continue;
        }
        part.leaf_of_col[col] = match survivor(tree, leaf) {
            Some(p) if collapse.contains(&p) => p,
            _ => NO_NODE,
        };
    }
    part.recount();
    for v in postorder(tree, root) {
        if deleted.contains(&v) {
            tree.remove_node(v);
        }
    }
    PruneReport {
        pruned,
        collapsed: collapse.len(),
        root_collapsed: collapse.contains(&root),
    }
}

/// Tangent-minorizer update of every vertex's conditionals from the stored counts;
/// vertices are independent and updated in parallel.
pub fn m_step(tree: &mut TaxTree, rho: EMProb, omega: EMProb) {
    let groups: Vec<(Vec<NodeID>, Vec<EMProb>, Vec<EMProb>)> = preorder(tree)
        .into_iter()
        .filter_map(|v| {
            let ch = children(tree, v);
            if ch.is_empty() {
                return None;
            }
            let n = ch.iter().map(|&c| count(tree, c)).collect();
            let q = ch.iter().map(|&c| weight(tree, c)).collect();
            Some((ch, n, q))
        })
        .collect();
    let updates: Vec<Vec<EMProb>> = groups
        .par_iter()
        .map(|(_, n, q)| vertex_update(n, q, rho, omega))
        .collect();
    for ((ch, _, _), q_new) in groups.iter().zip(updates) {
        for (&c, q) in ch.iter().zip(q_new) {
            set_weight(tree, c, q);
        }
    }
}

/// New conditionals of one vertex from its children's counts `n` and current `q`:
/// `q_c = N_c / (λ + w_c)` with `w_c = ρ / (ω + q_c)`, renormalised over the children with
/// reads; a single such child takes 1; `ρ = 0` is the count ratio.
fn vertex_update(n: &[EMProb], q: &[EMProb], rho: EMProb, omega: EMProb) -> Vec<EMProb> {
    let alive: Vec<usize> = (0..n.len()).filter(|&i| n[i] > 0.0).collect();
    let mut out = vec![0.0; n.len()];
    match alive.len() {
        0 => {}
        1 => out[alive[0]] = 1.0,
        _ => {
            let na: Vec<EMProb> = alive.iter().map(|&i| n[i]).collect();
            let vals: Vec<EMProb> = if rho == 0.0 {
                let total: EMProb = na.iter().sum();
                na.iter().map(|&x| x / total).collect()
            } else {
                let w: Vec<EMProb> = alive.iter().map(|&i| rho / (omega + q[i])).collect();
                let lambda = solve_lambda(&na, &w);
                let raw: Vec<EMProb> = na
                    .iter()
                    .zip(&w)
                    .map(|(&x, &wc)| x / (lambda + wc))
                    .collect();
                let s: EMProb = raw.iter().sum();
                raw.iter().map(|&x| x / s).collect()
            };
            for (&i, v) in alive.iter().zip(vals) {
                out[i] = v;
            }
        }
    }
    out
}

/// `ρ Σ log(1 + q_{c|v}/ω)` over every conditional in the tree.
pub fn penalty(tree: &TaxTree, rho: EMProb, omega: EMProb) -> EMProb {
    if rho == 0.0 {
        return 0.0;
    }
    let root = tree.get_root_id();
    rho * tree
        .get_node_ids()
        .filter(|&v| v != root)
        .map(|v| (weight(tree, v) / omega).ln_1p())
        .sum::<EMProb>()
}

/// The multiplier `λ` with `Σ_c n_c / (λ + w_c) = 1`, unique on `(−min w, Σ n)`, by
/// Newton–Raphson from `λ_0 = n_k − w_k`, `k = argmin w`: there `g(λ_0) ≥ 1`, and `g` is
/// convex and decreasing, so the iterates climb monotonically to the root without
/// overshooting. Every `n_c` must be positive.
pub fn solve_lambda(n: &[EMProb], w: &[EMProb]) -> EMProb {
    debug_assert_eq!(n.len(), w.len());
    let total: EMProb = n.iter().sum();
    if n.len() == 1 {
        return total - w[0];
    }
    if w.iter().all(|&x| x == 0.0) {
        return total;
    }
    let g = |l: EMProb| {
        n.iter()
            .zip(w)
            .map(|(&ni, &wi)| ni / (l + wi))
            .sum::<EMProb>()
            - 1.0
    };
    let dg = |l: EMProb| {
        -n.iter()
            .zip(w)
            .map(|(&ni, &wi)| ni / ((l + wi) * (l + wi)))
            .sum::<EMProb>()
    };
    let (k, &wmin) = w
        .iter()
        .enumerate()
        .min_by(|a, b| a.1.total_cmp(b.1))
        .expect("at least two children");
    let mut lambda = n[k] - wmin;
    for _ in 0..100 {
        let gv = g(lambda);
        if gv.abs() < 1e-12 {
            break;
        }
        let next = lambda - gv / dg(lambda);
        debug_assert!(next.is_finite() && next > -wmin, "Newton left the bracket");
        lambda = next;
    }
    lambda
}

/// Descend from the root to the child with the largest subtree posterior while it
/// exceeds `gamma`. `hits` are the read's posteriors per leaf. Returns the vertex and its
/// subtree posterior.
pub fn assign_read(tree: &TaxTree, hits: &[(NodeID, EMProb)], gamma: EMProb) -> (NodeID, EMProb) {
    let mut h: HashMap<NodeID, EMProb> = HashMap::new();
    for &(leaf, p) in hits {
        if let Ok(path) = tree.node_to_root_ids(leaf) {
            for id in path {
                *h.entry(id).or_insert(0.0) += p;
            }
        }
    }
    let mut v = tree.get_root_id();
    let mut hv = h.get(&v).copied().unwrap_or(0.0);
    loop {
        let best = children(tree, v)
            .into_iter()
            .map(|c| (c, h.get(&c).copied().unwrap_or(0.0)))
            .max_by(|a, b| a.1.total_cmp(&b.1));
        match best {
            Some((c, hc)) if hc > 0.0 && hc > gamma => {
                v = c;
                hv = hc;
            }
            _ => break,
        }
    }
    (v, hv)
}

/// Pre-order intervals for O(1) subtree membership tests.
pub struct SubtreeIndex {
    /// Entry and exit clock of every vertex by id; `usize::MAX` for an id not in the tree.
    tin: Vec<usize>,
    tout: Vec<usize>,
}

impl SubtreeIndex {
    pub fn new(tree: &TaxTree) -> Self {
        let n = arena_len(tree);
        let mut tin = vec![usize::MAX; n];
        let mut tout = vec![usize::MAX; n];
        let mut clock = 0;
        let mut stack: Vec<(NodeID, bool)> = vec![(tree.get_root_id(), false)];
        while let Some((v, done)) = stack.pop() {
            if done {
                tout[v] = clock;
                continue;
            }
            tin[v] = clock;
            clock += 1;
            stack.push((v, true));
            for c in children(tree, v) {
                stack.push((c, false));
            }
        }
        Self { tin, tout }
    }

    /// True when `node` is `ancestor` or lies in its subtree.
    pub fn contains(&self, ancestor: NodeID, node: NodeID) -> bool {
        let at = |v: &[usize], id: NodeID| v.get(id).copied().unwrap_or(usize::MAX);
        let (a_in, a_out, n_in) = (
            at(&self.tin, ancestor),
            at(&self.tout, ancestor),
            at(&self.tin, node),
        );
        a_in != usize::MAX && n_in != usize::MAX && a_in <= n_in && n_in < a_out
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::{build_index_from_bytes, load_index};
    use crate::RefIndex;

    const FASTA: &[u8] = b">s0\nACGTACGTACGTACGTACGTACGTACGT\n>s1\nTTAGCCTTAGCCTTAGCCTTAGCCTTAG\n>s2\nGGGACCCTTTGGGACCCTTTGGGACCCT\n>s3\nCATCATCATCATCATCATCATCATCATC\n";

    fn index() -> RefIndex {
        let (bytes, _) = build_index_from_bytes(FASTA).expect("build");
        load_index(&bytes, "test.fmidx").expect("load")
    }

    fn grouping(pairs: &[(&str, u32)]) -> Grouping {
        let map: HashMap<String, u32> = pairs.iter().map(|(a, b)| (a.to_string(), *b)).collect();
        Grouping::from_pairs(&index(), &map)
    }

    const NODES: &str = "1\t|\t1\t|\tno rank\t|\n10\t|\t1\t|\tgenus\t|\n11\t|\t10\t|\tspecies\t|\n12\t|\t10\t|\tspecies\t|\n20\t|\t1\t|\tgenus\t|\n21\t|\t20\t|\tspecies\t|\n";

    fn nodes() -> TaxNodes {
        TaxNodes::parse(NODES).unwrap()
    }

    fn seqs(ids: &[u32]) -> Vec<SeqId> {
        ids.iter().map(|&i| SeqId::new(i)).collect()
    }

    fn by_label(tree: &TaxTree, l: &str) -> NodeID {
        tree.get_node_ids()
            .find(|&id| tree.get_node(id).unwrap().get_taxa().map(String::as_str) == Some(l))
            .unwrap_or_else(|| panic!("no node labelled {l}"))
    }

    /// Lineage labels from a node up to the root.
    fn lineage(tree: &TaxTree, id: NodeID) -> Vec<String> {
        tree.node_to_root_ids(id)
            .unwrap()
            .map(|v| label(tree, v).to_string())
            .collect()
    }

    // --- TaxNodes ---------------------------------------------------------------------

    #[test]
    fn parse_pipes_and_plain_tsv() {
        let t = nodes();
        assert_eq!(t.len(), 6);
        assert_eq!(t.parent(11), Some(10));
        assert_eq!(t.parent(1), None);
        assert_eq!(t.rank(10), "genus");
        assert_eq!(t.rank(99), "-");
        let p = TaxNodes::parse("# c\n1 1\n2 1 kingdom\n\n3\t2\n").unwrap();
        assert_eq!(p.parent(3), Some(2));
        assert_eq!(p.rank(2), "kingdom");
        assert_eq!(p.rank(3), "-");
    }

    #[test]
    fn parse_rejects_garbage() {
        assert!(TaxNodes::parse("1\n")
            .unwrap_err()
            .to_string()
            .contains("line 1"));
        assert!(TaxNodes::parse("a 1\n").is_err());
        assert!(TaxNodes::parse("2 1\n2 3\n")
            .unwrap_err()
            .to_string()
            .contains("parents"));
        assert!(TaxNodes::parse("0 1\n").is_err());
    }

    // --- construction -------------------------------------------------------------------

    #[test]
    fn reference_tree_puts_every_sequence_under_its_taxid() {
        // s0, s1 -> 11; s2 -> 12 (siblings under genus 10); s3 unmapped
        let g = grouping(&[("s0", 11), ("s1", 11), ("s2", 12)]);
        let rt = reference_tree(&g, Some(&nodes())).unwrap();
        let t = &rt.tree;
        assert_eq!(rt.n_unknown, 0);
        assert_eq!(rt.n_unmapped, 1);
        assert_eq!(rt.n_tips(), 4);
        assert_eq!(t.num_nodes(), 5 + 4, "root, 1, 10, 11, 12 + 4 sequences");
        let sp11 = by_label(t, "11");
        assert_eq!(children(t, sp11).len(), 2, "s0 and s1 are leaves of 11");
        assert_eq!(
            lineage(t, rt.leaf_of_seq[0]),
            ["s0", "11", "10", "1", "root"]
        );
        assert_eq!(
            lineage(t, rt.leaf_of_seq[2]),
            ["s2", "12", "10", "1", "root"]
        );
        assert_eq!(lineage(t, rt.leaf_of_seq[3]), ["s3", "root"]);
        for &l in &rt.leaf_of_seq {
            assert!(t.is_leaf(l));
            assert!(is_sequence(t, l));
        }
        assert!(!is_sequence(t, sp11));
        assert!(
            t.get_leaf_ids().all(|l| rt.leaf_of_seq.contains(&l)),
            "every tip is a sequence"
        );
    }

    #[test]
    fn without_nodes_taxids_hang_under_root() {
        let g = grouping(&[("s0", 11), ("s1", 11), ("s2", 21)]);
        let rt = reference_tree(&g, None).unwrap();
        let t = &rt.tree;
        assert_eq!(rt.n_unknown, 0, "no nodes file: nothing is unknown");
        assert_eq!(lineage(t, rt.leaf_of_seq[0]), ["s0", "11", "root"]);
        assert_eq!(lineage(t, rt.leaf_of_seq[1]), ["s1", "11", "root"]);
        assert_eq!(lineage(t, rt.leaf_of_seq[2]), ["s2", "21", "root"]);
        assert_eq!(lineage(t, rt.leaf_of_seq[3]), ["s3", "root"]);
        // identity grouping: the star over the sequences
        let star = reference_tree(&Grouping::identity(&index()), None).unwrap();
        assert_eq!(star.tree.num_nodes(), 5);
        assert!(star
            .leaf_of_seq
            .iter()
            .all(|&l| star.tree.get_node_parent_id(l) == Some(ROOT)));
    }

    #[test]
    fn unknown_taxid_hangs_under_root() {
        let g = grouping(&[("s0", 99)]);
        let rt = reference_tree(&g, Some(&nodes())).unwrap();
        assert_eq!(rt.n_unknown, 1);
        let n = by_label(&rt.tree, "99");
        assert_eq!(rt.tree.get_node_parent_id(n), Some(ROOT));
        assert_eq!(rt.tree.get_node_parent_id(rt.leaf_of_seq[0]), Some(n));
    }

    #[test]
    fn contraction_keeps_root_paths_and_unary_vertices() {
        let g = grouping(&[("s0", 11), ("s1", 12), ("s2", 21)]);
        let rt = reference_tree(&g, Some(&nodes())).unwrap();
        // only s0 has reads: root -> 1 -> 10 -> 11 -> s0 survive, nothing else
        let (t, part) = contract_to_observed(&rt, &seqs(&[0])).unwrap();
        assert_eq!(t.get_root_id(), ROOT);
        assert_eq!(t.num_nodes(), 5);
        assert_eq!(
            lineage(&t, part.leaf_of(0)),
            ["s0", "11", "10", "1", "root"]
        );
        assert_eq!(part.size_of(part.leaf_of(0)), 1);
        assert_eq!(part.leaf_of(0), rt.leaf_of_seq[0], "ids preserved");
        // s0 and s2: genus 10 and 20 both kept, 12 and s1/s3 gone
        let (t, part) = contract_to_observed(&rt, &seqs(&[0, 2])).unwrap();
        assert_eq!(t.num_nodes(), 8);
        assert_eq!(part.n_cols(), 2);
        assert_eq!(
            lineage(&t, part.leaf_of(1)),
            ["s2", "21", "20", "1", "root"]
        );
        assert!(t.get_node(rt.leaf_of_seq[1]).is_none());
        assert!(contract_to_observed(&rt, &[]).is_err());
    }

    // --- EM pieces ----------------------------------------------------------------------

    #[test]
    fn marginals_multiply_along_the_path_and_contig_pi_divides_by_size() {
        let g = grouping(&[("s0", 11), ("s1", 11), ("s2", 12)]);
        let rt = reference_tree(&g, Some(&nodes())).unwrap();
        let (mut t, part) = contract_to_observed(&rt, &seqs(&[0, 1, 2])).unwrap();
        let (sp11, sp12) = (by_label(&t, "11"), by_label(&t, "12"));
        set_weight(&mut t, sp11, 0.25);
        set_weight(&mut t, sp12, 0.75);
        set_weight(&mut t, part.leaf_of(0), 0.5);
        set_weight(&mut t, part.leaf_of(1), 0.5);
        set_weight(&mut t, part.leaf_of(2), 1.0);
        for v in [by_label(&t, "1"), by_label(&t, "10")] {
            set_weight(&mut t, v, 1.0);
        }
        let pi = leaf_marginals(&t);
        assert_eq!(pi[t.get_root_id()], 1.0);
        assert_eq!(pi[sp11], 0.25);
        assert_eq!(pi[part.leaf_of(0)], 0.125);
        let cp = contig_pi(&part, &pi);
        assert_eq!(cp, vec![0.125, 0.125, 0.75]);
    }

    #[test]
    fn uniform_init_splits_every_vertex_evenly() {
        let g = grouping(&[("s0", 11), ("s1", 11), ("s2", 12)]);
        let rt = reference_tree(&g, Some(&nodes())).unwrap();
        let (mut t, part) = contract_to_observed(&rt, &seqs(&[0, 1, 2])).unwrap();
        init_uniform(&mut t);
        let pi = leaf_marginals(&t);
        assert_eq!(pi[by_label(&t, "10")], 1.0, "unary chain above the genus");
        assert_eq!(pi[by_label(&t, "11")], 0.5);
        assert_eq!(pi[part.leaf_of(0)], 0.25);
        assert_eq!(pi[part.leaf_of(1)], 0.25);
        assert_eq!(pi[part.leaf_of(2)], 0.5);
    }

    #[test]
    fn m_step_without_penalty_is_the_count_ratio() {
        let g = grouping(&[("s0", 11), ("s1", 12), ("s2", 21)]);
        let rt = reference_tree(&g, Some(&nodes())).unwrap();
        let (mut t, part) = contract_to_observed(&rt, &seqs(&[0, 1, 2])).unwrap();
        subtree_counts(&mut t, &part, &[30.0, 10.0, 60.0]);
        m_step(&mut t, 0.0, 1e-10);
        let pi = leaf_marginals(&t);
        assert!((pi[part.leaf_of(0)] - 0.3).abs() < 1e-12);
        assert!((pi[part.leaf_of(1)] - 0.1).abs() < 1e-12);
        assert!((pi[part.leaf_of(2)] - 0.6).abs() < 1e-12);
        assert!(
            (weight(&t, by_label(&t, "11")) - 0.75).abs() < 1e-12,
            "q_{{11|10}}"
        );
        assert_eq!(weight(&t, part.leaf_of(0)), 1.0, "unary edge");
    }

    #[test]
    fn penalised_m_step_shrinks_a_weak_child_geometrically() {
        let g = grouping(&[("s0", 11), ("s1", 12)]);
        let rt = reference_tree(&g, Some(&nodes())).unwrap();
        let (mut t, part) = contract_to_observed(&rt, &seqs(&[0, 1])).unwrap();
        subtree_counts(&mut t, &part, &[1000.0, 20.0]);
        conditionals_from_counts(&mut t);
        let weak = by_label(&t, "12");
        let mut prev = weight(&t, weak);
        for _ in 0..10 {
            m_step(&mut t, 150.0, 1e-10);
            let q = weight(&t, weak);
            assert!(q < prev, "{q} !< {prev}");
            prev = q;
        }
        assert!(prev < 1e-6);
        let s: EMProb = children(&t, by_label(&t, "10"))
            .iter()
            .map(|&c| weight(&t, c))
            .sum();
        assert!((s - 1.0).abs() < 1e-12, "renormalised: {s}");
    }

    #[test]
    fn m_step_renormalises_after_the_penalised_update() {
        let q = vertex_update(&[300.0, 5.0, 40.0], &[0.8, 0.02, 0.18], 150.0, 1e-10);
        assert!((q.iter().sum::<EMProb>() - 1.0).abs() < 1e-12, "{q:?}");
        // a child without reads keeps 0, a lone child takes 1
        assert_eq!(
            vertex_update(&[0.0, 7.0], &[0.5, 0.5], 150.0, 1e-10),
            vec![0.0, 1.0]
        );
        assert_eq!(
            vertex_update(&[0.0, 0.0], &[0.5, 0.5], 1.0, 1e-10),
            vec![0.0, 0.0]
        );
    }

    #[test]
    fn prune_deletes_subtree_and_collapse_merges_contigs() {
        // genus 10 {11: s0, 12: s1}, genus 20 {21: s2}
        let g = grouping(&[("s0", 11), ("s1", 12), ("s2", 21)]);
        let rt = reference_tree(&g, Some(&nodes())).unwrap();
        let (mut t, mut part) = contract_to_observed(&rt, &seqs(&[0, 1, 2])).unwrap();
        // 12 is weak: 12 and s1 pruned; 10 keeps 11
        subtree_counts(&mut t, &part, &[300.0, 5.0, 300.0]);
        let r = prune_and_collapse(&mut t, &mut part, 50.0);
        assert_eq!((r.pruned, r.collapsed, r.root_collapsed), (2, 0, false));
        assert_eq!(part.leaf_of(1), NO_NODE);
        assert_eq!(t.num_nodes(), 8, "root, 1, 10, 11, s0, 20, 21, s2");
        assert!(t.get_node(by_label(&rt.tree, "12")).is_none());
        // every child of the root weak: the whole sample collapses into the root
        subtree_counts(&mut t, &part, &[20.0, 0.0, 20.0]);
        let r = prune_and_collapse(&mut t, &mut part, 50.0);
        assert_eq!((r.pruned, r.collapsed, r.root_collapsed), (0, 1, true));
        let root = t.get_root_id();
        assert!(t.is_leaf(root));
        assert_eq!(t.num_nodes(), 1);
        assert_eq!(part.leaf_of(0), root);
        assert_eq!(part.leaf_of(2), root);
        assert_eq!(part.size_of(root), 2);
        assert!(!prune_and_collapse(&mut t, &mut part, 50.0).changed());
    }

    #[test]
    fn collapse_merges_children_into_an_internal_vertex() {
        let g = grouping(&[("s0", 11), ("s1", 11), ("s2", 21)]);
        let rt = reference_tree(&g, Some(&nodes())).unwrap();
        let (mut t, mut part) = contract_to_observed(&rt, &seqs(&[0, 1, 2])).unwrap();
        // species 11 has 60 reads split 20/40 over its sequences, below rho = 50
        subtree_counts(&mut t, &part, &[20.0, 40.0, 300.0]);
        let r = prune_and_collapse(&mut t, &mut part, 50.0);
        assert_eq!((r.pruned, r.collapsed, r.root_collapsed), (0, 1, false));
        let sp = by_label(&t, "11");
        assert!(t.is_leaf(sp));
        assert!(!is_sequence(&t, sp));
        assert_eq!(part.leaf_of(0), sp);
        assert_eq!(part.leaf_of(1), sp);
        assert_eq!(part.size_of(sp), 2);
        assert_eq!(t.num_nodes(), 7, "root, 1, 10, 11, 20, 21, s2");
    }

    // --- solve_lambda ---------------------------------------------------------------------

    fn q_of(n: &[EMProb], w: &[EMProb]) -> Vec<EMProb> {
        let l = solve_lambda(n, w);
        n.iter().zip(w).map(|(&ni, &wi)| ni / (l + wi)).collect()
    }

    #[test]
    fn lambda_without_weights_is_the_total() {
        let q = q_of(&[30.0, 10.0, 60.0], &[0.0, 0.0, 0.0]);
        assert_eq!(q, vec![0.3, 0.1, 0.6]);
    }

    #[test]
    fn lambda_sums_to_one_and_huge_weight_kills_a_child() {
        let q = q_of(&[300.0, 5.0], &[150.0 / 0.98, 150.0 / 0.02]);
        assert!((q.iter().sum::<EMProb>() - 1.0).abs() < 1e-9, "{q:?}");
        assert!(q[1] < 0.02);
        let q = q_of(&[1.0, 1.0, 1.0], &[0.0, 1e6, 1e9]);
        assert!((q.iter().sum::<EMProb>() - 1.0).abs() < 1e-9, "{q:?}");
        assert!(q[0] > 0.999);
        // a wide sibling group with very unequal weights
        let n: Vec<EMProb> = (1..=200).map(|i| i as EMProb).collect();
        let w: Vec<EMProb> = (1..=200)
            .map(|i| 150.0 / (1e-10 + 1.0 / i as EMProb))
            .collect();
        let q = q_of(&n, &w);
        assert!((q.iter().sum::<EMProb>() - 1.0).abs() < 1e-9);
        assert!(q.iter().all(|&x| x > 0.0 && x.is_finite()));
    }

    #[test]
    fn lambda_can_be_negative() {
        // large weights relative to counts: λ ∈ (−min w, 0)
        let n = [2.0, 3.0];
        let w = [10.0, 10.0];
        let l = solve_lambda(&n, &w);
        assert!(l < 0.0 && l > -10.0, "{l}");
        let q = q_of(&n, &w);
        assert!((q.iter().sum::<EMProb>() - 1.0).abs() < 1e-9);
    }

    // --- assignment -----------------------------------------------------------------------

    #[test]
    fn assignment_stops_where_the_posterior_no_longer_exceeds_gamma() {
        let g = grouping(&[("s0", 11), ("s1", 12), ("s2", 10)]);
        let rt = reference_tree(&g, Some(&nodes())).unwrap();
        let (t, part) = contract_to_observed(&rt, &seqs(&[0, 1, 2])).unwrap();
        let (l0, l1, l2) = (part.leaf_of(0), part.leaf_of(1), part.leaf_of(2));
        let genus = by_label(&t, "10");
        // split evenly over the species: stays at the genus for gamma 0.6
        let (v, h) = assign_read(&t, &[(l0, 0.5), (l1, 0.5)], 0.6);
        assert_eq!(v, genus);
        assert!((h - 1.0).abs() < 1e-12);
        // confident read descends through the unary species vertex to its sequence
        let (v, h) = assign_read(&t, &[(l0, 0.9), (l1, 0.1)], 0.6);
        assert_eq!((v, h), (l0, 0.9));
        // a sequence that sits directly under the genus is a sibling of the species
        let (v, h) = assign_read(&t, &[(l2, 0.8), (l1, 0.2)], 0.6);
        assert_eq!((v, h), (l2, 0.8));
        // a child must exceed gamma, not merely reach it
        let (v, h) = assign_read(&t, &[(l0, 0.6), (l1, 0.4)], 0.6);
        assert_eq!((v, h), (genus, 1.0));
        // the unary vertices root -> 1 -> 10 carry the whole mass; the split below the
        // genus is too even for gamma 0.7
        let (v, h) = assign_read(&t, &[(l0, 0.4), (l2, 0.6)], 0.7);
        assert_eq!((v, h), (genus, 1.0));
        // no hit at all: the synthetic root
        let (v, _) = assign_read(&t, &[], 0.6);
        assert!(is_synthetic_root(&t, v));
        let idx = SubtreeIndex::new(&t);
        assert!(idx.contains(genus, l0));
        assert!(idx.contains(l0, l0));
        assert!(!idx.contains(l0, l1));
    }
}
