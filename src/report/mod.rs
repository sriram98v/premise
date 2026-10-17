//! The run report: a NEXUS file holding the fitted taxonomy (down to the index
//! sequences) with the fitted proportions, the sequence → vertex map, every alignment
//! log-likelihood of every read, the EM trace and the run metadata. Only what the EM
//! cannot regenerate is stored: the posteriors follow from the likelihoods and the
//! proportions ([`posterior`]), so `.matches` for any `gamma` can be regenerated without
//! the index or the reads ([`tables`]).

pub mod export;
pub mod nexus;
pub mod posterior;
pub mod read;
pub mod tables;
pub mod write;

use crate::em::EmFit;
use crate::tree::{
    arena_len, children, conditional, count, label, new_node, preorder, rank, Partition, TaxTree,
    NO_NODE, ROOT,
};
use crate::utils::EMProb;
use crate::SeqId;
use anyhow::{bail, Result};
use phylo::node::NodeID;
use phylo::prelude::*;
use std::collections::{HashMap, HashSet};

/// Format version written in the PREMISE block.
pub const FORMAT_VERSION: u32 = 4;

/// What the EM did to a vertex of the contracted tree.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum VertexStatus {
    /// In the final tree (an internal vertex, or a sequence tip).
    Active,
    /// A final-tree leaf that holds the sequences beneath it.
    Collapsed,
    /// Beneath a collapsed vertex.
    Pooled,
    /// Deleted by the penalty, with no collapsed ancestor.
    Pruned,
}

impl VertexStatus {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Active => "active",
            Self::Collapsed => "collapsed",
            Self::Pooled => "pooled",
            Self::Pruned => "pruned",
        }
    }

    pub fn parse(s: &str) -> Option<Self> {
        match s {
            "active" => Some(Self::Active),
            "collapsed" => Some(Self::Collapsed),
            "pooled" => Some(Self::Pooled),
            "pruned" => Some(Self::Pruned),
            _ => None,
        }
    }

    /// In the final EM tree.
    pub fn is_final(self) -> bool {
        matches!(self, Self::Active | Self::Collapsed)
    }
}

/// One vertex of the stored tree.
#[derive(Clone, Debug, PartialEq)]
pub struct ReportVertex {
    pub id: NodeID,
    pub parent: Option<NodeID>,
    pub label: String,
    pub rank: String,
    pub status: VertexStatus,
    /// Sequences the vertex holds as an EM leaf (0 for internal, pooled, pruned).
    pub contigs: usize,
    /// `N_v` (final-tree vertices only).
    pub count: Option<EMProb>,
    /// `π_v` for final-tree vertices; `π̃_s` for a pooled sequence; 0 otherwise.
    pub proportion: EMProb,
    /// `q_{v|parent}` (final-tree vertices only).
    pub conditional: Option<EMProb>,
}

/// The contracted tree (every sequence with reads at a tip) with the EM's verdict on
/// each vertex, in pre-order.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ReportTree {
    pub vertices: Vec<ReportVertex>,
}

impl ReportTree {
    /// Overlay the fitted tree on `stored`, the clone taken before the EM; `tips_of_col`
    /// is the tip of every likelihood column in that clone.
    pub fn from_fit(
        stored: &TaxTree,
        tips_of_col: &[NodeID],
        fitted: &TaxTree,
        part: &Partition,
        fit: &EmFit,
    ) -> Self {
        let alive: HashSet<NodeID> = fitted.get_node_ids().collect();
        let col_of_tip: HashMap<NodeID, usize> = tips_of_col
            .iter()
            .enumerate()
            .map(|(c, &t)| (t, c))
            .collect();
        let mut status_of: HashMap<NodeID, VertexStatus> = HashMap::new();
        let mut vertices = Vec::with_capacity(arena_len(stored));
        for v in preorder(stored) {
            let parent = stored.get_node_parent_id(v);
            let status = if alive.contains(&v) {
                if children(fitted, v).is_empty() && !children(stored, v).is_empty() {
                    VertexStatus::Collapsed
                } else {
                    VertexStatus::Active
                }
            } else {
                let under_collapsed = parent.is_some_and(|p| {
                    matches!(
                        status_of.get(&p),
                        Some(VertexStatus::Collapsed) | Some(VertexStatus::Pooled)
                    )
                });
                if under_collapsed {
                    VertexStatus::Pooled
                } else {
                    VertexStatus::Pruned
                }
            };
            status_of.insert(v, status);
            let (contigs, count_v, proportion, cond) = match status {
                VertexStatus::Active | VertexStatus::Collapsed => (
                    part.size_of(v),
                    Some(count(fitted, v)),
                    fit.node_pi.get(v).copied().unwrap_or(0.0),
                    Some(conditional(fitted, v)),
                ),
                VertexStatus::Pooled => {
                    let p = col_of_tip
                        .get(&v)
                        .and_then(|&c| fit.pi_tilde.get(c))
                        .copied()
                        .unwrap_or(0.0);
                    (0, None, p, None)
                }
                VertexStatus::Pruned => (0, None, 0.0, None),
            };
            vertices.push(ReportVertex {
                id: v,
                parent,
                label: label(stored, v).to_string(),
                rank: rank(stored, v).to_string(),
                status,
                contigs,
                count: count_v,
                proportion,
                conditional: cond,
            });
        }
        Self { vertices }
    }

    pub fn root(&self) -> NodeID {
        self.vertices.first().map_or(0, |v| v.id)
    }

    /// `(leaf, π_ℓ)` of the final tree's leaves with mass, sorted by label: the
    /// proportions the EM ends with. A final vertex is a leaf when no final vertex hangs
    /// below it, and its `proportion` is the EM's `π_ℓ`, so the report need not store
    /// these separately.
    pub fn final_leaf_props(&self) -> Vec<(NodeID, EMProb)> {
        let status: HashMap<NodeID, VertexStatus> =
            self.vertices.iter().map(|v| (v.id, v.status)).collect();
        let has_final_child: HashSet<NodeID> = self
            .vertices
            .iter()
            .filter(|v| v.status.is_final())
            .filter_map(|v| v.parent)
            .filter(|p| status.get(p).is_some_and(|s| s.is_final()))
            .collect();
        let mut leaves: Vec<&ReportVertex> = self
            .vertices
            .iter()
            .filter(|v| v.status.is_final() && !has_final_child.contains(&v.id))
            .filter(|v| v.proportion > 0.0)
            .collect();
        leaves.sort_by(|a, b| a.label.cmp(&b.label));
        leaves.into_iter().map(|v| (v.id, v.proportion)).collect()
    }

    /// Label of every vertex by id.
    pub fn labels(&self) -> HashMap<NodeID, &str> {
        self.vertices
            .iter()
            .map(|v| (v.id, v.label.as_str()))
            .collect()
    }

    /// The final EM tree: the `active` and `collapsed` vertices with their ids,
    /// children order, counts and conditionals, so that read assignment on it is the
    /// assignment the run made.
    pub fn to_tax_tree(&self) -> TaxTree {
        let root = self.root();
        let mut tree = TaxTree::new(root);
        for v in &self.vertices {
            if !v.status.is_final() {
                continue;
            }
            let mut node = new_node(v.id, &v.label, &v.rank);
            node.set_weight(v.conditional);
            node.set_zeta(v.count);
            match v.parent {
                None => tree.set_node(node),
                Some(p) => tree.add_child(p, node),
            }
        }
        tree
    }

    /// The contracted tree as it was before the EM: every vertex, fresh node defaults
    /// (no counts, no conditionals, root weight 1), so a refit starts where `query` did.
    /// The ids are kept; the root must be [`ROOT`], which the EM assumes.
    pub fn to_stored_tax_tree(&self) -> Result<TaxTree> {
        if self.vertices.is_empty() {
            bail!("the report holds no tree");
        }
        let root = self.root();
        if root != ROOT {
            bail!(
                "the stored tree's root is vertex {} (expected {})",
                root,
                ROOT
            );
        }
        let mut tree = TaxTree::new(root);
        for v in &self.vertices {
            let mut node = new_node(v.id, &v.label, &v.rank);
            match v.parent {
                None => {
                    node.set_weight(Some(1.0));
                    tree.set_node(node);
                }
                Some(p) => tree.add_child(p, node),
            }
        }
        Ok(tree)
    }

    /// The EM leaf that owns each vertex: the nearest collapsed ancestor-or-self of a
    /// pooled vertex, else the vertex itself ([`NO_NODE`] when pruned).
    pub fn em_leaf_of(&self) -> HashMap<NodeID, NodeID> {
        let mut owner: HashMap<NodeID, NodeID> = HashMap::with_capacity(self.vertices.len());
        for v in &self.vertices {
            let leaf = match v.status {
                VertexStatus::Pooled => v
                    .parent
                    .and_then(|p| owner.get(&p).copied())
                    .unwrap_or(NO_NODE),
                VertexStatus::Pruned => NO_NODE,
                _ => v.id,
            };
            owner.insert(v.id, leaf);
        }
        owner
    }
}

/// Which stored-tree tip each likelihood column (sequence) sits at; sorted by [`SeqId`].
#[derive(Clone, Debug, Default, PartialEq)]
pub struct RefMap {
    pub entries: Vec<(SeqId, NodeID)>,
}

impl RefMap {
    pub fn new(seqs: &[SeqId], tips: &[NodeID]) -> Self {
        let mut entries: Vec<(SeqId, NodeID)> =
            seqs.iter().copied().zip(tips.iter().copied()).collect();
        entries.sort_unstable_by_key(|e| e.0.index());
        Self { entries }
    }

    pub fn tip_of(&self, seq: SeqId) -> Option<NodeID> {
        self.entries
            .binary_search_by_key(&seq.index(), |e| e.0.index())
            .ok()
            .map(|i| self.entries[i].1)
    }
}

/// The parameters of a run.
#[derive(Clone, Debug, PartialEq)]
pub struct RunParams {
    pub mem: usize,
    pub eps_1: EMProb,
    pub eps_2: EMProb,
    pub iter: usize,
    pub rho: EMProb,
    pub omega: EMProb,
    pub em_threshold: EMProb,
    pub init: String,
    pub penalty: bool,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct TaxonomyStats {
    pub ref_vertices: usize,
    pub ref_tips: usize,
    pub unknown_taxids: usize,
    pub unmapped_seqs: usize,
    pub contracted_vertices: usize,
    pub contracted_tips: usize,
    pub final_tips: usize,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct EmStats {
    pub iterations: usize,
    pub pruned: usize,
    pub collapsed: usize,
    pub root_collapsed: bool,
    pub final_ll: EMProb,
}

/// Seconds spent before the report is written.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct PhaseTimes {
    pub load_s: f64,
    pub align_s: f64,
    pub em_s: f64,
}

/// The PREMISE block. Input files are named by their basename only.
#[derive(Clone, Debug, PartialEq)]
pub struct ReportMeta {
    pub version: u32,
    pub tool: String,
    pub timestamp: String,
    /// `RefID` or `TaxID`: the second column of `.matches`.
    pub label_header: String,
    pub index: String,
    /// SHA-256 of the index file, lowercase hex; `None` when it could not be computed.
    pub index_sha256: Option<String>,
    pub r1: String,
    pub r2: String,
    pub seq2taxid: Option<String>,
    pub nodes: Option<String>,
    /// The report an `analyze` refit started from (basename); `None` for a `query` run.
    pub source: Option<String>,
    pub params: RunParams,
    pub threads: usize,
    pub reads_total: usize,
    pub reads_aligned: usize,
    pub taxonomy: TaxonomyStats,
    pub em: EmStats,
    pub runtime: PhaseTimes,
}

/// One EM iteration of the EM_TRACE block.
#[derive(Clone, Debug, PartialEq)]
pub struct EmIter {
    pub ll: EMProb,
    /// `(leaf, π_ℓ)` by node id.
    pub leaves: Vec<(NodeID, EMProb)>,
}

impl ReportMeta {
    /// The EM that produced the proportions, on one line for the head of `.props` and
    /// `.matches`: its settings (as in PARAMS), the iterations it ran, its final penalized
    /// log-likelihood, and for a refit the report it started from. Holds no `[` or `]`, so
    /// it fits in a NEXUS comment.
    pub fn em_summary(&self) -> String {
        let (p, e) = (&self.params, &self.em);
        let final_ll = if e.final_ll.is_finite() {
            format!("{:e}", e.final_ll)
        } else {
            "-".to_string()
        };
        let source = self.source.as_deref().map_or_else(String::new, |s| {
            format!(", refit of {}", s.replace(['[', ']'], "_"))
        });
        format!(
            "{} EM: init={} iter={} rho={:e} omega={:e} em_threshold={:e} penalty={}; {} iterations, final penalized log-likelihood {}{}",
            self.tool.replace(['[', ']'], "_"),
            p.init,
            p.iter,
            p.rho,
            p.omega,
            p.em_threshold,
            u8::from(p.penalty),
            e.iterations,
            final_ll,
            source
        )
    }
}

/// Everything in the report before the per-read blocks.
#[derive(Clone, Debug, PartialEq)]
pub struct ReportHead {
    pub meta: ReportMeta,
    pub em_trace: Vec<EmIter>,
    pub tree: ReportTree,
    pub refmap: RefMap,
    /// `(node, π_ℓ)` of the final-tree leaves with mass, sorted by label: always
    /// [`ReportTree::final_leaf_props`] of `tree`, never written to the report.
    pub props: Vec<(NodeID, EMProb)>,
}

/// The PROFILE block, known only once everything is written.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Profile {
    pub write_s: f64,
    pub total_s: f64,
    pub peak_rss_kb: Option<u64>,
}

/// One alignment of a read: a row of ALIGNMENT_LIKELIHOODS.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AlignRow {
    pub seq: SeqId,
    /// `full_match_ll`, the alignment log-likelihood (finite).
    pub loglik: EMProb,
    pub pos: (usize, usize),
}

/// One read of ALIGNMENT_LIKELIHOODS, its rows in [`SeqId`] order.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ReadRows {
    pub read_id: String,
    pub rows: Vec<AlignRow>,
}

#[cfg(test)]
pub(crate) mod fixtures {
    //! A small stored tree and its fitted counterpart shared by the report tests.
    use super::*;
    use crate::tree::{new_tree, prune_and_collapse, subtree_counts};

    /// root 0 → genus 1 → species 3 → seqs 6, 7; species 4 → seq 8;
    ///        → genus 2 → species 5 → seq 9. Columns 0..4 sit at tips 6, 7, 8, 9.
    pub fn stored() -> (TaxTree, Vec<NodeID>) {
        let mut t = new_tree();
        t.add_child(0, new_node(1, "A", "genus"));
        t.add_child(0, new_node(2, "B", "genus"));
        t.add_child(1, new_node(3, "sp1", "species"));
        t.add_child(1, new_node(4, "sp2", "species"));
        t.add_child(2, new_node(5, "sp3", "species"));
        t.add_child(3, new_node(6, "s0", "sequence"));
        t.add_child(3, new_node(7, "s 1", "sequence"));
        t.add_child(4, new_node(8, "s2|x", "sequence"));
        t.add_child(5, new_node(9, "s3", "sequence"));
        (t, vec![6, 7, 8, 9])
    }

    /// Run the pruning step with counts 10, 10, 30, 5 per column at rho 15: species 3
    /// collapses over seqs 6, 7; genus 2 is pruned with everything beneath it.
    pub fn fitted() -> (TaxTree, Vec<NodeID>, TaxTree, Partition, EmFit) {
        let (stored, tips) = stored();
        let mut fitted = stored.clone();
        let mut part = Partition::new(tips.clone());
        subtree_counts(&mut fitted, &part, &[10.0, 10.0, 30.0, 5.0]);
        let report = prune_and_collapse(&mut fitted, &mut part, 15.0);
        assert_eq!((report.pruned, report.collapsed), (3, 1));
        let mut node_pi = vec![0.0; 10];
        for (v, p) in [(0, 1.0), (1, 1.0), (3, 0.4), (4, 0.6), (8, 0.6)] {
            node_pi[v] = p;
        }
        let fit = EmFit {
            node_pi,
            pi_tilde: vec![0.2, 0.2, 0.6, 0.0],
            ll_trace: vec![-10.0, -8.0],
            pi_trace: vec![
                vec![(6, 0.2), (7, 0.2), (8, 0.5), (9, 0.1)],
                vec![(3, 0.4), (8, 0.6)],
            ],
            pruned: 1,
            collapsed: 1,
            root_collapsed: false,
        };
        (stored, tips, fitted, part, fit)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fixtures::fitted;

    fn status(tree: &ReportTree, id: NodeID) -> VertexStatus {
        tree.vertices.iter().find(|v| v.id == id).unwrap().status
    }

    #[test]
    fn from_fit_labels_every_vertex() {
        let (stored, tips, fit_tree, part, fit) = fitted();
        let rt = ReportTree::from_fit(&stored, &tips, &fit_tree, &part, &fit);
        let ids: Vec<NodeID> = rt.vertices.iter().map(|v| v.id).collect();
        assert_eq!(ids, vec![0, 1, 3, 6, 7, 4, 8, 2, 5, 9]);
        use VertexStatus::*;
        for (id, s) in [
            (0, Active),
            (1, Active),
            (3, Collapsed),
            (6, Pooled),
            (7, Pooled),
            (4, Active),
            (8, Active),
            (2, Pruned),
            (5, Pruned),
            (9, Pruned),
        ] {
            assert_eq!(status(&rt, id), s, "vertex {id}");
        }
        let v3 = &rt.vertices[2];
        assert_eq!((v3.contigs, v3.count, v3.proportion), (2, Some(20.0), 0.4));
        let v6 = &rt.vertices[3];
        assert_eq!(
            (v6.contigs, v6.count, v6.proportion, v6.conditional),
            (0, None, 0.2, None)
        );
        let v9 = &rt.vertices[9];
        assert_eq!((v9.contigs, v9.proportion), (0, 0.0));
        assert_eq!(rt.vertices[0].parent, None);
        assert_eq!(rt.vertices[0].conditional, Some(1.0));
        assert_eq!(rt.vertices[5].parent, Some(1));
    }

    #[test]
    fn to_tax_tree_is_the_fitted_tree() {
        let (stored, tips, fit_tree, part, fit) = fitted();
        let rt = ReportTree::from_fit(&stored, &tips, &fit_tree, &part, &fit);
        let rebuilt = rt.to_tax_tree();
        let mut want: Vec<NodeID> = fit_tree.get_node_ids().collect();
        want.sort_unstable();
        let mut got: Vec<NodeID> = rebuilt.get_node_ids().collect();
        got.sort_unstable();
        assert_eq!(got, want);
        assert_eq!(preorder(&rebuilt), preorder(&fit_tree));
        for v in preorder(&fit_tree) {
            assert_eq!(children(&rebuilt, v), children(&fit_tree, v));
            assert_eq!(label(&rebuilt, v), label(&fit_tree, v));
            assert_eq!(rank(&rebuilt, v), rank(&fit_tree, v));
            assert_eq!(count(&rebuilt, v), count(&fit_tree, v));
            assert_eq!(conditional(&rebuilt, v), conditional(&fit_tree, v));
        }
        let mut leaves: Vec<NodeID> = rebuilt.get_leaf_ids().collect();
        leaves.sort_unstable();
        assert_eq!(leaves, vec![3, 8]);
    }

    #[test]
    fn to_stored_tax_tree_is_the_pre_em_tree() {
        let (stored, tips, fit_tree, part, fit) = fitted();
        let rt = ReportTree::from_fit(&stored, &tips, &fit_tree, &part, &fit);
        let rebuilt = rt.to_stored_tax_tree().unwrap();
        assert_eq!(preorder(&rebuilt), preorder(&stored));
        for v in preorder(&stored) {
            assert_eq!(children(&rebuilt, v), children(&stored, v));
            assert_eq!(label(&rebuilt, v), label(&stored, v));
            assert_eq!(rank(&rebuilt, v), rank(&stored, v));
            assert_eq!(count(&rebuilt, v), 0.0);
            assert_eq!(
                rebuilt.get_node(v).unwrap().get_weight(),
                if v == ROOT { Some(1.0) } else { None }
            );
        }
        // A tree whose root is not vertex 0 cannot feed the EM.
        let mut shifted = rt.clone();
        for v in shifted.vertices.iter_mut() {
            v.id += 1;
            v.parent = v.parent.map(|p| p + 1);
        }
        assert!(shifted.to_stored_tax_tree().is_err());
        assert!(ReportTree::default().to_stored_tax_tree().is_err());
    }

    #[test]
    fn em_leaf_of_pools_under_collapsed() {
        let (stored, tips, fit_tree, part, fit) = fitted();
        let rt = ReportTree::from_fit(&stored, &tips, &fit_tree, &part, &fit);
        let owner = rt.em_leaf_of();
        assert_eq!(owner[&6], 3);
        assert_eq!(owner[&7], 3);
        assert_eq!(owner[&3], 3);
        assert_eq!(owner[&8], 8);
        assert_eq!(owner[&9], NO_NODE);
        assert_eq!(owner[&0], 0);
    }

    #[test]
    fn refmap_lookup_is_by_seqid() {
        let m = RefMap::new(&[SeqId::new(5), SeqId::new(2)], &[9, 8]);
        assert_eq!(m.entries, vec![(SeqId::new(2), 8), (SeqId::new(5), 9)]);
        assert_eq!(m.tip_of(SeqId::new(5)), Some(9));
        assert_eq!(m.tip_of(SeqId::new(3)), None);
    }
}
