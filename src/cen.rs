// Census of MMCS subtrees: dead subtrees and repeats of completed restricted states.
//
// Legend:
//   ev, un, cr, ca   as in srch: edges, uncovered edges, critical edges per owner, candidates
//   b                remaining cardinality budget (None when unbounded)
//   U'               min {ev[e] ∩ ca : e ∈ un}      O_j  min {ev[e] ∩ ca : e ∈ cr[j]}
//   Tab              completed keys (sx exact, sw WL) and open repeat roots on the stack (ox, ow)
//
// The completions T of a node are the minimal transversals of U' that, for every owner j, miss
// some member of O_j. So the state is keyed by (b, U', multiset of O_j), dropping owners with an
// empty member (always safe). Owner identities, |s| and vertices outside ca matter only for
// labelled output. Keys enter the table when their subtree is complete, so an open ancestor never
// matches; a repeat nested inside an open repeat is not counted again.

use crate::bits::Bs;
use std::collections::HashSet;
use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};

#[derive(Default, Debug)]
pub struct Cen {
    pub dead: u64,
    pub ex_roots: u64,
    pub ex_nodes: u64,
    pub ex_dead: u64,
    pub wl_roots: u64,
    pub wl_nodes: u64,
    pub wl_dead: u64,
    pub ex_keys: u64,
    pub wl_keys: u64,
}

#[derive(Default)]
pub(crate) struct Tab {
    sx: HashSet<u128>,
    sw: HashSet<u128>,
    ox: u32,
    ow: u32,
    c: Cen,
}

impl Tab {
    /// Enters a node with keys (exact, WL); returns whether each is a counted repeat root.
    pub(crate) fn open(&mut self, (kx, kw): (u128, u128)) -> (bool, bool) {
        let rx = self.ox == 0 && self.sx.contains(&kx);
        let rw = self.ow == 0 && self.sw.contains(&kw);
        self.ox += rx as u32;
        self.ow += rw as u32;
        (rx, rw)
    }

    /// Leaves a completed node whose subtree has sub nodes and emitted nothing when dead.
    pub(crate) fn close(&mut self, (kx, kw): (u128, u128), (rx, rw): (bool, bool), sub: u64, dead: bool) {
        let c = &mut self.c;
        c.dead += dead as u64;
        let dn = if dead { sub } else { 0 };
        if rx {
            self.ox -= 1;
            c.ex_roots += 1;
            c.ex_nodes += sub;
            c.ex_dead += dn;
        }
        if rw {
            self.ow -= 1;
            c.wl_roots += 1;
            c.wl_nodes += sub;
            c.wl_dead += dn;
        }
        self.sx.insert(kx);
        self.sw.insert(kw);
    }

    pub(crate) fn finish(self) -> Cen {
        Cen { ex_keys: self.sx.len() as u64, wl_keys: self.sw.len() as u64, ..self.c }
    }
}

/// (exact, WL) keys of the restricted state of a node over n vertices.
pub(crate) fn keys(n: usize, ev: &[Bs], un: &Bs, cr: &[Bs], ca: &Bs, b: Option<usize>) -> (u128, u128) {
    let u = minimal(un.iter().map(|e| ev[e].and(ca)).collect(), n);
    let mut ow: Vec<Vec<Bs>> = vec![];
    for c in cr {
        let o: Vec<Bs> = c.iter().map(|e| ev[e].and(ca)).collect();
        if !o.iter().any(Bs::is0) {
            ow.push(minimal(o, n));
        }
    }
    ow.sort();
    (h128(&(b, &u, &ow)), wl(b, &u, &ow, n))
}

/// Inclusion-minimal members, deduplicated and sorted.
fn minimal(mut v: Vec<Bs>, n: usize) -> Vec<Bs> {
    v.sort_by_key(|b| (b.cnt(), b.clone()));
    v.dedup();
    if v.first().is_some_and(Bs::is0) {
        return vec![v.swap_remove(0)];
    }
    let mut by: Vec<Vec<usize>> = vec![vec![]; n];
    let mut o: Vec<Bs> = vec![];
    for b in v {
        if b.iter().any(|x| by[x].iter().any(|&k| o[k].sub(&b))) {
            continue;
        }
        by[b.first().unwrap()].push(o.len());
        o.push(b);
    }
    o.sort();
    o
}

fn h64<T: Hash>(t: &T, salt: u8) -> u64 {
    let mut h = DefaultHasher::new();
    salt.hash(&mut h);
    t.hash(&mut h);
    h.finish()
}

fn h128<T: Hash>(t: &T) -> u128 {
    (h64(t, 0) as u128) << 64 | h64(t, 1) as u128
}

/// Colour-refinement invariant of the coloured incidence structure: vertices, U' sets, owners and
/// owner sets. Equal invariants are necessary for isomorphism, so WL repeats bound isomorphic
/// repeats from above.
fn wl(b: Option<usize>, u: &[Bs], ow: &[Vec<Bs>], n: usize) -> u128 {
    // Node ids: 0..n vertices, then U' sets, then per owner: the owner node and its sets.
    let mut adj: Vec<Vec<usize>> = vec![vec![]; n];
    let mut col: Vec<u64> = vec![0; n];
    let add = |adj: &mut Vec<Vec<usize>>, col: &mut Vec<u64>, c: u64| {
        adj.push(vec![]);
        col.push(c);
        adj.len() - 1
    };
    for s in u {
        let x = add(&mut adj, &mut col, 1);
        for v in s.iter() {
            adj[x].push(v);
            adj[v].push(x);
        }
    }
    for o in ow {
        let y = add(&mut adj, &mut col, 2);
        for s in o {
            let x = add(&mut adj, &mut col, 3);
            adj[x].push(y);
            adj[y].push(x);
            for v in s.iter() {
                adj[x].push(v);
                adj[v].push(x);
            }
        }
    }
    // Vertices outside every set carry no information.
    let live: Vec<usize> = (0..adj.len()).filter(|&x| x >= n || !adj[x].is_empty()).collect();
    let classes = |col: &Vec<u64>| live.iter().map(|&x| col[x]).collect::<HashSet<_>>().len();
    let mut k0 = classes(&col);
    for _ in 0..live.len() {
        col = (0..adj.len())
            .map(|x| {
                let mut ns: Vec<u64> = adj[x].iter().map(|&y| col[y]).collect();
                ns.sort_unstable();
                h64(&(col[x], ns), 2)
            })
            .collect();
        let k1 = classes(&col);
        if k1 == k0 {
            break;
        }
        k0 = k1;
    }
    let mut ms: Vec<u64> = live.iter().map(|&x| col[x]).collect();
    ms.sort_unstable();
    h128(&(b, ms))
}
