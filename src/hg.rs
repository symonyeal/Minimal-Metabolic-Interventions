// Normalization, static twins and expansion.
//
// Legend:
//   rows   raw edges as sorted dense vertex lists, eligibility applied
//   H      hypergraph: n vertices, ev[e] vertex set of edge e
//   Tw     twin reduction: cls[c] dense members of class c, h hypergraph over classes
//   n, m   vertex and edge counts
//
// Tr(H) is unchanged by deleting duplicate and superset edges. An empty edge makes Tr(H) empty;
// no edges makes Tr(H) = {∅}. Twins (equal incidence) never both occur in a minimal hitting set,
// and swapping one for another preserves minimality, so Tr(H) is the product expansion of
// Tr(H / twins).

use crate::bits::Bs;
use std::collections::{HashMap, HashSet};

pub struct H {
    pub n: usize,
    pub ev: Vec<Bs>,
}

/// Inclusion-minimal distinct rows, in order of first raw occurrence, and the number of rows
/// merged as duplicates. Superset removal files kept edges under their least vertex; worst case
/// O(m² n / 64).
pub fn norm(n: usize, rows: &[Vec<usize>]) -> (H, usize) {
    let mut seen: HashSet<Bs> = HashSet::new();
    let mut ue: Vec<Bs> = vec![];
    for row in rows {
        let b = Bs::of(n, row.iter().copied());
        if seen.insert(b.clone()) {
            ue.push(b);
        }
    }
    let merged = rows.len() - ue.len();
    if let Some(z) = ue.iter().find(|b| b.is0()) {
        return (H { n, ev: vec![z.clone()] }, merged);
    }
    let mut ord: Vec<usize> = (0..ue.len()).collect();
    ord.sort_by_key(|&j| (ue[j].cnt(), j));
    let mut by: Vec<Vec<usize>> = vec![vec![]; n];
    let mut keep = vec![false; ue.len()];
    for &j in &ord {
        if !ue[j].iter().any(|v| by[v].iter().any(|&k| ue[k].sub(&ue[j]))) {
            keep[j] = true;
            by[ue[j].first().unwrap()].push(j);
        }
    }
    let ev = ue.into_iter().zip(keep).filter_map(|(b, k)| k.then_some(b)).collect();
    (H { n, ev }, merged)
}

pub struct Tw {
    pub cls: Vec<Vec<usize>>,
    pub h: H,
}

/// Static twin classes (singletons when `on` is false). Vertices in no edge are dropped: they
/// occur in no minimal hitting set.
pub fn twins(h: &H, on: bool) -> Tw {
    let m = h.ev.len();
    let mut inc = vec![Bs::new(m); h.n];
    for (e, b) in h.ev.iter().enumerate() {
        for v in b.iter() {
            inc[v].set(e);
        }
    }
    let mut cls: Vec<Vec<usize>> = vec![];
    let mut ix: HashMap<&Bs, usize> = HashMap::new();
    let mut cof = vec![usize::MAX; h.n];
    for v in 0..h.n {
        if inc[v].is0() {
            continue;
        }
        let c = match ix.get(&inc[v]) {
            Some(&c) if on => c,
            _ => {
                ix.insert(&inc[v], cls.len());
                cls.push(vec![]);
                cls.len() - 1
            }
        };
        cls[c].push(v);
        cof[v] = c;
    }
    let k = cls.len();
    let ev = h.ev.iter().map(|b| Bs::of(k, b.iter().map(|v| cof[v]))).collect();
    Tw { cls, h: H { n: k, ev } }
}

/// Calls f on every expansion of a class-level solution s, each as sorted dense vertex ids.
pub fn expand(cls: &[Vec<usize>], s: &[usize], f: &mut dyn FnMut(&[usize])) {
    let mut ix = vec![0usize; s.len()];
    let mut o: Vec<usize> = Vec::with_capacity(s.len());
    loop {
        o.clear();
        o.extend(s.iter().zip(&ix).map(|(&c, &i)| cls[c][i]));
        o.sort_unstable();
        f(&o);
        let mut p = 0;
        loop {
            if p == s.len() {
                return;
            }
            ix[p] += 1;
            if ix[p] < cls[s[p]].len() {
                break;
            }
            ix[p] = 0;
            p += 1;
        }
    }
}
