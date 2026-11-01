// MMCS (Murakami–Uno 2014) over dense bitsets, with an optional census of its subtrees.
//
// Legend:
//   ev[e]  vertex set of edge e          inc[v]  edges containing vertex v
//   s      selected vertices (stack)     cr[j]   critical (private) edges of s[j]
//   un     uncovered edges               ca      candidate vertices
//   k      optional cardinality bound    Ord     candidate order: Idx ascending, Unhit by
//                                                increasing unhit degree (Hörmann–Schirneck 2026)
//   dl     optional deadline             St      counters; done is false when dl stopped the run
//
// The subtree of (s, ca) emits exactly the minimal hitting sets M with s ⊆ M ⊆ s ∪ ca, each once.

use crate::bits::Bs;
use crate::cen::{self, Cen, Tab};
use std::time::Instant;

#[derive(Clone, Copy, PartialEq, Debug)]
pub enum Ord {
    Idx,
    Unhit,
}

pub struct Cfg {
    pub ord: Ord,
    pub k: Option<usize>,
    pub census: bool,
    pub dl: Option<Instant>,
}

pub struct St {
    pub nodes: u64,
    pub sols: u64,
    pub done: bool,
    pub cen: Cen,
}

struct Sr<'a> {
    n: usize,
    ev: &'a [Bs],
    inc: Vec<Bs>,
    s: Vec<usize>,
    cr: Vec<Bs>,
    un: Bs,
    cfg: &'a Cfg,
    st: St,
    tab: Tab,
    emit: &'a mut dyn FnMut(&[usize]),
}

/// Enumerates Tr(H) of the hypergraph (n, ev); emit receives each solution as sorted vertex ids.
pub fn run(n: usize, ev: &[Bs], cfg: &Cfg, emit: &mut dyn FnMut(&[usize])) -> St {
    let m = ev.len();
    let mut inc = vec![Bs::new(m); n];
    for (e, b) in ev.iter().enumerate() {
        for v in b.iter() {
            inc[v].set(e);
        }
    }
    let st = St { nodes: 0, sols: 0, done: true, cen: Cen::default() };
    let mut r = Sr { n, ev, inc, s: vec![], cr: vec![], un: Bs::of(m, 0..m), cfg, st, tab: Tab::default(), emit };
    r.go(&mut Bs::of(n, 0..n));
    r.st.cen = r.tab.finish();
    r.st
}

impl Sr<'_> {
    // Node cost ranges from nanoseconds to seconds (wide inputs, checked output), so the deadline
    // is read at every node and after every emitted solution; a clock read is far below node cost.
    fn late(&mut self) -> bool {
        if self.cfg.dl.is_some_and(|d| Instant::now() >= d) {
            self.st.done = false;
        }
        !self.st.done
    }

    fn out(&mut self) -> bool {
        self.st.sols += 1;
        let mut o = self.s.clone();
        o.sort_unstable();
        (self.emit)(&o);
        !self.late()
    }

    /// Removes inc[v] from every cr[j]; returns the removed parts and whether some owner died.
    fn hit(&mut self, v: usize) -> (Vec<(usize, Bs)>, bool) {
        let mut sv = vec![];
        for j in 0..self.s.len() {
            if self.cr[j].meets(&self.inc[v]) {
                let rm = self.cr[j].and(&self.inc[v]);
                self.cr[j].andnot_with(&self.inc[v]);
                let z = self.cr[j].is0();
                sv.push((j, rm));
                if z {
                    return (sv, true);
                }
            }
        }
        (sv, false)
    }

    fn go(&mut self, ca: &mut Bs) -> bool {
        self.st.nodes += 1;
        if self.late() {
            return false;
        }
        let b = self.cfg.k.map(|k| k - self.s.len());
        let key = self.cfg.census.then(|| cen::keys(self.n, self.ev, &self.un, &self.cr, ca, b));
        let rep = key.map(|k| self.tab.open(k));
        let (n0, s0) = (self.st.nodes, self.st.sols);
        if !self.body(ca) {
            return false;
        }
        if let (Some(k), Some(r)) = (key, rep) {
            self.tab.close(k, r, self.st.nodes - n0 + 1, self.st.sols == s0);
        }
        true
    }

    fn body(&mut self, ca: &mut Bs) -> bool {
        if self.un.is0() {
            return self.out();
        }
        if self.cfg.k.is_some_and(|k| self.s.len() >= k) {
            return true;
        }
        // Branch on the uncovered edge with fewest candidates; an edge with none ends the node.
        let (mut f, mut bc) = (usize::MAX, u32::MAX);
        for e in self.un.iter() {
            let c = self.ev[e].and_cnt(ca);
            if c < bc {
                (f, bc) = (e, c);
                if c == 0 {
                    return true;
                }
            }
        }
        let mut cs: Vec<usize> = self.ev[f].and(ca).iter().collect();
        if self.cfg.ord == Ord::Unhit {
            cs.sort_by_key(|&v| (self.inc[v].and_cnt(&self.un), v));
        }
        for &v in &cs {
            ca.clr(v);
        }
        for &v in &cs {
            let (sv, viol) = self.hit(v);
            let mut ok = true;
            if !viol {
                let nc = self.inc[v].and(&self.un);
                self.un.andnot_with(&nc);
                self.s.push(v);
                self.cr.push(nc);
                ok = self.go(ca);
                let nc = self.cr.pop().unwrap();
                self.s.pop();
                self.un.or_with(&nc);
            }
            for (j, rm) in sv {
                self.cr[j].or_with(&rm);
            }
            if !ok {
                return false;
            }
            // A violator stays out of later siblings' candidates: it violates every superset of s.
            if !viol {
                ca.set(v);
            }
        }
        for &v in &cs {
            ca.set(v);
        }
        true
    }
}
