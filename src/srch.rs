// MMCS and RS (Murakami–Uno 2014) over dense bitsets, with the G2 census.
//
// Legend:
//   ev[e]  vertex set of edge e          inc[v]  edges containing vertex v
//   s      selected vertices (stack)     cr[j]   critical (private) edges of s[j]
//   un     uncovered edges               ca      MMCS candidate vertices
//   k      optional cardinality bound    Ord     candidate order: Idx ascending, Unhit by
//                                                increasing unhit degree (Hörmann–Schirneck 2026)
//   St     counters and stop state       Cen     census counters
//
// MMCS subtree of (s, ca) emits exactly the minimal hitting sets M with s ⊆ M ⊆ s ∪ ca. Its
// completions T are the minimal transversals of U' = {e ∩ ca : e ∈ un} that, for every owner
// s[j], miss some member of O_j = {e ∩ ca : e ∈ cr[j]}. So the census keys a state by
// (remaining budget, min U', multiset of min O_j), dropping owners with an empty member (always
// safe). Owner identities, |s| and vertices outside ca matter only for labelled output.

use crate::bits::Bs;
use std::collections::HashSet;
use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};
use std::time::Instant;

#[derive(Clone, Copy, PartialEq, Debug)]
pub enum Ord {
    Idx,
    Unhit,
}
#[derive(Clone, Copy, PartialEq, Debug)]
pub enum Alg {
    Mmcs,
    Rs,
}
#[derive(Clone, Copy, PartialEq, Debug)]
pub enum Stop {
    Done,
    Cap,
    Time,
}

pub struct Cfg {
    pub alg: Alg,
    pub ord: Ord,
    pub k: Option<usize>,
    pub census: bool,
    pub cap: Option<u64>,
    pub dl: Option<Instant>,
}

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

pub struct St {
    pub nodes: u64,
    pub sols: u64,
    pub stop: Stop,
    pub cen: Cen,
}

struct Sr<'a> {
    ev: &'a [Bs],
    inc: Vec<Bs>,
    n: usize,
    s: Vec<usize>,
    cr: Vec<Bs>,
    un: Bs,
    cfg: &'a Cfg,
    st: St,
    emit: &'a mut dyn FnMut(&[usize]) -> bool,
    sx: HashSet<u128>,
    sw: HashSet<u128>,
    ix: u32,
    iw: u32,
}

pub fn run(n: usize, ev: &[Bs], cfg: &Cfg, emit: &mut dyn FnMut(&[usize]) -> bool) -> St {
    let m = ev.len();
    let mut inc = vec![Bs::new(m); n];
    for (e, b) in ev.iter().enumerate() {
        for v in b.iter() {
            inc[v].set(e);
        }
    }
    let mut un = Bs::new(m);
    (0..m).for_each(|e| un.set(e));
    let st = St { nodes: 0, sols: 0, stop: Stop::Done, cen: Cen::default() };
    let mut r = Sr { ev, inc, n, s: vec![], cr: vec![], un, cfg, st, emit, sx: HashSet::new(), sw: HashSet::new(), ix: 0, iw: 0 };
    match cfg.alg {
        Alg::Mmcs => {
            let mut ca = Bs::of(n, 0..n);
            r.go(&mut ca);
        }
        Alg::Rs => {
            r.rs();
        }
    }
    r.st.cen.ex_keys = r.sx.len() as u64;
    r.st.cen.wl_keys = r.sw.len() as u64;
    r.st
}

impl Sr<'_> {
    // Node cost ranges from nanoseconds to seconds (wide inputs, checked output), so the deadline
    // is read at every node and after every emitted solution; a clock read is far below node cost.
    fn late(&mut self) -> bool {
        if self.cfg.dl.is_some_and(|d| Instant::now() >= d) {
            self.st.stop = Stop::Time;
            return true;
        }
        false
    }

    fn tick(&mut self) -> bool {
        self.st.nodes += 1;
        !self.late()
    }

    fn out(&mut self) -> bool {
        self.st.sols += 1;
        let mut o = self.s.clone();
        o.sort_unstable();
        if !(self.emit)(&o) || self.cfg.cap.is_some_and(|c| self.st.sols >= c) {
            self.st.stop = Stop::Cap;
            return false;
        }
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
        if !self.tick() {
            return false;
        }
        if !self.cfg.census {
            return self.body(ca);
        }
        let (n0, s0) = (self.st.nodes, self.st.sols);
        let (kx, kw) = self.keys(ca);
        let rx = self.ix == 0 && self.sx.contains(&kx);
        let rw = self.iw == 0 && self.sw.contains(&kw);
        self.ix += rx as u32;
        self.iw += rw as u32;
        if !self.body(ca) {
            return false;
        }
        let sub = self.st.nodes - n0 + 1;
        let dead = self.st.sols == s0;
        let c = &mut self.st.cen;
        c.dead += dead as u64;
        if rx {
            self.ix -= 1;
            c.ex_roots += 1;
            c.ex_nodes += sub;
            c.ex_dead += if dead { sub } else { 0 };
        }
        if rw {
            self.iw -= 1;
            c.wl_roots += 1;
            c.wl_nodes += sub;
            c.wl_dead += if dead { sub } else { 0 };
        }
        // Keys enter the table when their subtree is complete, so an open ancestor never matches.
        self.sx.insert(kx);
        self.sw.insert(kw);
        true
    }

    fn body(&mut self, ca: &mut Bs) -> bool {
        if self.un.is0() {
            return self.out();
        }
        if self.cfg.k.is_some_and(|k| self.s.len() >= k) {
            return true;
        }
        let (mut f, mut bc) = (usize::MAX, u32::MAX);
        for e in self.un.iter() {
            let c = self.ev[e].and_cnt(ca);
            if c < bc {
                (f, bc) = (e, c);
                if c == 0 {
                    break;
                }
            }
        }
        if bc == 0 {
            return true;
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

    fn rs(&mut self) -> bool {
        if !self.tick() {
            return false;
        }
        if self.un.is0() {
            return self.out();
        }
        if self.cfg.k.is_some_and(|k| self.s.len() >= k) {
            return true;
        }
        let i = self.un.first().unwrap();
        let vs: Vec<usize> = self.ev[i].iter().collect();
        for v in vs {
            let (sv, viol) = self.hit(v);
            let mut ok = true;
            // Child test: every earlier owner keeps a critical edge of index below i.
            if !viol && self.cr.iter().all(|c| c.first().is_some_and(|t| t < i)) {
                let nc = self.inc[v].and(&self.un);
                self.un.andnot_with(&nc);
                self.s.push(v);
                self.cr.push(nc);
                ok = self.rs();
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
        }
        true
    }

    fn keys(&self, ca: &Bs) -> (u128, u128) {
        let u = minimal(self.un.iter().map(|e| self.ev[e].and(ca)).collect(), self.n);
        let mut ow: Vec<Vec<Bs>> = vec![];
        for c in &self.cr {
            let o: Vec<Bs> = c.iter().map(|e| self.ev[e].and(ca)).collect();
            if o.iter().any(Bs::is0) {
                continue;
            }
            ow.push(minimal(o, self.n));
        }
        ow.sort();
        let b = self.cfg.k.map(|k| k - self.s.len());
        (h128(&(b, &u, &ow)), wl(b, &u, &ow, self.n))
    }
}

/// Inclusion-minimal members, deduplicated and sorted.
pub fn minimal(mut v: Vec<Bs>, n: usize) -> Vec<Bs> {
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
        if let Some(x) = b.first() {
            by[x].push(o.len());
        }
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

/// Color-refinement invariant of the colored incidence structure: vertices, U' sets, owners and
/// owner sets. Equal invariants are necessary for isomorphism, so WL repeats upper-bound
/// isomorphic repeats.
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
        let nc: Vec<u64> = (0..adj.len())
            .map(|x| {
                let mut ns: Vec<u64> = adj[x].iter().map(|&y| col[y]).collect();
                ns.sort_unstable();
                h64(&(col[x], ns), 2)
            })
            .collect();
        col = nc;
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
