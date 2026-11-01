// Metabolic network → irreversible integer flux cone, elementary modes and cut-set scenarios.
//
// Legend:
//   Net    mets metabolite ids, rx reactions
//   Rxn    id, name, lo/hi flux bounds, st stoichiometry (metabolite index, coefficient)
//   Cone   cols (reaction, direction ±1) in reaction order, forward first; a integer matrix over
//          the metabolites in some column; pruned reactions removed before enumeration
//   sup    reaction support of a mode, sorted reaction indices
//
// The cone is {v : N v = 0, sign constraints from the bounds}: forward allowed iff hi > 0, reverse
// iff lo < 0; finite bound magnitudes are ignored, as in elementary-mode analysis. Reversible
// reactions are split into forward and reverse columns, which adds one two-cycle per reversible
// reaction; modes drops them. Each column of N is scaled by the lcm of the denominators of its
// decimal coefficients, which leaves supports unchanged.
//
// Blocked reactions. A direction of reaction j is infeasible when some metabolite i it touches has
// no other reaction able to contribute the opposite sign: steady state would fail at i. This rule,
// run to a fixpoint, is sound but not complete; reactions left with no feasible direction are
// pruned. Every extreme ray of the split cone other than a two-cycle is a flux vector, and every
// feasible flux decomposes into such rays, so the exactly blocked reactions are those in no mode.

use crate::dd::gcd;

pub struct Rxn {
    pub id: String,
    pub name: Option<String>,
    pub lo: f64,
    pub hi: f64,
    pub st: Vec<(usize, f64)>,
}

pub struct Net {
    pub mets: Vec<String>,
    pub rx: Vec<Rxn>,
}

pub struct Cone {
    pub cols: Vec<(usize, i8)>,
    pub a: Vec<Vec<i128>>,
    pub pruned: Vec<usize>,
}

/// Exact value p/q of the shortest decimal that round-trips to x.
fn dec(x: f64) -> Result<(i128, i128), String> {
    let s = format!("{x}");
    let (ip, fp) = s.trim_start_matches('-').split_once('.').unwrap_or((s.trim_start_matches('-'), ""));
    let bad = || format!("coefficient {s} exceeds i128");
    let p: i128 = format!("{ip}{fp}").parse().map_err(|_| bad())?;
    let q = 10i128.checked_pow(fp.len() as u32).ok_or_else(bad)?;
    let g = gcd(p, q);
    let p = if x < 0.0 { -p } else { p };
    Ok((p / g, q / g))
}

pub fn cone(net: &Net) -> Result<Cone, String> {
    let mut by: Vec<Vec<(usize, f64)>> = vec![vec![]; net.mets.len()];
    for (j, r) in net.rx.iter().enumerate() {
        for &(i, c) in r.st.iter().filter(|x| x.1 != 0.0) {
            by[i].push((j, c));
        }
    }
    // live[j] = [forward, reverse] feasible so far.
    let mut live: Vec<[bool; 2]> = net.rx.iter().map(|r| [r.hi > 0.0, r.lo < 0.0]).collect();
    let mut ch = true;
    while ch {
        ch = false;
        for inc in &by {
            for &(j, c) in inc {
                for (k, s) in [(0, c), (1, -c)] {
                    let opp = |&(o, co): &(usize, f64)| {
                        o != j && ((live[o][0] && co * s < 0.0) || (live[o][1] && co * s > 0.0))
                    };
                    if live[j][k] && !inc.iter().any(opp) {
                        live[j][k] = false;
                        ch = true;
                    }
                }
            }
        }
    }
    let pruned: Vec<usize> = (0..net.rx.len()).filter(|&j| live[j] == [false, false]).collect();
    let mut cols = vec![];
    for (j, r) in net.rx.iter().enumerate() {
        if live[j] != [false, false] {
            if r.hi > 0.0 {
                cols.push((j, 1));
            }
            if r.lo < 0.0 {
                cols.push((j, -1));
            }
        }
    }
    let mut a = vec![vec![0i128; cols.len()]; net.mets.len()];
    let mut used = vec![false; net.mets.len()];
    for (k, &(j, d)) in cols.iter().enumerate() {
        let st = &net.rx[j].st;
        let pq = st.iter().map(|&(_, c)| dec(c)).collect::<Result<Vec<_>, _>>()?;
        let l = pq.iter().try_fold(1i128, |l, &(_, q)| (l / gcd(l, q)).checked_mul(q)).ok_or("lcm exceeds i128")?;
        for (&(i, _), &(p, q)) in st.iter().zip(&pq) {
            a[i][k] = (p * (l / q)) * d as i128;
            used[i] = true;
        }
    }
    let a = a.into_iter().zip(used).filter_map(|(r, u)| u.then_some(r)).collect();
    Ok(Cone { cols, a, pruned })
}

pub struct Modes {
    /// Columns of each mode, sorted.
    pub cols: Vec<Vec<usize>>,
    pub sup: Vec<Vec<usize>>,
    pub two_cycles: usize,
    pub blocked: Vec<usize>,
}

/// Modes from the extreme rays of the cone, two-cycles dropped.
pub fn modes(net: &Net, c: &Cone, rays: Vec<Vec<usize>>) -> Modes {
    let mut m = Modes { cols: vec![], sup: vec![], two_cycles: 0, blocked: vec![] };
    let mut hit = vec![false; net.rx.len()];
    for r in rays {
        let mut s: Vec<usize> = r.iter().map(|&k| c.cols[k].0).collect();
        s.dedup();
        if s.len() == 1 && r.len() == 2 {
            m.two_cycles += 1;
            continue;
        }
        s.iter().for_each(|&j| hit[j] = true);
        m.cols.push(r);
        m.sup.push(s);
    }
    m.blocked = (0..net.rx.len()).filter(|&j| !hit[j]).collect();
    m
}

#[derive(Clone, Copy, PartialEq, Debug)]
pub enum Pol {
    /// Every reaction but the target.
    AllButTarget,
    /// No exchange, demand, sink, biomass or ATP maintenance reaction, and no reaction whose
    /// bounds exclude zero.
    Core,
}

/// Cut-set hypergraph: edges are the supports of the modes using t (every mode when t is None);
/// returns the edges, their vertices (sorted) and the eligible vertices.
pub fn scenario<'a>(
    net: &Net,
    sup: &'a [Vec<usize>],
    t: Option<usize>,
    pol: Pol,
) -> (Vec<&'a Vec<usize>>, Vec<usize>, Vec<usize>) {
    let es: Vec<&Vec<usize>> = sup.iter().filter(|s| t.is_none_or(|t| s.binary_search(&t).is_ok())).collect();
    let mut vs: Vec<usize> = es.iter().flat_map(|s| s.iter().copied()).collect();
    vs.sort_unstable();
    vs.dedup();
    let ok = |&j: &usize| {
        let r = &net.rx[j];
        match pol {
            Pol::AllButTarget => Some(j) != t,
            Pol::Core => {
                !["EX_", "DM_", "SK_", "BIOMASS"].iter().any(|p| r.id.starts_with(p))
                    && r.id != "ATPM"
                    && r.lo <= 0.0
                    && 0.0 <= r.hi
            }
        }
    };
    let el = vs.iter().copied().filter(ok).collect();
    (es, vs, el)
}
