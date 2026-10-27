// Input, normalization and static twins.
//
// Legend:
//   Raw    parsed input: ids (dense vertex order), lab labels, rows (raw edges as dense vertex
//          lists, eligibility applied), dup (repeated ids removed inside edges), fmt (input format)
//   H      hypergraph: n vertices, ev[e] vertex set of edge e, src[e] raw rows equal to e,
//          dom (raw row, kept edge) for rows dropped as strict supersets, empty (an edge is empty)
//   Tw     twin reduction: cls[c] dense members of class c, h hypergraph over classes
//   n, m   vertex and edge counts
//
// Tr(H) is unchanged by deleting duplicate and superset edges. An empty edge makes Tr(H) empty;
// no edges makes Tr(H) = {∅}. Twins (equal incidence) never both occur in a minimal hitting set,
// and swapping one for another preserves minimality, so Tr(H) is the product expansion of
// Tr(H / twins).

use crate::bits::Bs;
use serde_json::Value;
use std::collections::{BTreeSet, HashMap};

pub const SCHEMA: &str = "mcs-hg/1";

pub struct Raw {
    pub ids: Vec<String>,
    pub lab: Vec<Option<String>>,
    pub rows: Vec<Vec<usize>>,
    pub dup: usize,
    pub fmt: &'static str,
}

fn strs(v: &Value, what: &str) -> Result<Vec<String>, String> {
    v.as_array()
        .ok_or(format!("{what}: expected an array"))?
        .iter()
        .map(|x| x.as_str().map(str::to_owned).ok_or(format!("{what}: expected string ids")))
        .collect()
}

pub fn parse(b: &[u8]) -> Result<Raw, String> {
    let v: Value = serde_json::from_slice(b).map_err(|e| format!("malformed JSON: {e}"))?;
    let (ids, lab, er, el, fmt) = if let Some(s) = v.get("schema") {
        if s.as_str() != Some(SCHEMA) {
            return Err(format!("unsupported schema {s}; expected {SCHEMA}"));
        }
        let mut ids: Vec<String> = vec![];
        let mut lab: Vec<Option<String>> = vec![];
        let mut ix: HashMap<String, usize> = HashMap::new();
        for x in v["vertices"].as_array().ok_or("vertices: expected an array")? {
            let id = x["id"].as_str().ok_or("vertex without string id")?.to_owned();
            let l = match x.get("label") {
                None | Some(Value::Null) => None,
                Some(l) => Some(l.as_str().ok_or("label: expected a string")?.to_owned()),
            };
            match ix.get(&id) {
                Some(&i) if lab[i] != l => return Err(format!("conflicting labels for vertex {id}")),
                Some(_) => {}
                None => {
                    ix.insert(id.clone(), ids.len());
                    ids.push(id);
                    lab.push(l);
                }
            }
        }
        let look = |s: &String| ix.get(s).copied().ok_or(format!("unknown vertex id {s}"));
        let mut er = vec![];
        for e in v["edges"].as_array().ok_or("edges: expected an array")? {
            er.push(strs(e, "edge")?.iter().map(look).collect::<Result<Vec<_>, _>>()?);
        }
        let el = match v.get("eligible") {
            None | Some(Value::Null) => None,
            Some(x) => Some(strs(x, "eligible")?.iter().map(look).collect::<Result<BTreeSet<_>, _>>()?),
        };
        (ids, lab, er, el, "mcs-hg/1")
    } else if let Some(s) = v.get("sets") {
        // Publisher format: positive integer ids, optional {"id": label} map.
        let mut rows: Vec<Vec<u64>> = vec![];
        for e in s.as_array().ok_or("sets: expected an array")? {
            let r = e.as_array().ok_or("set: expected an array")?;
            rows.push(r.iter().map(|x| x.as_u64().ok_or(format!("malformed vertex number {x}"))).collect::<Result<_, _>>()?);
        }
        let mut ls: HashMap<u64, String> = HashMap::new();
        if let Some(l) = v.get("labels") {
            for (k, x) in l.as_object().ok_or("labels: expected an object")? {
                let k: u64 = k.parse().map_err(|_| format!("malformed label key {k}"))?;
                ls.insert(k, x.as_str().ok_or("label: expected a string")?.to_owned());
            }
        }
        let all: BTreeSet<u64> = rows.iter().flatten().copied().chain(ls.keys().copied()).collect();
        let ix: HashMap<u64, usize> = all.iter().enumerate().map(|(i, &x)| (x, i)).collect();
        let ids = all.iter().map(u64::to_string).collect();
        let lab = all.iter().map(|x| ls.get(x).cloned()).collect();
        let er = rows.iter().map(|r| r.iter().map(|x| ix[x]).collect()).collect();
        (ids, lab, er, None, "sets")
    } else {
        return Err("input has neither a schema field nor publisher sets".into());
    };
    let mut dup = 0;
    let rows = er
        .into_iter()
        .map(|r: Vec<usize>| {
            let s: BTreeSet<usize> = r.iter().copied().collect();
            dup += r.len() - s.len();
            s.into_iter().filter(|i| el.as_ref().is_none_or(|el| el.contains(i))).collect()
        })
        .collect();
    Ok(Raw { ids, lab, rows, dup, fmt })
}

#[derive(Clone)]
pub struct H {
    pub n: usize,
    pub ev: Vec<Bs>,
    pub src: Vec<Vec<usize>>,
    pub dom: Vec<(usize, usize)>,
    pub empty: bool,
}

pub fn norm(r: &Raw) -> H {
    let n = r.ids.len();
    let ez: Vec<usize> = (0..r.rows.len()).filter(|&i| r.rows[i].is_empty()).collect();
    if !ez.is_empty() {
        let dom = (0..r.rows.len()).filter(|i| !r.rows[*i].is_empty()).map(|i| (i, 0)).collect();
        return H { n, ev: vec![Bs::new(n)], src: vec![ez], dom, empty: true };
    }
    // Unique edges in order of first raw occurrence.
    let mut ix: HashMap<Bs, usize> = HashMap::new();
    let mut ue: Vec<Bs> = vec![];
    let mut us: Vec<Vec<usize>> = vec![];
    for (i, row) in r.rows.iter().enumerate() {
        let b = Bs::of(n, row.iter().copied());
        match ix.get(&b) {
            Some(&j) => us[j].push(i),
            None => {
                ix.insert(b.clone(), ue.len());
                ue.push(b);
                us.push(vec![i]);
            }
        }
    }
    // Minimalize: test each edge against kept edges filed under their least vertex.
    let mut ord: Vec<usize> = (0..ue.len()).collect();
    ord.sort_by_key(|&j| (ue[j].cnt(), j));
    let mut by: Vec<Vec<usize>> = vec![vec![]; n];
    let mut keep = vec![false; ue.len()];
    let mut dby = vec![usize::MAX; ue.len()];
    for &j in &ord {
        let d = ue[j].iter().find_map(|v| by[v].iter().copied().find(|&k| ue[k].sub(&ue[j])));
        match d {
            Some(k) => dby[j] = k,
            None => {
                keep[j] = true;
                by[ue[j].first().unwrap()].push(j);
            }
        }
    }
    let kj: Vec<usize> = (0..ue.len()).filter(|&j| keep[j]).collect();
    let pos: HashMap<usize, usize> = kj.iter().enumerate().map(|(p, &j)| (j, p)).collect();
    let mut dom = vec![];
    for j in 0..ue.len() {
        if !keep[j] {
            for &i in &us[j] {
                dom.push((i, pos[&dby[j]]));
            }
        }
    }
    H { n, ev: kj.iter().map(|&j| ue[j].clone()).collect(), src: kj.iter().map(|&j| us[j].clone()).collect(), dom, empty: false }
}

pub struct Tw {
    pub cls: Vec<Vec<usize>>,
    pub h: H,
}

/// Static twin classes. Vertices in no edge are dropped: they occur in no minimal hitting set.
pub fn twins(h: &H, on: bool) -> Tw {
    let m = h.ev.len();
    let mut inc = vec![Bs::new(m); h.n];
    for (e, b) in h.ev.iter().enumerate() {
        for v in b.iter() {
            inc[v].set(e);
        }
    }
    let mut cls: Vec<Vec<usize>> = vec![];
    let mut ix: HashMap<Bs, usize> = HashMap::new();
    for v in 0..h.n {
        if inc[v].is0() {
            continue;
        }
        let c = if on { ix.get(&inc[v]).copied() } else { None };
        match c {
            Some(c) => cls[c].push(v),
            None => {
                ix.insert(inc[v].clone(), cls.len());
                cls.push(vec![v]);
            }
        }
    }
    let k = cls.len();
    let mut cof = vec![usize::MAX; h.n];
    for (c, ms) in cls.iter().enumerate() {
        for &v in ms {
            cof[v] = c;
        }
    }
    let ev = h.ev.iter().map(|b| Bs::of(k, b.iter().map(|v| cof[v]))).collect();
    Tw { cls, h: H { n: k, ev, src: h.src.clone(), dom: h.dom.clone(), empty: h.empty } }
}

/// Calls f on every expansion of a class-level solution s (sorted dense vertex ids).
pub fn expand(cls: &[Vec<usize>], s: &[usize], f: &mut dyn FnMut(&[usize]) -> bool) -> bool {
    let mut ix = vec![0usize; s.len()];
    let mut o: Vec<usize> = Vec::with_capacity(s.len());
    loop {
        o.clear();
        o.extend(s.iter().zip(&ix).map(|(&c, &i)| cls[c][i]));
        o.sort_unstable();
        if !f(&o) {
            return false;
        }
        let mut p = 0;
        loop {
            if p == s.len() {
                return true;
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
