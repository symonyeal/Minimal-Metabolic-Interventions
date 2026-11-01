// Input formats: hypergraphs (mcs-hg/1 or publisher sets), solution files and BiGG models.
//
// Legend:
//   Raw    parsed hypergraph: ids (dense vertex order), rows (raw edges as sorted dense vertex
//          lists, eligibility applied), dup (repeated ids removed inside edges), fmt (format)

use crate::die;
use mcs::net::{Net, Rxn};
use serde_json::{Value, json};
use std::collections::{BTreeSet, HashMap};

pub const SCHEMA: &str = "mcs-hg/1";

pub struct Raw {
    pub ids: Vec<String>,
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

/// Labels are validated (a repeated vertex must repeat its label) but not kept.
pub fn parse(b: &[u8]) -> Result<Raw, String> {
    let v: Value = serde_json::from_slice(b).map_err(|e| format!("malformed JSON: {e}"))?;
    let (ids, er, el, fmt) = if let Some(s) = v.get("schema") {
        if s.as_str() != Some(SCHEMA) {
            return Err(format!("unsupported schema {s}; expected {SCHEMA}"));
        }
        let mut ids: Vec<String> = vec![];
        let mut lab: Vec<Option<&str>> = vec![];
        let mut ix: HashMap<String, usize> = HashMap::new();
        for x in v["vertices"].as_array().ok_or("vertices: expected an array")? {
            let id = x["id"].as_str().ok_or("vertex without string id")?.to_owned();
            let l = match x.get("label") {
                None | Some(Value::Null) => None,
                Some(l) => Some(l.as_str().ok_or("label: expected a string")?),
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
        (ids, er, el, SCHEMA)
    } else if let Some(s) = v.get("sets") {
        // Publisher format: integer ids, optional {"id": label} map.
        let mut rows: Vec<Vec<u64>> = vec![];
        for e in s.as_array().ok_or("sets: expected an array")? {
            let r = e.as_array().ok_or("set: expected an array")?;
            rows.push(
                r.iter().map(|x| x.as_u64().ok_or(format!("malformed vertex number {x}"))).collect::<Result<_, _>>()?,
            );
        }
        let mut ls: Vec<u64> = vec![];
        if let Some(l) = v.get("labels") {
            for (k, x) in l.as_object().ok_or("labels: expected an object")? {
                ls.push(k.parse().map_err(|_| format!("malformed label key {k}"))?);
                x.as_str().ok_or("label: expected a string")?;
            }
        }
        let all: BTreeSet<u64> = rows.iter().flatten().copied().chain(ls).collect();
        let ix: HashMap<u64, usize> = all.iter().enumerate().map(|(i, &x)| (x, i)).collect();
        let er = rows.iter().map(|r| r.iter().map(|x| ix[x]).collect()).collect();
        (all.iter().map(u64::to_string).collect(), er, None, "sets")
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
    Ok(Raw { ids, rows, dup, fmt })
}

/// Reads and parses a hypergraph; invalid input is reported as JSON with exit status 3.
pub fn load(p: &str) -> (Vec<u8>, Raw) {
    let b = std::fs::read(p).unwrap_or_else(|e| die(&format!("{p}: {e}")));
    match parse(&b) {
        Ok(r) => (b, r),
        Err(e) => {
            println!("{}", json!({"search_status": "invalid_input", "error": e}));
            std::process::exit(3)
        }
    }
}

/// Solutions file: a JSON array of arrays, or JSON lines of arrays; ids as strings or integers.
pub fn sols(p: &str, ix: &HashMap<&str, usize>) -> Vec<Vec<usize>> {
    let t = std::fs::read_to_string(p).unwrap_or_else(|e| die(&format!("{p}: {e}")));
    let rows: Vec<Value> = match serde_json::from_str::<Value>(&t) {
        Ok(Value::Array(a)) if a.iter().all(Value::is_array) => a,
        _ => t
            .lines()
            .filter(|l| !l.trim().is_empty())
            .map(|l| serde_json::from_str(l).unwrap_or_else(|e| die(&format!("{p}: {e}"))))
            .collect(),
    };
    let id = |x: &Value| {
        let s = match x {
            Value::String(s) => s.clone(),
            Value::Number(n) => n.to_string(),
            _ => die("solution id: expected a string or integer"),
        };
        *ix.get(s.as_str()).unwrap_or_else(|| die(&format!("solution uses unknown id {s}")))
    };
    rows.iter()
        .map(|r| {
            let mut c: Vec<usize> =
                r.as_array().unwrap_or_else(|| die("solution: expected an array")).iter().map(id).collect();
            c.sort_unstable();
            c
        })
        .collect()
}

/// BiGG JSON model: metabolite ids and reactions with stoichiometry and flux bounds.
pub fn bigg(b: &[u8]) -> Result<Net, String> {
    let m: Value = serde_json::from_slice(b).map_err(|e| format!("malformed JSON: {e}"))?;
    let arr = |k: &str| m[k].as_array().ok_or(format!("{k}: expected an array"));
    let mets: Vec<String> = arr("metabolites")?
        .iter()
        .map(|x| x["id"].as_str().map(str::to_owned).ok_or("metabolite without id"))
        .collect::<Result<_, _>>()?;
    let mi: HashMap<&str, usize> = mets.iter().enumerate().map(|(i, k)| (k.as_str(), i)).collect();
    let mut rx = vec![];
    for r in arr("reactions")? {
        let id = r["id"].as_str().ok_or("reaction without id")?.to_owned();
        let f = |k: &str| r[k].as_f64().ok_or(format!("{id}: {k} missing"));
        let mut st = vec![];
        for (k, c) in r["metabolites"].as_object().ok_or(format!("{id}: metabolites missing"))? {
            st.push((
                *mi.get(k.as_str()).ok_or(format!("{id}: unknown metabolite {k}"))?,
                c.as_f64().ok_or(format!("{id}: malformed coefficient"))?,
            ));
        }
        let name = r["name"].as_str().filter(|s| !s.is_empty()).map(str::to_owned);
        rx.push(Rxn { lo: f("lower_bound")?, hi: f("upper_bound")?, id, name, st });
    }
    Ok(Net { mets, rx })
}

#[cfg(test)]
mod tests {
    use super::parse;
    use mcs::hg;

    #[test]
    fn input_contract() {
        assert!(parse(br#"{"schema":"mcs-hg/2","vertices":[],"edges":[]}"#).is_err());
        assert!(parse(br#"{"schema":"mcs-hg/1","vertices":[{"id":"a"}],"edges":[["b"]]}"#).is_err());
        assert!(
            parse(br#"{"schema":"mcs-hg/1","vertices":[{"id":"a","label":"x"},{"id":"a","label":"y"}],"edges":[]}"#)
                .is_err()
        );
        assert!(parse(br#"{"sets":[[1,2.5]]}"#).is_err());
        assert!(parse(br#"{"sets":[[1]],"labels":{"x":"y"}}"#).is_err());
        let r = parse(br#"{"schema":"mcs-hg/1","vertices":[{"id":"a"},{"id":"b"},{"id":"c"}],"edges":[["a","a","b"],["c"]],"eligible":["a","b"]}"#).unwrap();
        assert_eq!(r.dup, 1);
        assert_eq!(r.rows, vec![vec![0, 1], vec![]]);
        let (h, _) = hg::norm(r.ids.len(), &r.rows);
        assert!(h.ev.len() == 1 && h.ev[0].is0());
        let r = parse(br#"{"sets":[[3,1],[1],[1,3]],"labels":{"7":"x"}}"#).unwrap();
        assert_eq!(r.ids, vec!["1", "3", "7"]);
        let (h, merged) = hg::norm(r.ids.len(), &r.rows);
        assert_eq!((h.ev.len(), merged), (1, 1));
    }
}
