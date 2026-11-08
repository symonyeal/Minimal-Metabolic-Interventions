// enumerate and check.
//
// Legend:
//   r     parsed input (inp::Raw)        h    normalized hypergraph     tw  twin reduction
//   rep   class-level solution            nexp expanded solution count
//   sm    summary JSON                    t0   start, t1 end of load and normalization
//
// Outputs are write-once. The summary accepts the solution bytes by SHA-256 and is published last.

use crate::inp::{self, Raw};
use crate::output::Pending;
use crate::{Args, die, hex};
use mcs::chk;
use mcs::hg::{self, H};
use mcs::srch::{self, Cfg, Ord, St};
use serde_json::{Value, json};
use std::collections::{HashMap, HashSet};
use std::io::{BufWriter, Write};
use std::time::{Duration, Instant};

/// Input description; sha256_normalized hashes the sorted original ids of each normalized edge.
pub fn info(p: &str, b: &[u8], r: &Raw, h: &H, merged: usize) -> Value {
    let es: Vec<Vec<&str>> = h.ev.iter().map(|e| e.iter().map(|v| r.ids[v].as_str()).collect()).collect();
    json!({
        "path": p, "format": r.fmt, "sha256_bytes": hex(b), "sha256_normalized": hex(serde_json::to_string(&es).unwrap().as_bytes()),
        "vertices": r.ids.len(), "raw_edges": r.rows.len(), "dup_ids_in_edges": r.dup,
        "empty_edge": r.rows.iter().any(Vec::is_empty), "normalized_edges": h.ev.len(),
        "duplicate_rows_merged": merged, "superset_rows_dropped": r.rows.len() - merged - h.ev.len(),
    })
}

pub fn status(st: &St, k: Option<usize>) -> &'static str {
    match (st.done, k) {
        (false, _) => "resource_limited",
        (true, Some(_)) => "bounded_complete",
        (true, None) => "complete",
    }
}

pub fn ord(x: Option<&str>) -> Ord {
    match x {
        None | Some("idx") => Ord::Idx,
        Some("unhit") => Ord::Unhit,
        Some(x) => die(&format!("unknown --ord {x}")),
    }
}

/// Deadline secs after t; secs must be finite and nonnegative.
pub fn deadline(t: Instant, secs: f64) -> Instant {
    if !(secs.is_finite() && secs >= 0.0) {
        die(&format!("time limit {secs}: expected nonnegative seconds"));
    }
    Duration::try_from_secs_f64(secs)
        .ok()
        .and_then(|d| t.checked_add(d))
        .unwrap_or_else(|| die("time limit exceeds the supported clock range"))
}

/// Expanded size of a class-level solution.
pub fn size(cls: &[Vec<usize>], rep: &[usize]) -> Option<u128> {
    rep.iter().try_fold(1u128, |n, &c| n.checked_mul(cls[c].len() as u128))
}

pub fn enumerate(a: &[String]) {
    let t0 = Instant::now();
    let x = Args::new(
        a,
        1,
        &["--ord", "--k", "--time", "--out", "--summary"],
        &["--twins", "--count", "--census", "--check"],
    );
    let p = &x.pos[0];
    let ord = ord(x.get("--ord"));
    let k: Option<usize> = x.num("--k");
    let secs: Option<f64> = x.num("--time");
    let (tws, cnt, cen, chk_on) = (x.has("--twins"), x.has("--count"), x.has("--census"), x.has("--check"));
    if cnt && x.get("--out").is_some() {
        die("--count and --out are exclusive");
    }
    let (b, r) = inp::load(p);
    let (h, merged) = hg::norm(r.ids.len(), &r.rows);
    let tw = hg::twins(&h, tws);
    let t1 = Instant::now();
    let cfg = Cfg { ord, k, census: cen, dl: secs.map(|s| deadline(t1, s)) };
    let pending = |p| Pending::new(p).unwrap_or_else(|e| die(&format!("{p}: {e}")));
    let mut output = x.get("--out").map(pending);
    let mut summary = x.get("--summary").map(pending);
    if output.as_ref().zip(summary.as_ref()).is_some_and(|(a, b)| a.dest == b.dest) {
        die("solutions and summary need different filenames");
    }
    let mut w = output.as_mut().map(|o| BufWriter::new(&mut o.file));
    let (mut nexp, mut bad, mut dups) = (Some(0u128), 0u64, 0u64);
    let mut first_bad: Option<Value> = None;
    let mut seen: HashSet<Vec<usize>> = HashSet::new();
    let mut t_first: Option<f64> = None;
    let cls = &tw.cls;
    let st = srch::run(tw.h.n, &tw.h.ev, &cfg, &mut |rep| {
        t_first.get_or_insert_with(|| t1.elapsed().as_secs_f64());
        if cnt && !chk_on {
            nexp = nexp.and_then(|n| n.checked_add(size(cls, rep)?));
            return;
        }
        hg::expand(cls, rep, &mut |c| {
            nexp = nexp.and_then(|n| n.checked_add(1));
            if chk_on {
                if let Err(e) = chk::check(&r.rows, c) {
                    bad += 1;
                    first_bad.get_or_insert(json!({"solution": c.iter().map(|&v| &r.ids[v]).collect::<Vec<_>>(), "error": format!("{e:?}")}));
                }
                dups += !seen.insert(c.to_vec()) as u64;
            }
            if let Some(w) = w.as_mut() {
                let ids: Vec<&str> = c.iter().map(|&v| r.ids[v].as_str()).collect();
                writeln!(w, "{}", serde_json::to_string(&ids).unwrap()).unwrap_or_else(|e| die(&format!("write: {e}")));
            }
        })
    });
    if let Some(mut w) = w.take() {
        w.flush().unwrap_or_else(|e| die(&format!("write: {e}")));
    }
    drop(w);
    let t2 = Instant::now();
    let verification_status = match (chk_on, bad + dups) {
        (false, _) => "not_run",
        (true, 0) => "passed",
        (true, _) => "failed",
    };
    let c = &st.cen;
    let n = st.nodes.max(1) as f64;
    let census = cen.then(|| {
        json!({
            "nodes": st.nodes, "dead_nodes": c.dead, "dead_share": c.dead as f64 / n,
            "exact": {"distinct_keys": c.ex_keys, "repeat_roots": c.ex_roots, "skippable_nodes": c.ex_nodes,
                      "skippable_share": c.ex_nodes as f64 / n, "skippable_dead_nodes": c.ex_dead},
            "wl_upper_bound": {"distinct_keys": c.wl_keys, "repeat_roots": c.wl_roots, "skippable_nodes": c.wl_nodes,
                               "skippable_share": c.wl_nodes as f64 / n, "skippable_dead_nodes": c.wl_dead},
            "note": "exact keys are 128-bit hashes of the restricted residual state; WL equality is necessary for isomorphism, so its share bounds isomorphic reuse from above",
        })
    });
    let receipt = output.as_mut().map(|o| {
        let hash = o.digest().unwrap_or_else(|e| die(&format!("solutions: {e}")));
        json!({"path": o.dest, "sha256": hash})
    });
    let sm = json!({
        "tool": "mcs", "version": env!("CARGO_PKG_VERSION"),
        "input": info(p, &b, &r, &h, merged),
        "method": {"ord": format!("{ord:?}").to_lowercase(), "twins": tws, "k": k, "count_only": cnt,
                   "census": cen, "check": chk_on, "time_limit_s": secs},
        "reduction": {"classes": cls.len(), "nontrivial_classes": cls.iter().filter(|c| c.len() > 1).count()},
        "search_status": status(&st, k),
        "verification_status": verification_status,
        "verification": chk_on.then(|| json!({"failed": bad, "duplicates": dups, "first_failure": first_bad})),
        "counts": {"representative_solutions": st.sols, "expanded_solutions": nexp.map(|n| n.to_string()), "nodes": st.nodes},
        "count_status": if nexp.is_some() {"exact"} else {"overflow"},
        "output": receipt,
        "elapsed_s": {"load_normalize": (t1 - t0).as_secs_f64(), "search": (t2 - t1).as_secs_f64(),
                      "first_solution": t_first, "total": t0.elapsed().as_secs_f64()},
        "census": census,
    });
    let s = serde_json::to_string_pretty(&sm).unwrap();
    if let Some(o) = summary.as_mut() {
        o.file.write_all(s.as_bytes()).unwrap_or_else(|e| die(&format!("summary: {e}")));
        o.file.sync_all().unwrap_or_else(|e| die(&format!("summary: {e}")));
    }
    if let Some(o) = output {
        o.publish().unwrap_or_else(|e| die(&format!("solutions: {e}")));
    }
    if let Some(o) = summary {
        o.publish().unwrap_or_else(|e| die(&format!("summary: {e}")));
    } else {
        println!("{s}");
    }
    if verification_status == "failed" {
        std::process::exit(1);
    }
}

pub fn check(a: &[String]) {
    let x = Args::new(a, 2, &["--ref", "--summary"], &[]);
    let (bytes, r) = inp::load(&x.pos[0]);
    let ix: HashMap<&str, usize> = r.ids.iter().enumerate().map(|(i, s)| (s.as_str(), i)).collect();
    let sols = inp::sols(&x.pos[1], &ix);
    let fails: Vec<(&Vec<usize>, chk::Bad)> =
        sols.iter().filter_map(|c| chk::check(&r.rows, c).err().map(|e| (c, e))).collect();
    let first: Vec<Value> = fails
        .iter()
        .take(10)
        .map(|(c, e)| json!({"solution": c.iter().map(|&v| &r.ids[v]).collect::<Vec<_>>(), "error": format!("{e:?}")}))
        .collect();
    let duplicates = sols.len() - sols.iter().collect::<HashSet<_>>().len();
    let mut out =
        json!({"solutions": sols.len(), "failed": fails.len(), "duplicates": duplicates, "first_failures": first});
    let mut ok = fails.is_empty() && duplicates == 0;
    if let Some(path) = x.get("--summary") {
        let read = |p| std::fs::read(p).unwrap_or_else(|e| die(&format!("{p}: {e}")));
        let summary: Value = serde_json::from_slice(&read(path)).unwrap_or_else(|e| die(&format!("{path}: {e}")));
        let matched = summary["input"]["sha256_bytes"] == hex(&bytes)
            && summary["output"]["sha256"] == hex(&read(&x.pos[1]))
            && summary["counts"]["expanded_solutions"].as_str().and_then(|n| n.parse::<u128>().ok())
                == Some(sols.len() as u128)
            && summary["count_status"] == "exact"
            && summary["verification_status"] != "failed"
            && matches!(summary["search_status"].as_str(), Some("complete" | "bounded_complete" | "resource_limited"));
        ok &= matched;
        out["summary"] = json!({"path": path, "matches": matched, "search_status": summary["search_status"]});
    }
    if let Some(rf) = x.get("--ref") {
        let refs = inp::sols(rf, &ix);
        let c = chk::cmp(&sols, &refs);
        let exact = c.missing == 0 && c.extra == 0 && c.dup_a == 0;
        ok &= exact;
        out["reference"] = json!({"path": rf, "solutions": refs.len(), "missing": c.missing, "extra": c.extra,
                                  "duplicates_tested": c.dup_a, "duplicates_reference": c.dup_b, "exact_family_match": exact});
    }
    out["verification_status"] = json!(if ok { "passed" } else { "failed" });
    println!("{}", serde_json::to_string_pretty(&out).unwrap());
    if !ok {
        std::process::exit(1);
    }
}
