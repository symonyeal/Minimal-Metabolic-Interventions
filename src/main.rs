// mcs: enumerate inclusion-minimal hitting sets (minimal intervention sets) of an explicit
// hypergraph.
//
//   mcs inspect FILE
//   mcs enumerate FILE [--alg mmcs|rs] [--ord idx|unhit] [--twins] [--k K] [--count]
//                 [--out SOLS] [--summary JSON] [--census] [--time SECS] [--cap N] [--check]
//   mcs check FILE SOLS [--ref SOLS]
//
// Legend:
//   r    parsed input (hg::Raw)      h   normalized hypergraph     tw  twin reduction
//   rep  class-level solutions        exp expanded (labelled) solutions
//   sm   summary JSON                 t0  start instant
//
// The summary is written last, through a temporary file and rename: a missing summary means the
// run did not finish. A solutions file is valid only with its summary.

use mcs::chk;
use mcs::hg::{self, Raw};
use mcs::srch::{self, Alg, Cfg, Ord, Stop};
use serde_json::{Value, json};
use sha2::{Digest, Sha256};
use std::collections::{HashMap, HashSet};
use std::io::{BufWriter, Write};
use std::time::{Duration, Instant};

fn die(s: &str) -> ! {
    eprintln!("mcs: {s}");
    std::process::exit(2)
}

/// Peak working set of this process in bytes (Windows); None elsewhere until implemented.
#[cfg(windows)]
fn peak_mem() -> Option<u64> {
    #[repr(C)]
    #[derive(Default)]
    struct Pmc {
        cb: u32,
        faults: u32,
        peak_ws: usize,
        rest: [usize; 7],
    }
    unsafe extern "system" {
        fn GetCurrentProcess() -> isize;
        fn K32GetProcessMemoryInfo(h: isize, p: *mut Pmc, cb: u32) -> i32;
    }
    let mut p = Pmc { cb: size_of::<Pmc>() as u32, ..Default::default() };
    // SAFETY: p is a correctly sized PROCESS_MEMORY_COUNTERS; the pseudo-handle needs no close.
    let ok = unsafe { K32GetProcessMemoryInfo(GetCurrentProcess(), &mut p, p.cb) };
    (ok != 0).then_some(p.peak_ws as u64)
}
#[cfg(not(windows))]
fn peak_mem() -> Option<u64> {
    None
}

fn hex(b: &[u8]) -> String {
    Sha256::digest(b).iter().map(|x| format!("{x:02x}")).collect()
}

fn opt(a: &[String], k: &str) -> Option<String> {
    a.iter().position(|x| x == k).map(|i| a.get(i + 1).cloned().unwrap_or_else(|| die(&format!("{k} needs a value"))))
}

fn flag(a: &[String], k: &str) -> bool {
    a.iter().any(|x| x == k)
}

fn load(p: &str) -> (Vec<u8>, Raw) {
    let b = std::fs::read(p).unwrap_or_else(|e| die(&format!("{p}: {e}")));
    match hg::parse(&b) {
        Ok(r) => (b, r),
        Err(e) => {
            println!("{}", json!({"search_status": "invalid_input", "error": e}));
            std::process::exit(3)
        }
    }
}

/// Canonical normalized form: sorted original ids per edge, edges in normalized order.
fn norm_hash(r: &Raw, h: &hg::H) -> String {
    let es: Vec<Vec<&str>> = h.ev.iter().map(|b| b.iter().map(|v| r.ids[v].as_str()).collect()).collect();
    hex(serde_json::to_string(&es).unwrap().as_bytes())
}

fn input_info(p: &str, b: &[u8], r: &Raw, h: &hg::H) -> Value {
    let merged: usize = h.src.iter().map(|s| s.len() - 1).sum();
    json!({
        "path": p, "format": r.fmt, "sha256_bytes": hex(b), "sha256_normalized": norm_hash(r, h),
        "vertices": r.ids.len(), "raw_edges": r.rows.len(), "dup_ids_in_edges": r.dup,
        "empty_edge": h.empty, "normalized_edges": h.ev.len(),
        "duplicate_rows_merged": merged, "superset_rows_dropped": h.dom.len(),
    })
}

/// Solutions file: a JSON array of arrays, or JSON lines of arrays; ids as strings or integers.
fn read_sols(p: &str, ix: &HashMap<&str, usize>) -> Vec<Vec<usize>> {
    let t = std::fs::read_to_string(p).unwrap_or_else(|e| die(&format!("{p}: {e}")));
    let rows: Vec<Value> = match serde_json::from_str::<Value>(&t) {
        Ok(Value::Array(a)) if a.iter().all(Value::is_array) => a,
        _ => t.lines().filter(|l| !l.trim().is_empty()).map(|l| serde_json::from_str(l).unwrap_or_else(|e| die(&format!("{p}: {e}")))).collect(),
    };
    rows.iter()
        .map(|r| {
            let mut c: Vec<usize> = r
                .as_array()
                .unwrap_or_else(|| die("solution: expected an array"))
                .iter()
                .map(|x| {
                    let s = match x {
                        Value::String(s) => s.clone(),
                        Value::Number(n) => n.to_string(),
                        _ => die("solution id: expected a string or integer"),
                    };
                    *ix.get(s.as_str()).unwrap_or_else(|| die(&format!("solution uses unknown id {s}")))
                })
                .collect();
            c.sort_unstable();
            c
        })
        .collect()
}

fn main() {
    let a: Vec<String> = std::env::args().skip(1).collect();
    match a.first().map(String::as_str) {
        Some("inspect") => inspect(&a[1..]),
        Some("enumerate") => enumerate(&a[1..]),
        Some("check") => check(&a[1..]),
        _ => die("usage: mcs inspect|enumerate|check FILE ..."),
    }
}

fn inspect(a: &[String]) {
    let p = a.first().unwrap_or_else(|| die("inspect needs FILE"));
    let (b, r) = load(p);
    let h = hg::norm(&r);
    let tw = hg::twins(&h, true);
    let big: Vec<usize> = tw.cls.iter().map(Vec::len).filter(|&l| l > 1).collect();
    let inc = tw.cls.iter().map(Vec::len).sum::<usize>();
    println!(
        "{}",
        serde_json::to_string_pretty(&json!({
            "input": input_info(p, &b, &r, &h),
            "vertices_in_edges": inc, "twin_classes": tw.cls.len(), "nontrivial_class_sizes": big,
        }))
        .unwrap()
    );
}

fn enumerate(a: &[String]) {
    let t0 = Instant::now();
    let p = a.first().unwrap_or_else(|| die("enumerate needs FILE"));
    let alg = match opt(a, "--alg").as_deref() {
        None | Some("mmcs") => Alg::Mmcs,
        Some("rs") => Alg::Rs,
        Some(x) => die(&format!("unknown --alg {x}")),
    };
    let ord = match opt(a, "--ord").as_deref() {
        None | Some("idx") => Ord::Idx,
        Some("unhit") => Ord::Unhit,
        Some(x) => die(&format!("unknown --ord {x}")),
    };
    let k = opt(a, "--k").map(|x| x.parse::<usize>().unwrap_or_else(|_| die("--k: expected an integer")));
    let cap = opt(a, "--cap").map(|x| x.parse::<u64>().unwrap_or_else(|_| die("--cap: expected an integer")));
    let secs = opt(a, "--time").map(|x| x.parse::<f64>().unwrap_or_else(|_| die("--time: expected seconds")));
    let (tws, cnt, cen, chk_on) = (flag(a, "--twins"), flag(a, "--count"), flag(a, "--census"), flag(a, "--check"));
    if cen && alg != Alg::Mmcs {
        die("--census is defined for MMCS states only");
    }
    if cnt && a.iter().any(|x| x == "--out") {
        die("--count and --out are exclusive");
    }
    let (b, r) = load(p);
    let h = hg::norm(&r);
    let tw = hg::twins(&h, tws);
    let t1 = Instant::now();

    let mut w = opt(a, "--out").map(|o| BufWriter::new(std::fs::File::create(&o).unwrap_or_else(|e| die(&format!("{o}: {e}")))));
    let (mut nexp, mut bad, mut dups) = (0u128, 0u64, 0u64);
    let mut first_bad: Option<Value> = None;
    let mut seen: HashSet<Vec<usize>> = HashSet::new();
    let mut t_first: Option<f64> = None;
    let cls = &tw.cls;
    let cfg = Cfg { alg, ord, k, census: cen, cap, dl: secs.map(|s| t1 + Duration::from_secs_f64(s)) };
    let mut emit = |rep: &[usize]| -> bool {
        if t_first.is_none() {
            t_first = Some(t1.elapsed().as_secs_f64());
        }
        if cnt && !chk_on {
            nexp += rep.iter().map(|&c| cls[c].len() as u128).product::<u128>();
            return true;
        }
        hg::expand(cls, rep, &mut |c| {
            nexp += 1;
            if chk_on {
                if let Err(e) = chk::check(&r.rows, c) {
                    bad += 1;
                    first_bad.get_or_insert(json!({"solution": c.iter().map(|&v| &r.ids[v]).collect::<Vec<_>>(), "error": format!("{e:?}")}));
                }
                if !seen.insert(c.to_vec()) {
                    dups += 1;
                }
            }
            if let Some(w) = w.as_mut() {
                let ids: Vec<&str> = c.iter().map(|&v| r.ids[v].as_str()).collect();
                writeln!(w, "{}", serde_json::to_string(&ids).unwrap()).unwrap_or_else(|e| die(&format!("write: {e}")));
            }
            true
        })
    };
    let st = srch::run(tw.h.n, &tw.h.ev, &cfg, &mut emit);
    if let Some(mut w) = w {
        w.flush().unwrap_or_else(|e| die(&format!("write: {e}")));
    }
    let t2 = Instant::now();
    let search_status = match st.stop {
        Stop::Done if k.is_some() => "bounded_complete",
        Stop::Done => "complete",
        Stop::Cap | Stop::Time => "resource_limited",
    };
    let verification_status = if !chk_on {
        "not_run"
    } else if bad == 0 && dups == 0 {
        "passed"
    } else {
        "failed"
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
    let sm = json!({
        "tool": "mcs", "version": env!("CARGO_PKG_VERSION"),
        "input": input_info(p, &b, &r, &h),
        "method": {"alg": format!("{alg:?}").to_lowercase(), "ord": format!("{ord:?}").to_lowercase(), "twins": tws,
                   "k": k, "count_only": cnt, "census": cen, "check": chk_on, "cap": cap, "time_limit_s": secs},
        "reduction": {"classes": tw.cls.len(), "nontrivial_classes": tw.cls.iter().filter(|c| c.len() > 1).count()},
        "search_status": search_status,
        "limit": match st.stop { Stop::Cap => Some("output_cap"), Stop::Time => Some("time"), Stop::Done => None },
        "verification_status": verification_status,
        "verification": chk_on.then(|| json!({"failed": bad, "duplicates": dups, "first_failure": first_bad})),
        "counts": {"representative_solutions": st.sols, "expanded_solutions": nexp.to_string(), "nodes": st.nodes},
        "elapsed_s": {"load_normalize": (t1 - t0).as_secs_f64(), "search": (t2 - t1).as_secs_f64(),
                      "first_solution": t_first, "total": t0.elapsed().as_secs_f64()},
        "census": census,
        "peak_memory": {"bytes": peak_mem(), "method": if cfg!(windows) { "peak working set, K32GetProcessMemoryInfo" } else { "not measured" }},
    });
    let s = serde_json::to_string_pretty(&sm).unwrap();
    match opt(a, "--summary") {
        Some(o) => {
            let tmp = format!("{o}.tmp");
            std::fs::write(&tmp, &s).and_then(|_| std::fs::rename(&tmp, &o)).unwrap_or_else(|e| die(&format!("{o}: {e}")));
        }
        None => println!("{s}"),
    }
    if verification_status == "failed" {
        std::process::exit(1);
    }
}

fn check(a: &[String]) {
    let (p, q) = match (a.first(), a.get(1)) {
        (Some(p), Some(q)) => (p, q),
        _ => die("check needs FILE SOLS"),
    };
    let (_, r) = load(p);
    let ix: HashMap<&str, usize> = r.ids.iter().enumerate().map(|(i, s)| (s.as_str(), i)).collect();
    let sols = read_sols(q, &ix);
    let mut bad: Vec<Value> = vec![];
    for c in &sols {
        if let Err(e) = chk::check(&r.rows, c) {
            if bad.len() < 10 {
                bad.push(json!({"solution": c.iter().map(|&v| &r.ids[v]).collect::<Vec<_>>(), "error": format!("{e:?}")}));
            }
        }
    }
    let nbad = sols.iter().filter(|c| chk::check(&r.rows, c).is_err()).count();
    let mut out = json!({"solutions": sols.len(), "failed": nbad, "first_failures": bad});
    let mut ok = nbad == 0;
    if let Some(rf) = opt(a, "--ref") {
        let refs = read_sols(&rf, &ix);
        let c = chk::cmp(&sols, &refs);
        ok &= c.missing == 0 && c.extra == 0 && c.dup_a == 0;
        out["reference"] = json!({"path": rf, "solutions": refs.len(), "missing": c.missing, "extra": c.extra,
                                  "duplicates_tested": c.dup_a, "duplicates_reference": c.dup_b, "exact_family_match": c.missing == 0 && c.extra == 0 && c.dup_a == 0});
    }
    out["verification_status"] = json!(if ok { "passed" } else { "failed" });
    println!("{}", serde_json::to_string_pretty(&out).unwrap());
    if !ok {
        std::process::exit(1);
    }
}
