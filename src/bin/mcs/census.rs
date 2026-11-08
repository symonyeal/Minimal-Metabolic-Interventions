// census: dead and repeated MMCS subtrees per instance, for both candidate orders with and
// without static twins. Writes OUT/census.json and OUT/census.md; shares are deterministic,
// search times are not.
//
// Legend:
//   dead  share of nodes whose subtree emits nothing
//   ex    share of nodes inside subtrees rooted at an exact repeat of a completed state
//   wl    same with colour-refinement equality: an upper bound on isomorphic repeats
//   gap   wl − ex: the most isomorphism can add over exact memoization

use crate::enumerate::{deadline, size, status};
use crate::output::Pending;
use crate::{Args, die, inp};
use mcs::hg;
use mcs::srch::{self, Cfg, Ord};
use serde_json::{Value, json};
use std::io::Write;
use std::path::Path;
use std::time::Instant;

fn pct(x: f64) -> String {
    format!("{:.1}%", 100.0 * x)
}

/// Integer with comma thousands separators.
fn sep(x: u64) -> String {
    let s = x.to_string();
    let mut o = String::new();
    for (i, ch) in s.chars().enumerate() {
        if i > 0 && (s.len() - i).is_multiple_of(3) {
            o.push(',');
        }
        o.push(ch);
    }
    o
}

pub fn census(a: &[String]) {
    let x = Args::new(a, 3, &[], &[]);
    let out = Path::new(&x.pos[0]);
    let secs: f64 = x.pos[1].parse().unwrap_or_else(|_| die("census: SECS must be a number"));
    deadline(Instant::now(), secs);
    let inputs: Vec<_> = x.pos[2..].iter().map(|p| (p, inp::load(p).1)).collect(); // before staging files
    std::fs::create_dir_all(out).unwrap_or_else(|e| die(&format!("{}: {e}", out.display())));
    let pending = |name| {
        let p = out.join(name);
        Pending::new(p.to_str().unwrap()).unwrap_or_else(|e| die(&format!("{}: {e}", p.display())))
    };
    let (mut data, mut table) = (pending("census.json"), pending("census.md"));
    let mut rows: Vec<Value> = vec![];
    for (p, r) in &inputs {
        let (h, _) = hg::norm(r.ids.len(), &r.rows);
        let stem = Path::new(p).file_stem().unwrap().to_string_lossy();
        for ord in [Ord::Idx, Ord::Unhit] {
            for tw in [false, true] {
                let t = hg::twins(&h, tw);
                let t1 = Instant::now();
                let cfg = Cfg { ord, k: None, census: true, dl: Some(deadline(t1, secs)) };
                let mut nexp = Some(0u128);
                let st = srch::run(t.h.n, &t.h.ev, &cfg, &mut |rep| {
                    nexp = nexp.and_then(|n| n.checked_add(size(&t.cls, rep)?));
                });
                let (c, n) = (&st.cen, st.nodes.max(1) as f64);
                rows.push(json!({
                    "instance": stem, "ord": format!("{ord:?}").to_lowercase(), "twins": tw, "status": status(&st, None),
                    "nodes": st.nodes, "solutions": nexp.map(|n| n.to_string()),
                    "count_status": if nexp.is_some() {"exact"} else {"overflow"},
                    "dead": c.dead as f64 / n, "ex": c.ex_nodes as f64 / n,
                    "ex_dead_nodes": c.ex_dead, "wl": c.wl_nodes as f64 / n, "search_s": t1.elapsed().as_secs_f64(),
                }));
                println!("{}", rows.last().unwrap());
            }
        }
    }
    let mut md = vec![
        "| instance | order | twins | status | nodes | dead | exact | WL bound | gap |".to_owned(),
        "|---|---|---|---|---:|---:|---:|---:|---:|".to_owned(),
    ];
    for r in &rows {
        let f = |k: &str| r[k].as_f64().unwrap();
        md.push(format!(
            "| {} | {} | {} | {} | {} | {} | {} | {} | {} |",
            r["instance"].as_str().unwrap(),
            r["ord"].as_str().unwrap(),
            if r["twins"].as_bool().unwrap() { "yes" } else { "no" },
            r["status"].as_str().unwrap(),
            sep(r["nodes"].as_u64().unwrap()),
            pct(f("dead")),
            pct(f("ex")),
            pct(f("wl")),
            pct(f("wl") - f("ex")),
        ));
    }
    data.file
        .write_all((serde_json::to_string_pretty(&rows).unwrap() + "\n").as_bytes())
        .and_then(|_| table.file.write_all((md.join("\n") + "\n").as_bytes()))
        .and_then(|_| data.file.sync_all())
        .and_then(|_| table.file.sync_all())
        .unwrap_or_else(|e| die(&format!("census: {e}")));
    table.publish().and_then(|_| data.publish()).unwrap_or_else(|e| die(&format!("census: {e}")));
}
