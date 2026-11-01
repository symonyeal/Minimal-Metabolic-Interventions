// modes and prepare: elementary modes of a BiGG model, and every input used in EVIDENCE.md.
//
// Legend:
//   A     lineage A: METATOOL E. coli elementary-mode supports from MHSGenerationAlgorithms @ REV
//   B     lineage B: BiGG e_coli_core, its elementary modes and three cut-set scenarios
//   grow  A growth instance: modes containing mue (id 1), with mue removed
//
// Downloads go through the system curl and are checked against pinned SHA-256 digests; a
// mismatch stops the run. data/ is not tracked: third-party redistribution terms are unresolved.

use crate::inp::{self, SCHEMA};
use crate::{Args, die, hex};
use mcs::net::{self, Net, Pol};
use serde_json::{Value, json};
use std::path::Path;
use std::process::Command;
use std::time::Instant;

const REPO: &str = "VeraLiconaResearchGroup/MHSGenerationAlgorithms";
const REV: &str = "f043cc3fae8088191307d59b987ce5fe9de4a326";
const A: [(&str, &str); 5] = [
    ("acetate", "bf5f50c582c8631af024dfcfe574f6892721bb72887c039e275f0ebf0b9b5717"),
    ("succinate", "8ad10e583218aa19b91530c6eb26d4f417b1fa6e7f9b6f3f2b639511ecea44fc"),
    ("glycerol", "a0b1afc2eec0852df9fbd5b98a57b2bbef25d78bbc00ae50a31d080ddadac413"),
    ("glucose", "a666e621bd597c6a8ce97931c0ebb49a9464063459267d0741be436d2284ef89"),
    ("combined", "141f94db708f88ff12781f5aac213ffe9aa7718c3a84dd03179428f7194a27db"),
];
const B_URL: &str = "https://bigg.ucsd.edu/static/models/e_coli_core.json";
const B_SHA: &str = "7bedec10576cfe935b19218dc881f3fb14f890a1871448fc19a9b4ee15b448d8";
const B_SCEN: [(Option<&str>, Pol, &str); 3] = [
    (Some("BIOMASS_Ecoli_core_w_GAM"), Pol::AllButTarget, "growth"),
    (Some("EX_etoh_e"), Pol::Core, "ethanol"),
    (None, Pol::AllButTarget, "all"),
];

fn write(p: &Path, b: &[u8]) {
    std::fs::write(p, b).unwrap_or_else(|e| die(&format!("{}: {e}", p.display())));
}

fn read(p: &Path) -> Vec<u8> {
    std::fs::read(p).unwrap_or_else(|e| die(&format!("{}: {e}", p.display())))
}

fn name(p: &Path) -> String {
    p.file_name().unwrap().to_string_lossy().into_owned()
}

/// Downloads u to out unless out already has digest sha.
fn get(u: &str, sha: &str, out: &Path) {
    if out.exists() && hex(&read(out)) == sha {
        return;
    }
    let tmp = out.with_extension("part");
    let ok =
        Command::new("curl").args(["-fsSL", "--retry", "3", "-o"]).arg(&tmp).arg(u).status().is_ok_and(|s| s.success());
    if !ok {
        die(&format!("{u}: download failed (curl)"));
    }
    let h = hex(&read(&tmp));
    if h != sha {
        die(&format!("{u}: sha256 {h}, expected {sha}"));
    }
    std::fs::rename(&tmp, out).unwrap_or_else(|e| die(&format!("{}: {e}", out.display())));
}

/// Writes the modes of the model as JSON lines and returns the network, the reaction supports and
/// the summary.
fn efm(model: &Path, out: &Path) -> (Net, Vec<Vec<usize>>, Value) {
    let t0 = Instant::now();
    let b = read(model);
    let net = inp::bigg(&b).unwrap_or_else(|e| die(&format!("{}: {e}", model.display())));
    let cone = net::cone(&net).unwrap_or_else(|e| die(&e));
    let dd = mcs::dd::rays(&cone.a, cone.cols.len()).unwrap_or_else(|e| die(&e));
    let m = net::modes(&net, &cone, dd.rays);
    let id = |j: usize| net.rx[j].id.as_str();
    let mut o = vec![];
    for (cols, sup) in m.cols.iter().zip(&m.sup) {
        let mut ids: Vec<&str> = sup.iter().map(|&j| id(j)).collect();
        ids.sort_unstable();
        let dirs: Vec<Value> = cols.iter().map(|&k| json!([id(cone.cols[k].0), cone.cols[k].1])).collect();
        o.push(json!({"support": ids, "columns": dirs}).to_string());
    }
    let o = o.join("\n") + "\n";
    write(out, o.as_bytes());
    let sm = json!({
        "tool": "mcs modes", "version": env!("CARGO_PKG_VERSION"),
        "model": name(model), "model_sha256": hex(&b), "modes_sha256": hex(o.as_bytes()),
        "reactions": net.rx.len(), "columns": cone.cols.len(), "rank": dd.rank, "free": dd.free,
        "modes": m.sup.len(), "two_cycles": m.two_cycles, "peak_rays": dd.peak,
        "pruned": rids(&net, &cone.pruned), "blocked": rids(&net, &m.blocked),
        "elapsed_s": t0.elapsed().as_secs_f64(),
    });
    (net, m.sup, sm)
}

fn rids<'a>(net: &'a Net, js: &[usize]) -> Vec<&'a str> {
    js.iter().map(|&j| net.rx[j].id.as_str()).collect()
}

pub fn modes(a: &[String]) {
    let x = Args::new(a, 2, &[], &[]);
    let (_, _, sm) = efm(Path::new(&x.pos[0]), Path::new(&x.pos[1]));
    println!("{}", serde_json::to_string_pretty(&sm).unwrap());
}

pub fn prepare(a: &[String]) {
    let x = Args::new(a, 1, &[], &["--skip-b"]);
    let dat = Path::new(&x.pos[0]);
    let (da, db) = (dat.join("a"), dat.join("b"));
    std::fs::create_dir_all(&da).unwrap_or_else(|e| die(&format!("{}: {e}", da.display())));
    for (n, sha) in A {
        // The JSON files are stored with Git LFS; media.githubusercontent.com serves their content.
        let p = da.join(format!("ecoli-{n}.json"));
        get(
            &format!("https://media.githubusercontent.com/media/{REPO}/{REV}/benchmark_inputs/ecoli/ecoli-{n}.json"),
            sha,
            &p,
        );
        let v: Value = serde_json::from_slice(&read(&p)).unwrap_or_else(|e| die(&format!("{}: {e}", p.display())));
        let s = v["sets"].as_array().unwrap_or_else(|| die("sets: expected an array"));
        let has1 = |e: &&Value| e.as_array().is_some_and(|e| e.iter().any(|x| x.as_u64() == Some(1)));
        let g: Vec<Vec<&Value>> = s
            .iter()
            .filter(has1)
            .map(|e| e.as_array().unwrap().iter().filter(|x| x.as_u64() != Some(1)).collect())
            .collect();
        let doc = json!({"sets": g, "scenario": "growth blocking: modes containing mue (id 1), mue removed", "source": format!("{REPO}@{REV}")});
        write(&da.join(format!("ecoli-{n}-growth.json")), doc.to_string().as_bytes());
        println!("{n}: {} modes, {} growth modes", s.len(), g.len());
    }
    if x.has("--skip-b") {
        return;
    }
    std::fs::create_dir_all(&db).unwrap_or_else(|e| die(&format!("{}: {e}", db.display())));
    let mdl = db.join("e_coli_core.json");
    get(B_URL, B_SHA, &mdl);
    let mo = db.join("modes.jsonl");
    let (net, sup, sm) = efm(&mdl, &mo);
    write(&db.join("modes-summary.json"), serde_json::to_string_pretty(&sm).unwrap().as_bytes());
    println!("modes: {}", sm);
    let id = |j: usize| net.rx[j].id.as_str();
    for (t, pol, nm) in B_SCEN {
        let tj = t
            .map(|t| net.rx.iter().position(|r| r.id == t).unwrap_or_else(|| die(&format!("target {t} not in model"))));
        let (es, mut vs, el) = net::scenario(&net, &sup, tj, pol);
        vs.sort_by_key(|&j| id(j));
        let edges: Vec<Vec<&str>> = es
            .iter()
            .map(|s| {
                let mut e: Vec<&str> = s.iter().map(|&j| id(j)).collect();
                e.sort_unstable();
                e
            })
            .collect();
        let doc = json!({
            "schema": SCHEMA,
            "vertices": vs.iter().map(|&j| json!({"id": id(j), "label": net.rx[j].name})).collect::<Vec<_>>(),
            "edges": edges,
            "eligible": el.iter().map(|&j| id(j)).collect::<Vec<_>>(),
            "meta": {"model": name(&mdl), "model_sha256": sm["model_sha256"], "modes": name(&mo), "modes_sha256": sm["modes_sha256"],
                     "target": t.unwrap_or("all"), "policy": if pol == Pol::Core { "core" } else { "all-but-target" }, "target_modes": es.len(), "all_modes": sup.len()},
        });
        write(&db.join(format!("ecc-{nm}.json")), doc.to_string().as_bytes());
        println!("{nm}: {} edges, {} vertices, {} eligible", es.len(), vs.len(), el.len());
    }
}
