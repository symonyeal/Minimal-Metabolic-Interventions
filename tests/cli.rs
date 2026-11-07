// The mcs command line against the tracked oracle families in results/oracle/.

use std::path::{Path, PathBuf};
use std::process::Command;

fn mcs(a: &[&str]) -> (i32, String) {
    let o = Command::new(env!("CARGO_BIN_EXE_mcs")).args(a).current_dir(env!("CARGO_MANIFEST_DIR")).output().unwrap();
    (o.status.code().unwrap(), String::from_utf8(o.stdout).unwrap())
}

fn tmp(f: &str) -> PathBuf {
    static NEXT: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
    let root = std::env::var_os("MCS_TEST_ROOT")
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from(env!("CARGO_TARGET_TMPDIR")));
    let dir = root.join(format!(
        "{}-{}-{}",
        std::process::id(),
        std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_nanos(),
        NEXT.fetch_add(1, std::sync::atomic::Ordering::Relaxed)
    ));
    std::fs::create_dir_all(&dir).unwrap();
    dir.join(f)
}

fn read(p: &Path) -> serde_json::Value {
    serde_json::from_slice(&std::fs::read(p).unwrap()).unwrap()
}

fn pairs(n: usize, triangle: bool) -> PathBuf {
    use serde_json::json;
    let extra = if triangle { 6 } else { 0 };
    let ids: Vec<_> = (0..2 * n + extra).map(|i| i.to_string()).collect();
    let mut edges: Vec<Vec<_>> = (0..n).map(|i| vec![ids[2 * i].clone(), ids[2 * i + 1].clone()]).collect();
    if triangle {
        // Three distinct twin classes, any two of which hit the triangle's three edges.
        for (a, b) in [(0, 1), (0, 2), (1, 2)] {
            edges.push([2 * a, 2 * a + 1, 2 * b, 2 * b + 1].iter().map(|i| ids[2 * n + i].clone()).collect());
        }
    }
    let p = tmp("pairs.json");
    std::fs::write(&p, json!({"schema":"mcs-hg/1", "vertices":ids.iter().map(|id|json!({"id":id})).collect::<Vec<_>>(), "edges":edges}).to_string()).unwrap();
    p
}

#[test]
fn expanded_count_never_wraps_products_or_sums() {
    for (n, triangle, exact) in
        [(127, false, Some(1u128 << 127)), (128, false, None), (129, false, None), (125, true, None)]
    {
        let p = pairs(n, triangle);
        let (code, text) = mcs(&["enumerate", p.to_str().unwrap(), "--twins", "--count"]);
        assert_eq!(code, 0, "{text}");
        let value: serde_json::Value = serde_json::from_str(&text).unwrap();
        assert_eq!(value["search_status"], "complete");
        assert_eq!(value["counts"]["representative_solutions"], if triangle { 3 } else { 1 });
        assert_eq!(value["count_status"], if exact.is_some() { "exact" } else { "overflow" });
        assert_eq!(value["counts"]["expanded_solutions"], serde_json::json!(exact.map(|n| n.to_string())));
    }
}

#[test]
fn reruns_and_invalid_options_preserve_accepted_outputs() {
    let p = pairs(1, false);
    let (out, summary) = (tmp("solutions.jsonl"), tmp("summary.json"));
    let args =
        ["enumerate", p.to_str().unwrap(), "--out", out.to_str().unwrap(), "--summary", summary.to_str().unwrap()];
    assert_eq!(mcs(&args).0, 0);
    let before = (std::fs::read(&out).unwrap(), std::fs::read(&summary).unwrap());
    for limit in ["NaN", "inf", "-1", "1e300", "1"] {
        let mut args = args.to_vec();
        args.extend(["--time", limit]);
        assert_eq!(mcs(&args).0, 2, "{limit}");
        assert_eq!((std::fs::read(&out).unwrap(), std::fs::read(&summary).unwrap()), before);
    }
    assert_eq!(
        mcs(&["check", p.to_str().unwrap(), out.to_str().unwrap(), "--summary", summary.to_str().unwrap()]).0,
        0
    );
    std::fs::write(&out, "").unwrap();
    let (code, text) =
        mcs(&["check", p.to_str().unwrap(), out.to_str().unwrap(), "--summary", summary.to_str().unwrap()]);
    assert_eq!(code, 1, "{text}");
    assert_eq!(serde_json::from_str::<serde_json::Value>(&text).unwrap()["summary"]["matches"], false);
}

#[test]
fn summary_failure_does_not_publish_solutions() {
    let p = pairs(1, false);
    let out = tmp("solutions.jsonl");
    let summary = out.parent().unwrap().join("missing/summary.json");
    assert_eq!(
        mcs(&[
            "enumerate",
            p.to_str().unwrap(),
            "--out",
            out.to_str().unwrap(),
            "--summary",
            summary.to_str().unwrap()
        ])
        .0,
        2
    );
    assert!(!out.exists());
    assert!(!summary.exists());
    assert_eq!(
        mcs(&["enumerate", p.to_str().unwrap(), "--out", out.to_str().unwrap(), "--summary", out.to_str().unwrap()]).0,
        2
    );
    assert!(!out.exists());
}

#[test]
fn summary_binds_input_and_solution_bytes() {
    let p = pairs(1, false);
    let (out, summary) = (tmp("solutions.jsonl"), tmp("summary.json"));
    assert_eq!(
        mcs(&[
            "enumerate",
            p.to_str().unwrap(),
            "--out",
            out.to_str().unwrap(),
            "--summary",
            summary.to_str().unwrap()
        ])
        .0,
        0
    );
    assert_eq!(read(&summary)["count_status"], "exact");
    let other = pairs(2, false);
    assert_eq!(
        mcs(&["check", other.to_str().unwrap(), out.to_str().unwrap(), "--summary", summary.to_str().unwrap()]).0,
        1
    );
}

#[test]
fn census_preserves_a_previous_run() {
    let p = pairs(1, false);
    let dir = tmp("census");
    let args = ["census", dir.to_str().unwrap(), "1", p.to_str().unwrap()];
    assert_eq!(mcs(&args).0, 0);
    let before = (std::fs::read(dir.join("census.json")).unwrap(), std::fs::read(dir.join("census.md")).unwrap());
    let data = read(&dir.join("census.json"));
    assert_eq!(data.as_array().unwrap().len(), 4);
    for row in data.as_array().unwrap() {
        assert_eq!(row["count_status"], "exact");
        assert_eq!(row["solutions"], "2");
    }
    assert_eq!(mcs(&args).0, 2);
    assert_eq!(
        (std::fs::read(dir.join("census.json")).unwrap(), std::fs::read(dir.join("census.md")).unwrap()),
        before
    );
}

#[test]
fn oracle_families() {
    for (nm, opts) in [("small-example", vec![]), ("fano", vec!["--twins", "--ord", "unhit"])] {
        let (inp, rf) = (format!("results/oracle/{nm}.json"), format!("results/oracle/{nm}.oracle.json"));
        let (out, sm) = (tmp(&format!("{nm}.jsonl")), tmp(&format!("{nm}.summary.json")));
        let mut a =
            vec!["enumerate", &inp, "--check", "--out", out.to_str().unwrap(), "--summary", sm.to_str().unwrap()];
        a.extend(&opts);
        assert_eq!(mcs(&a).0, 0, "{nm}");
        let s = std::fs::read_to_string(&sm).unwrap();
        assert!(
            s.contains(r#""search_status": "complete""#) && s.contains(r#""verification_status": "passed""#),
            "{s}"
        );
        let (code, o) = mcs(&["check", &inp, out.to_str().unwrap(), "--ref", &rf]);
        assert_eq!(code, 0, "{o}");
        assert!(o.contains(r#""exact_family_match": true"#), "{o}");
    }
}

#[test]
fn exit_status() {
    assert_eq!(mcs(&["enumerate", "results/oracle/fano.json", "--twin"]).0, 2, "unknown option");
    assert_eq!(mcs(&["enumerate", "results/oracle/fano.json", "--k"]).0, 2, "missing value");
    assert_eq!(mcs(&["enumerate", "results/oracle/fano.json", "--time", "-1"]).0, 2, "negative time");
    assert_eq!(mcs(&["enumerate", "results/oracle/fano.oracle.json"]).0, 3, "invalid input");
    let short = tmp("short.json");
    std::fs::write(&short, "[[1, 2], [1, 3], [1, 4], [2, 3]]").unwrap();
    let (code, o) = mcs(&[
        "check",
        "results/oracle/small-example.json",
        "results/oracle/small-example.oracle.json",
        "--ref",
        short.to_str().unwrap(),
    ]);
    assert_eq!(code, 1, "{o}");
    assert!(o.contains(r#""extra": 1"#), "{o}");
    let (code, o) = mcs(&["check", "results/oracle/small-example.json", short.to_str().unwrap()]);
    assert_eq!(code, 0, "a subfamily of valid solutions passes without --ref: {o}");
    let (code, _) = mcs(&["enumerate", "results/oracle/fano.json", "--k", "2", "--count"]);
    assert_eq!(code, 0);
}
