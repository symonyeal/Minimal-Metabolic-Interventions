// The mcs command line against the tracked oracle families in results/oracle/.

use std::path::{Path, PathBuf};
use std::process::Command;

fn mcs(a: &[&str]) -> (i32, String) {
    let o = Command::new(env!("CARGO_BIN_EXE_mcs")).args(a).current_dir(env!("CARGO_MANIFEST_DIR")).output().unwrap();
    (o.status.code().unwrap(), String::from_utf8(o.stdout).unwrap())
}

fn tmp(f: &str) -> PathBuf {
    Path::new(env!("CARGO_TARGET_TMPDIR")).join(f)
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
