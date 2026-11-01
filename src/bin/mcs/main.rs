// mcs: minimal cut sets of metabolic networks.
//
//   mcs enumerate FILE [--ord idx|unhit] [--twins] [--k K] [--count] [--census] [--check]
//                 [--time SECS] [--out SOLS] [--summary JSON]
//   mcs check FILE SOLS [--ref SOLS]
//   mcs census OUT SECS FILE...
//   mcs modes MODEL MODES
//   mcs prepare DIR [--skip-b]
//
// Exit status: 0 success, 1 verification failed, 2 usage or I/O error, 3 invalid input.

mod census;
mod enumerate;
mod inp;
mod prepare;

use sha2::{Digest, Sha256};
use std::collections::{HashMap, HashSet};
use std::str::FromStr;

fn die(s: &str) -> ! {
    eprintln!("mcs: {s}");
    std::process::exit(2)
}

fn hex(b: &[u8]) -> String {
    Sha256::digest(b).iter().map(|x| format!("{x:02x}")).collect()
}

/// Command line split into positionals, valued options and flags; unknown options are errors.
struct Args {
    pos: Vec<String>,
    kv: HashMap<&'static str, String>,
    fl: HashSet<&'static str>,
}

impl Args {
    fn new(a: &[String], np: usize, vals: &[&'static str], flags: &[&'static str]) -> Args {
        let mut r = Args { pos: vec![], kv: HashMap::new(), fl: HashSet::new() };
        let mut it = a.iter();
        while let Some(x) = it.next() {
            if let Some(&k) = vals.iter().find(|&&k| k == x) {
                r.kv.insert(k, it.next().unwrap_or_else(|| die(&format!("{k} needs a value"))).clone());
            } else if let Some(&k) = flags.iter().find(|&&k| k == x) {
                r.fl.insert(k);
            } else if x.starts_with("--") {
                die(&format!("unknown option {x}"));
            } else {
                r.pos.push(x.clone());
            }
        }
        if r.pos.len() < np {
            die("missing arguments; see the usage in README.md");
        }
        r
    }

    fn get(&self, k: &str) -> Option<&str> {
        self.kv.get(k).map(String::as_str)
    }

    fn num<T: FromStr>(&self, k: &str) -> Option<T> {
        self.get(k).map(|x| x.parse().unwrap_or_else(|_| die(&format!("{k}: malformed number {x}"))))
    }

    fn has(&self, k: &str) -> bool {
        self.fl.contains(k)
    }
}

fn main() {
    let a: Vec<String> = std::env::args().skip(1).collect();
    let r = a.get(1..).unwrap_or_default();
    match a.first().map(String::as_str) {
        Some("enumerate") => enumerate::enumerate(r),
        Some("check") => enumerate::check(r),
        Some("census") => census::census(r),
        Some("modes") => prepare::modes(r),
        Some("prepare") => prepare::prepare(r),
        _ => die("usage: mcs enumerate|check|census|modes|prepare ..."),
    }
}
