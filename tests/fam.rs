// Exact family checks against a brute-force oracle on original vertex sets.
//
// Legend:
//   rows  raw edges (dense ids)    n  vertex count    bf  brute-force Tr(rows)
//   go    run one configuration and return its expanded, sorted family

use mcs::chk;
use mcs::hg::{self, Raw};
use mcs::srch::{self, Alg, Cfg, Ord};

fn raw(n: usize, rows: &[Vec<usize>]) -> Raw {
    Raw { ids: (0..n).map(|i| i.to_string()).collect(), lab: vec![None; n], rows: rows.iter().map(|e| { let mut e = e.clone(); e.sort_unstable(); e.dedup(); e }).collect(), dup: 0, fmt: "test" }
}

fn bf(n: usize, rows: &[Vec<usize>], k: Option<usize>) -> Vec<Vec<usize>> {
    let mut o = vec![];
    for m in 0u32..1 << n {
        let c: Vec<usize> = (0..n).filter(|&i| m >> i & 1 == 1).collect();
        if k.is_some_and(|k| c.len() > k) {
            continue;
        }
        let hits = |c: &[usize]| rows.iter().all(|e| e.iter().any(|v| c.contains(v)));
        if hits(&c) && c.iter().all(|v| !hits(&c.iter().copied().filter(|x| x != v).collect::<Vec<_>>())) {
            o.push(c);
        }
    }
    o.sort();
    o
}

fn go(n: usize, rows: &[Vec<usize>], alg: Alg, ord: Ord, tw: bool, k: Option<usize>, census: bool) -> Vec<Vec<usize>> {
    let r = raw(n, rows);
    let h = hg::norm(&r);
    let t = hg::twins(&h, tw);
    let cfg = Cfg { alg, ord, k, census, cap: None, dl: None };
    let mut o: Vec<Vec<usize>> = vec![];
    srch::run(t.h.n, &t.h.ev, &cfg, &mut |s| hg::expand(&t.cls, s, &mut |c| {
        assert!(chk::check(&r.rows, c).is_ok(), "checker rejected {c:?} {alg:?} {ord:?} tw={tw} k={k:?} rows={rows:?}");
        o.push(c.to_vec());
        true
    }));
    let l = o.len();
    o.sort();
    o.dedup();
    assert_eq!(l, o.len(), "duplicate solutions");
    o
}

struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }
}

fn cases() -> Vec<(usize, Vec<Vec<usize>>)> {
    let mut g = Rng(0x9e3779b97f4a7c15);
    let mut cs = vec![
        (3, vec![]),
        (3, vec![vec![]]),
        (3, vec![vec![0, 1], vec![]]),
        (5, vec![vec![0, 1, 4], vec![1, 2, 3], vec![0, 2]]),
        (7, vec![vec![0, 1, 2], vec![0, 4, 5], vec![0, 6, 3], vec![1, 3, 4], vec![1, 5, 6], vec![2, 3, 5], vec![2, 4, 6]]),
        // Planted twins with distinct labels: 0~1 and 3~4~5.
        (7, vec![vec![0, 1, 2], vec![0, 1, 3, 4, 5], vec![2, 3, 4, 5, 6], vec![6]]),
    ];
    for _ in 0..400 {
        let n = 2 + (g.next() % 9) as usize;
        let m = (g.next() % 12) as usize;
        let rows = (0..m)
            .map(|_| {
                let p = 1 + g.next() % 3;
                (0..n).filter(|_| g.next() % 4 < p).collect::<Vec<_>>()
            })
            .collect();
        cs.push((n, rows));
    }
    cs
}

#[test]
fn degenerate() {
    assert_eq!(go(3, &[], Alg::Mmcs, Ord::Idx, false, None, false), vec![Vec::<usize>::new()]);
    assert!(go(3, &[vec![]], Alg::Mmcs, Ord::Idx, false, None, false).is_empty());
    assert!(go(3, &[vec![]], Alg::Rs, Ord::Idx, true, None, false).is_empty());
}

#[test]
fn all_configs_match_brute_force() {
    for (n, rows) in cases() {
        let want = bf(n, &rows, None);
        for alg in [Alg::Mmcs, Alg::Rs] {
            for ord in [Ord::Idx, Ord::Unhit] {
                if alg == Alg::Rs && ord == Ord::Unhit {
                    continue;
                }
                for tw in [false, true] {
                    assert_eq!(go(n, &rows, alg, ord, tw, None, false), want, "{alg:?} {ord:?} tw={tw} n={n} rows={rows:?}");
                }
            }
        }
        // Census instrumentation must not change the family (F09).
        assert_eq!(go(n, &rows, Alg::Mmcs, Ord::Unhit, true, None, true), want);
    }
}

#[test]
fn bounded_matches_brute_force() {
    for (n, rows) in cases().into_iter().take(150) {
        for k in 0..4 {
            let want = bf(n, &rows, Some(k));
            assert_eq!(go(n, &rows, Alg::Mmcs, Ord::Idx, true, Some(k), false), want, "k={k}");
            assert_eq!(go(n, &rows, Alg::Rs, Ord::Idx, false, Some(k), false), want, "rs k={k}");
        }
    }
}

#[test]
fn checker_rejects_bad_solutions() {
    let rows = vec![vec![0, 1, 4], vec![1, 2, 3], vec![0, 2]];
    assert_eq!(chk::check(&rows, &[0, 1]).map(|_| ()), Ok(()));
    assert_eq!(chk::check(&rows, &[0]), Err(chk::Bad::Uncovered(1)));
    assert_eq!(chk::check(&rows, &[0, 1, 2]), Err(chk::Bad::NoWitness(0)));
    let want = bf(5, &rows, None);
    let mut omit = want.clone();
    omit.pop();
    let c = chk::cmp(&omit, &want);
    assert_eq!((c.missing, c.extra), (1, 0));
    let mut extra = want.clone();
    extra.push(vec![0, 1, 2]);
    assert_eq!(chk::cmp(&extra, &want).extra, 1);
    let mut dup = want.clone();
    dup.push(want[0].clone());
    assert_eq!(chk::cmp(&dup, &want).dup_a, 1);
}

#[test]
fn input_contract() {
    assert!(hg::parse(br#"{"schema":"mcs-hg/2","vertices":[],"edges":[]}"#).is_err());
    assert!(hg::parse(br#"{"schema":"mcs-hg/1","vertices":[{"id":"a"}],"edges":[["b"]]}"#).is_err());
    assert!(hg::parse(br#"{"schema":"mcs-hg/1","vertices":[{"id":"a","label":"x"},{"id":"a","label":"y"}],"edges":[]}"#).is_err());
    assert!(hg::parse(br#"{"sets":[[1,2.5]]}"#).is_err());
    let r = hg::parse(br#"{"schema":"mcs-hg/1","vertices":[{"id":"a"},{"id":"b"},{"id":"c"}],"edges":[["a","a","b"],["c"]],"eligible":["a","b"]}"#).unwrap();
    assert_eq!(r.dup, 1);
    assert_eq!(r.rows, vec![vec![0, 1], vec![]]);
    assert!(hg::norm(&r).empty);
    let r = hg::parse(br#"{"sets":[[3,1],[1],[1,3]],"labels":{"7":"x"}}"#).unwrap();
    assert_eq!(r.ids, vec!["1", "3", "7"]);
    let h = hg::norm(&r);
    assert_eq!((h.ev.len(), h.dom.len()), (1, 2));
}
