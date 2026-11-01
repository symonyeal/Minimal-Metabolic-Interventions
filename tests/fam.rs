// Exact family checks of MMCS against a brute-force oracle on original vertex sets.
//
// Legend:
//   rows  raw edges (sorted dense ids)    n  vertex count    bf  brute-force Tr(rows)
//   go    run one configuration and return its expanded, sorted family

use mcs::chk;
use mcs::hg;
use mcs::srch::{self, Cfg, Ord};

fn bf(n: usize, rows: &[Vec<usize>], k: Option<usize>) -> Vec<Vec<usize>> {
    let hits = |c: &[usize]| rows.iter().all(|e| e.iter().any(|v| c.contains(v)));
    let mut o = vec![];
    for m in 0u32..1 << n {
        let c: Vec<usize> = (0..n).filter(|&i| m >> i & 1 == 1).collect();
        if k.is_none_or(|k| c.len() <= k)
            && hits(&c)
            && c.iter().all(|v| !hits(&c.iter().copied().filter(|x| x != v).collect::<Vec<_>>()))
        {
            o.push(c);
        }
    }
    o.sort();
    o
}

fn go(n: usize, rows: &[Vec<usize>], ord: Ord, tw: bool, k: Option<usize>, census: bool) -> Vec<Vec<usize>> {
    let (h, _) = hg::norm(n, rows);
    let t = hg::twins(&h, tw);
    let cfg = Cfg { ord, k, census, dl: None };
    let mut o: Vec<Vec<usize>> = vec![];
    srch::run(t.h.n, &t.h.ev, &cfg, &mut |s| {
        hg::expand(&t.cls, s, &mut |c| {
            assert!(chk::check(rows, c).is_ok(), "checker rejected {c:?} {ord:?} tw={tw} k={k:?} rows={rows:?}");
            o.push(c.to_vec());
        })
    });
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
        (
            7,
            vec![
                vec![0, 1, 2],
                vec![0, 4, 5],
                vec![0, 3, 6],
                vec![1, 3, 4],
                vec![1, 5, 6],
                vec![2, 3, 5],
                vec![2, 4, 6],
            ],
        ),
        // Planted twins: 0~1 and 3~4~5.
        (7, vec![vec![0, 1, 2], vec![0, 1, 3, 4, 5], vec![2, 3, 4, 5, 6], vec![6]]),
    ];
    for _ in 0..400 {
        let n = 2 + (g.next() % 9) as usize;
        let m = (g.next() % 12) as usize;
        let rows = (0..m)
            .map(|_| {
                let p = 1 + g.next() % 3;
                (0..n).filter(|_| g.next() % 4 < p).collect()
            })
            .collect();
        cs.push((n, rows));
    }
    cs
}

#[test]
fn degenerate() {
    assert_eq!(go(3, &[], Ord::Idx, false, None, false), vec![Vec::<usize>::new()]);
    assert!(go(3, &[vec![]], Ord::Idx, false, None, false).is_empty());
    assert!(go(3, &[vec![]], Ord::Unhit, true, Some(1), true).is_empty());
}

#[test]
fn all_configs_match_brute_force() {
    for (n, rows) in cases() {
        let want = bf(n, &rows, None);
        for ord in [Ord::Idx, Ord::Unhit] {
            for tw in [false, true] {
                for census in [false, true] {
                    assert_eq!(
                        go(n, &rows, ord, tw, None, census),
                        want,
                        "{ord:?} tw={tw} census={census} n={n} rows={rows:?}"
                    );
                }
            }
        }
    }
}

#[test]
fn bounded_matches_brute_force() {
    for (n, rows) in cases().into_iter().take(150) {
        for k in 0..4 {
            let want = bf(n, &rows, Some(k));
            assert_eq!(go(n, &rows, Ord::Idx, true, Some(k), false), want, "k={k}");
            assert_eq!(go(n, &rows, Ord::Unhit, false, Some(k), true), want, "census k={k}");
        }
    }
}

#[test]
fn census_counts_dead_subtrees() {
    // Triangle: root, children 0 and 1, then leaves {0,2}, {1,0}, {1,2}; every subtree emits.
    let cfg = Cfg { ord: Ord::Idx, k: None, census: true, dl: None };
    let (h, _) = hg::norm(3, &[vec![0, 1], vec![0, 2], vec![1, 2]]);
    let st = srch::run(h.n, &h.ev, &cfg, &mut |_| {});
    assert_eq!((st.nodes, st.sols, st.cen.dead, st.done), (6, 3, 0, true));
}

#[test]
fn normalization_counts() {
    let rows = vec![vec![0, 1], vec![0, 1], vec![1], vec![1, 2]];
    let (h, merged) = hg::norm(3, &rows);
    assert_eq!((h.ev.len(), merged), (1, 1));
    assert_eq!(h.ev[0].iter().collect::<Vec<_>>(), vec![1]);
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
