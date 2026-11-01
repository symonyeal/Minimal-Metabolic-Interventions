// Elementary modes: double description against a brute-force oracle, and network construction.
//
// Legend:
//   A     integer matrix (r x q)       S   candidate support (column subset)
//   null1 the null vector of A[:, S] when that null space has dimension exactly 1
//
// S supports an extreme ray of {A x = 0, x ≥ 0} iff A[:, S] has a one-dimensional null space
// spanned by a vector of one strict sign.

use mcs::dd;
use mcs::net::{self, Net, Pol, Rxn};

fn null1(a: &[Vec<i128>], s: &[usize]) -> Option<Vec<i128>> {
    let mut m: Vec<Vec<i128>> = a.iter().map(|r| s.iter().map(|&j| r[j]).collect()).collect();
    let k = s.len();
    let mut piv: Vec<usize> = vec![];
    for c in 0..k {
        let r0 = piv.len();
        let Some(p) = (r0..m.len()).find(|&r| m[r][c] != 0) else { continue };
        m.swap(r0, p);
        let pr = m[r0].clone();
        for (r, row) in m.iter_mut().enumerate() {
            if r != r0 && row[c] != 0 {
                let (x, y) = (pr[c], row[c]);
                row.iter_mut().zip(&pr).for_each(|(u, &w)| *u = *u * x - w * y);
            }
        }
        piv.push(c);
    }
    if k - piv.len() != 1 {
        return None;
    }
    let f = (0..k).find(|c| !piv.contains(c)).unwrap();
    let l: i128 = piv.iter().enumerate().map(|(r, &c)| m[r][c].abs()).product();
    let mut v = vec![0; k];
    v[f] = l;
    for (r, &c) in piv.iter().enumerate() {
        v[c] = -m[r][f] * l / m[r][c];
    }
    Some(v)
}

fn brute(a: &[Vec<i128>], q: usize) -> Vec<Vec<usize>> {
    let mut o = vec![];
    for mask in 1u32..1 << q {
        let s: Vec<usize> = (0..q).filter(|&j| mask >> j & 1 == 1).collect();
        if let Some(v) = null1(a, &s)
            && (v.iter().all(|&x| x > 0) || v.iter().all(|&x| x < 0))
        {
            o.push(s);
        }
    }
    o.sort();
    o
}

#[test]
fn rays_match_brute_force() {
    let mut g: u64 = 0x2545f4914f6cdd1d;
    let mut next = |b: u64| {
        g ^= g << 13;
        g ^= g >> 7;
        g ^= g << 17;
        g % b
    };
    let mut cases: Vec<(Vec<Vec<i128>>, usize)> =
        vec![(vec![vec![1, -1, 0], vec![0, 1, -1]], 3), (vec![], 2), (vec![vec![1, 1]], 2)];
    for _ in 0..400 {
        let (r, q) = (1 + next(4) as usize, 2 + next(7) as usize);
        cases.push(((0..r).map(|_| (0..q).map(|_| next(5) as i128 - 2).collect()).collect(), q));
    }
    let mut nonempty = 0;
    for (a, q) in cases {
        let mut got = dd::rays(&a, q).unwrap().rays;
        got.sort();
        let n = got.len();
        got.dedup();
        assert_eq!(n, got.len(), "duplicate rays {a:?}");
        assert_eq!(got, brute(&a, q), "A = {a:?}");
        nonempty += (n > 0) as usize;
    }
    assert!(nonempty > 100, "too few cases with rays: {nonempty}");
}

fn rx(id: &str, lo: f64, hi: f64, st: &[(usize, f64)]) -> Rxn {
    Rxn { id: id.into(), name: None, lo, hi, st: st.to_vec() }
}

#[test]
fn network_cone_and_modes() {
    // a_e ⇄ (exchange), a_e → a, a → 0.5 b, b →, and c → b with nothing producing c.
    let net = Net {
        mets: vec!["a_e".into(), "a".into(), "b".into(), "c".into()],
        rx: vec![
            rx("EX_a", -10.0, 10.0, &[(0, -1.0)]),
            rx("Ta", 0.0, 10.0, &[(0, -1.0), (1, 1.0)]),
            rx("R", 0.0, 10.0, &[(1, -1.0), (2, 0.5)]),
            rx("EX_b", 0.0, 10.0, &[(2, -1.0)]),
            rx("Rc", 0.0, 10.0, &[(3, -1.0), (2, 1.0)]),
        ],
    };
    let c = net::cone(&net).unwrap();
    assert_eq!(c.pruned, vec![4]);
    assert_eq!(c.cols, vec![(0, 1), (0, -1), (1, 1), (2, 1), (3, 1)]);
    assert_eq!(c.a.len(), 3);
    assert_eq!(c.a[1][3], -2, "R scaled by 2: a −2, b +1");
    let d = dd::rays(&c.a, c.cols.len()).unwrap();
    let m = net::modes(&net, &c, d.rays);
    assert_eq!((m.sup.clone(), m.two_cycles, m.blocked.clone()), (vec![vec![0, 1, 2, 3]], 1, vec![4]));
    let (es, vs, el) = net::scenario(&net, &m.sup, Some(3), Pol::Core);
    assert_eq!((es.len(), vs, el), (1, vec![0, 1, 2, 3], vec![1, 2]));
}
