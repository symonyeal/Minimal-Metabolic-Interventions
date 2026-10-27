// efm: elementary modes of an irreversible integer cone {x : A x = 0, x >= 0} by the nullspace
// double-description method with a combinatorial adjacency test (Fukuda–Prodon) over a bit-pattern
// tree. Exact i128 arithmetic; overflow aborts the run rather than rounding.
//
//   efm CONE.json MODES.jsonl SUMMARY.json
//
// Legend:
//   A      integer matrix (rows x q); cone columns are reaction directions from adapters/bigg_cone.py
//   K      kernel basis, identity-scaled on the d free columns (they start as processed constraints)
//   Z      zero pattern of a ray over processed columns (bitmask, q <= 256)
//   v      ray values on the remaining (unprocessed) pivot columns
//   rem    remaining pivot columns       i   column processed in this iteration
//   Bpt    bit-pattern tree for "some other ray r has Z(r) ⊇ Z" queries
//
// Rays of {Kλ : x_j >= 0 for processed j} p, n are adjacent iff no other ray's zero pattern
// contains Z(p) ∩ Z(n); a necessary count condition |Z(p) ∩ Z(n)| >= d - 2 filters first. A
// split reversible reaction yields one two-cycle ray (forward + reverse only); those are dropped.

use serde_json::{Value, json};
use std::io::Write;
use std::time::Instant;

type Z = [u64; 4];

fn zs(z: &mut Z, i: usize) {
    z[i >> 6] |= 1 << (i & 63);
}
fn zg(z: &Z, i: usize) -> bool {
    z[i >> 6] >> (i & 63) & 1 == 1
}
fn zand(a: &Z, b: &Z) -> Z {
    [a[0] & b[0], a[1] & b[1], a[2] & b[2], a[3] & b[3]]
}
fn zsub(a: &Z, b: &Z) -> bool {
    (0..4).all(|k| a[k] & !b[k] == 0)
}
fn zcnt(a: &Z) -> u32 {
    a.iter().map(|x| x.count_ones()).sum()
}

fn gcd(a: i128, b: i128) -> i128 {
    let (mut a, mut b) = (a.abs(), b.abs());
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a
}

fn die(s: &str) -> ! {
    eprintln!("efm: {s}");
    std::process::exit(2)
}

fn norm(v: &mut [i128]) {
    let g = v.iter().fold(0, |g, &x| gcd(g, x));
    if g > 1 {
        v.iter_mut().for_each(|x| *x /= g);
    }
}

struct Ray {
    z: Z,
    v: Vec<i128>,
}

enum Nd {
    Leaf(Vec<u32>),
    Split(usize, usize, usize),
}
struct Bpt {
    uni: Vec<Z>,
    nd: Vec<Nd>,
}

impl Bpt {
    fn new(rs: &[Ray], q: usize) -> Bpt {
        let mut t = Bpt { uni: vec![], nd: vec![] };
        t.build(rs, (0..rs.len() as u32).collect(), q);
        t
    }
    fn build(&mut self, rs: &[Ray], ix: Vec<u32>, q: usize) -> usize {
        let mut u = [0u64; 4];
        for &r in &ix {
            (0..4).for_each(|k| u[k] |= rs[r as usize].z[k]);
        }
        let me = self.nd.len();
        self.uni.push(u);
        self.nd.push(Nd::Leaf(vec![]));
        if ix.len() <= 24 {
            self.nd[me] = Nd::Leaf(ix);
            return me;
        }
        let mut c = vec![0u32; q];
        for &r in &ix {
            for b in 0..q {
                if zg(&rs[r as usize].z, b) {
                    c[b] += 1;
                }
            }
        }
        let h = ix.len() as i64 / 2;
        let b = (0..q).min_by_key(|&b| (c[b] as i64 - h).abs()).unwrap();
        if c[b] == 0 || c[b] as usize == ix.len() {
            self.nd[me] = Nd::Leaf(ix);
            return me;
        }
        let (has, not): (Vec<u32>, Vec<u32>) = ix.into_iter().partition(|&r| zg(&rs[r as usize].z, b));
        let l = self.build(rs, has, q);
        let r = self.build(rs, not, q);
        self.nd[me] = Nd::Split(b, l, r);
        me
    }
    /// Some ray other than a, b has a zero pattern containing z.
    fn sup(&self, rs: &[Ray], z: &Z, a: u32, b: u32, x: usize) -> bool {
        if !zsub(z, &self.uni[x]) {
            return false;
        }
        match &self.nd[x] {
            Nd::Leaf(ix) => ix.iter().any(|&r| r != a && r != b && zsub(z, &rs[r as usize].z)),
            Nd::Split(bit, l, r) => self.sup(rs, z, a, b, *l) || (!zg(z, *bit) && self.sup(rs, z, a, b, *r)),
        }
    }
}

fn main() {
    let t0 = Instant::now();
    let a: Vec<String> = std::env::args().skip(1).collect();
    if a.len() != 3 {
        die("usage: efm CONE.json MODES.jsonl SUMMARY.json");
    }
    let cone: Value = serde_json::from_slice(&std::fs::read(&a[0]).unwrap_or_else(|e| die(&e.to_string()))).unwrap_or_else(|e| die(&e.to_string()));
    let cols: Vec<(String, i64)> = cone["columns"].as_array().unwrap().iter().map(|c| (c["reaction"].as_str().unwrap().to_owned(), c["dir"].as_i64().unwrap())).collect();
    let mut m: Vec<Vec<i128>> = cone["A"].as_array().unwrap().iter().map(|r| r.as_array().unwrap().iter().map(|x| x.as_i64().unwrap() as i128).collect()).collect();
    let q = cols.len();
    if q > 256 {
        die("more than 256 cone columns");
    }

    // Reduced row echelon form, fraction-free with row gcd normalization.
    let mut piv: Vec<(usize, usize)> = vec![]; // (row, column)
    let mut r0 = 0;
    for c in 0..q {
        let Some(p) = (r0..m.len()).find(|&r| m[r][c] != 0) else { continue };
        m.swap(r0, p);
        for r in 0..m.len() {
            if r != r0 && m[r][c] != 0 {
                let (x, y) = (m[r0][c], m[r][c]);
                for k in 0..q {
                    m[r][k] = m[r][k].checked_mul(x).and_then(|u| m[r0][k].checked_mul(y).and_then(|w| u.checked_sub(w))).unwrap_or_else(|| die("overflow in elimination"));
                }
                norm(&mut m[r]);
            }
        }
        piv.push((r0, c));
        r0 += 1;
    }
    let isp: Vec<bool> = (0..q).map(|c| piv.iter().any(|&(_, pc)| pc == c)).collect();
    let free: Vec<usize> = (0..q).filter(|&c| !isp[c]).collect();
    let d = free.len();
    let mut rem: Vec<usize> = piv.iter().map(|&(_, c)| c).collect();

    // Kernel rays: x_f = L, x_p = -m[row_p][f] * L / m[row_p][p].
    let mut rs: Vec<Ray> = free
        .iter()
        .map(|&f| {
            let l = piv.iter().fold(1i128, |l, &(r, c)| {
                let x = m[r][c].abs();
                (l / gcd(l, x)).checked_mul(x).unwrap_or_else(|| die("overflow in kernel"))
            });
            let mut v: Vec<i128> = piv.iter().map(|&(r, c)| (-m[r][f]).checked_mul(l / m[r][c]).unwrap_or_else(|| die("overflow in kernel"))).collect();
            v.push(l);
            norm(&mut v);
            v.pop();
            let mut z = [0u64; 4];
            for &g in &free {
                if g != f {
                    zs(&mut z, g);
                }
            }
            Ray { z, v }
        })
        .collect();

    let mut peak = rs.len();
    let mut log = vec![];
    while !rem.is_empty() {
        // Next column: fewest candidate pairs.
        let (pi, _) = (0..rem.len())
            .map(|p| {
                let (mut a, mut b) = (0u64, 0u64);
                for r in &rs {
                    if r.v[p] > 0 {
                        a += 1
                    } else if r.v[p] < 0 {
                        b += 1
                    }
                }
                (p, a * b)
            })
            .min_by_key(|&(_, x)| x)
            .unwrap();
        let i = rem[pi];
        let t = Bpt::new(&rs, q);
        let pos: Vec<u32> = (0..rs.len() as u32).filter(|&r| rs[r as usize].v[pi] > 0).collect();
        let neg: Vec<u32> = (0..rs.len() as u32).filter(|&r| rs[r as usize].v[pi] < 0).collect();
        let mut new: Vec<Ray> = vec![];
        let mut tried = 0u64;
        for &p in &pos {
            for &n in &neg {
                let z = zand(&rs[p as usize].z, &rs[n as usize].z);
                if (zcnt(&z) as usize) + 2 < d {
                    continue;
                }
                tried += 1;
                if t.sup(&rs, &z, p, n, 0) {
                    continue;
                }
                let (rp, rn) = (&rs[p as usize], &rs[n as usize]);
                let (x, y) = (-rn.v[pi], rp.v[pi]);
                let mut v: Vec<i128> = rp
                    .v
                    .iter()
                    .zip(&rn.v)
                    .map(|(&u, &w)| u.checked_mul(x).and_then(|u| w.checked_mul(y).and_then(|w| u.checked_add(w))).unwrap_or_else(|| die("overflow in combination")))
                    .collect();
                norm(&mut v);
                let mut z = z;
                zs(&mut z, i);
                new.push(Ray { z, v });
            }
        }
        let (np, nn, nz) = (pos.len(), neg.len(), rs.len() - pos.len() - neg.len());
        let mut keep: Vec<Ray> = rs.into_iter().filter(|r| r.v[pi] >= 0).collect();
        for r in keep.iter_mut() {
            if r.v[pi] == 0 {
                zs(&mut r.z, i);
            }
        }
        keep.extend(new);
        for r in keep.iter_mut() {
            r.v.swap_remove(pi);
        }
        rem.swap_remove(pi);
        rs = keep;
        peak = peak.max(rs.len());
        log.push(json!({"column": i, "pos": np, "neg": nn, "zero": nz, "pairs_tested": tried, "rays": rs.len(), "s": t0.elapsed().as_secs_f64()}));
        eprintln!("col {i}: +{np} -{nn} 0{nz} tested {tried} -> {} rays ({:.1}s)", rs.len(), t0.elapsed().as_secs_f64());
    }

    // Supports by reaction; drop two-cycles.
    let mut w = std::io::BufWriter::new(std::fs::File::create(&a[1]).unwrap_or_else(|e| die(&e.to_string())));
    let (mut nm, mut cyc) = (0u64, 0u64);
    for r in &rs {
        let sup: Vec<usize> = (0..q).filter(|&c| !zg(&r.z, c)).collect();
        let mut ids: Vec<&str> = sup.iter().map(|&c| cols[c].0.as_str()).collect();
        ids.sort_unstable();
        ids.dedup();
        if ids.len() == 1 && sup.len() == 2 {
            cyc += 1;
            continue;
        }
        let dirs: Vec<Value> = sup.iter().map(|&c| json!([cols[c].0, cols[c].1])).collect();
        writeln!(w, "{}", json!({"support": ids, "columns": dirs})).unwrap();
        nm += 1;
    }
    w.flush().unwrap();
    let sm = json!({"tool": "efm", "cone": a[0], "columns": q, "rank": piv.len(), "free": d,
                    "modes": nm, "two_cycles": cyc, "peak_rays": peak, "elapsed_s": t0.elapsed().as_secs_f64(), "iterations": log});
    std::fs::write(&a[2], serde_json::to_string_pretty(&sm).unwrap()).unwrap();
    println!("{}", json!({"modes": nm, "two_cycles": cyc, "peak_rays": peak, "elapsed_s": t0.elapsed().as_secs_f64()}));
}
