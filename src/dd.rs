// Extreme rays of the pointed cone {x ∈ ℝ^q : A x = 0, x ≥ 0} by the nullspace double-description
// method with a combinatorial adjacency test (Fukuda–Prodon) over a bit-pattern tree. Exact i128
// arithmetic; overflow is an error rather than a rounding.
//
// Legend:
//   A      integer matrix, q columns (q ≤ 256)
//   d      kernel dimension: free columns of the echelon form, which start as processed constraints
//   Z      zero pattern of a ray over processed columns
//   v      ray values on the remaining (unprocessed) pivot columns
//   rem    remaining pivot columns       i   column processed in this iteration
//   Bpt    bit-pattern tree for "some other ray r has Z(r) ⊇ Z" queries
//
// Rays p, n of {Kλ : x_j ≥ 0 for processed j} are adjacent iff no other ray's zero pattern
// contains Z(p) ∩ Z(n); the necessary count condition |Z(p) ∩ Z(n)| ≥ d − 2 filters first.

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

pub(crate) fn gcd(a: i128, b: i128) -> i128 {
    let (mut a, mut b) = (a.abs(), b.abs());
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a
}

fn norm(v: &mut [i128]) {
    let g = v.iter().fold(0, |g, &x| gcd(g, x));
    if g > 1 {
        v.iter_mut().for_each(|x| *x /= g);
    }
}

const OVF: &str = "i128 overflow";

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

    /// Splits on the bit closest to halving the rays; leaves hold at most 24.
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
            for (b, cb) in c.iter_mut().enumerate() {
                *cb += zg(&rs[r as usize].z, b) as u32;
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

pub struct Dd {
    /// Supports of the extreme rays, as sorted column indices.
    pub rays: Vec<Vec<usize>>,
    pub rank: usize,
    pub free: usize,
    pub peak: usize,
}

pub fn rays(a: &[Vec<i128>], q: usize) -> Result<Dd, String> {
    if q > 256 {
        return Err(format!("{q} columns; at most 256 are supported"));
    }
    let mut m = a.to_vec();

    // Reduced row echelon form, fraction-free with row gcd normalization.
    let mut piv: Vec<(usize, usize)> = vec![]; // (row, column)
    for c in 0..q {
        let r0 = piv.len();
        let Some(p) = (r0..m.len()).find(|&r| m[r][c] != 0) else { continue };
        m.swap(r0, p);
        let pr = m[r0].clone();
        for (r, row) in m.iter_mut().enumerate() {
            if r != r0 && row[c] != 0 {
                let (x, y) = (pr[c], row[c]);
                for (u, &w) in row.iter_mut().zip(&pr) {
                    *u = u.checked_mul(x).and_then(|u| w.checked_mul(y).and_then(|w| u.checked_sub(w))).ok_or(OVF)?;
                }
                norm(row);
            }
        }
        piv.push((r0, c));
    }
    let free: Vec<usize> = (0..q).filter(|c| !piv.iter().any(|&(_, pc)| pc == *c)).collect();
    let d = free.len();
    let mut rem: Vec<usize> = piv.iter().map(|&(_, c)| c).collect();

    // Kernel rays: x_f = L, other free columns 0, x_p = −m[row_p][f] L / m[row_p][p], with L the
    // lcm of the pivots.
    let l = piv.iter().try_fold(1i128, |l, &(r, c)| {
        let x = m[r][c].abs();
        (l / gcd(l, x)).checked_mul(x).ok_or(OVF)
    })?;
    let mut rs: Vec<Ray> = vec![];
    for &f in &free {
        let mut v =
            piv.iter().map(|&(r, c)| (-m[r][f]).checked_mul(l / m[r][c]).ok_or(OVF)).collect::<Result<Vec<_>, _>>()?;
        v.push(l);
        norm(&mut v);
        v.pop();
        let mut z = [0u64; 4];
        free.iter().filter(|&&g| g != f).for_each(|&g| zs(&mut z, g));
        rs.push(Ray { z, v });
    }

    let mut peak = rs.len();
    while !rem.is_empty() {
        // Next column: fewest candidate pairs.
        let pi = (0..rem.len())
            .min_by_key(|&p| {
                let a = rs.iter().filter(|r| r.v[p] > 0).count();
                let b = rs.iter().filter(|r| r.v[p] < 0).count();
                a * b
            })
            .unwrap();
        let i = rem[pi];
        let t = Bpt::new(&rs, q);
        let pos: Vec<u32> = (0..rs.len() as u32).filter(|&r| rs[r as usize].v[pi] > 0).collect();
        let neg: Vec<u32> = (0..rs.len() as u32).filter(|&r| rs[r as usize].v[pi] < 0).collect();
        let mut new: Vec<Ray> = vec![];
        for &p in &pos {
            for &n in &neg {
                let (rp, rn) = (&rs[p as usize], &rs[n as usize]);
                let mut z = zand(&rp.z, &rn.z);
                if (zcnt(&z) as usize) + 2 < d || t.sup(&rs, &z, p, n, 0) {
                    continue;
                }
                let (x, y) = (-rn.v[pi], rp.v[pi]);
                let mut v =
                    rp.v.iter()
                        .zip(&rn.v)
                        .map(|(&u, &w)| {
                            u.checked_mul(x).and_then(|u| w.checked_mul(y).and_then(|w| u.checked_add(w))).ok_or(OVF)
                        })
                        .collect::<Result<Vec<_>, _>>()?;
                norm(&mut v);
                zs(&mut z, i);
                new.push(Ray { z, v });
            }
        }
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
    }
    let rays = rs.iter().map(|r| (0..q).filter(|&c| !zg(&r.z, c)).collect()).collect();
    Ok(Dd { rays, rank: piv.len(), free: d, peak })
}
