// Legend:
//   Bs  fixed-length bitset over 0..n, stored in 64-bit words w
//   o   the other operand of a binary operation (same length)

#[derive(Clone, PartialEq, Eq, Hash, PartialOrd, Ord, Debug)]
pub struct Bs {
    pub w: Vec<u64>,
}

impl Bs {
    pub fn new(n: usize) -> Self {
        Bs { w: vec![0; n.div_ceil(64)] }
    }
    pub fn of(n: usize, it: impl IntoIterator<Item = usize>) -> Self {
        let mut b = Bs::new(n);
        for i in it {
            b.set(i);
        }
        b
    }
    #[inline]
    pub fn set(&mut self, i: usize) {
        self.w[i >> 6] |= 1 << (i & 63);
    }
    #[inline]
    pub fn clr(&mut self, i: usize) {
        self.w[i >> 6] &= !(1 << (i & 63));
    }
    pub fn cnt(&self) -> u32 {
        self.w.iter().map(|x| x.count_ones()).sum()
    }
    pub fn is0(&self) -> bool {
        self.w.iter().all(|&x| x == 0)
    }
    pub fn and_cnt(&self, o: &Bs) -> u32 {
        self.w.iter().zip(&o.w).map(|(a, b)| (a & b).count_ones()).sum()
    }
    pub fn meets(&self, o: &Bs) -> bool {
        self.w.iter().zip(&o.w).any(|(a, b)| a & b != 0)
    }
    /// self ⊆ o
    pub fn sub(&self, o: &Bs) -> bool {
        self.w.iter().zip(&o.w).all(|(a, b)| a & !b == 0)
    }
    pub fn and(&self, o: &Bs) -> Bs {
        Bs { w: self.w.iter().zip(&o.w).map(|(a, b)| a & b).collect() }
    }
    pub fn andnot_with(&mut self, o: &Bs) {
        self.w.iter_mut().zip(&o.w).for_each(|(a, b)| *a &= !b);
    }
    pub fn or_with(&mut self, o: &Bs) {
        self.w.iter_mut().zip(&o.w).for_each(|(a, b)| *a |= b);
    }
    pub fn first(&self) -> Option<usize> {
        self.w.iter().enumerate().find(|(_, x)| **x != 0).map(|(i, x)| i * 64 + x.trailing_zeros() as usize)
    }
    pub fn iter(&self) -> impl Iterator<Item = usize> + '_ {
        self.w.iter().enumerate().flat_map(|(i, &x)| {
            let mut x = x;
            std::iter::from_fn(move || {
                if x == 0 {
                    return None;
                }
                let t = x.trailing_zeros() as usize;
                x &= x - 1;
                Some(i * 64 + t)
            })
        })
    }
}
