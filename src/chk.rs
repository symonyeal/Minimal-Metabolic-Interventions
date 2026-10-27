// Independent solution checks. Deliberately shares no state or bitsets with the search: it
// recomputes intersections from the raw input rows (eligibility applied, before normalization).
//
// Legend:
//   rows  raw edges as sorted dense vertex lists     c  candidate solution, sorted dense ids
//   w[i]  index of a raw row meeting c only at c[i] (private witness)

use std::collections::BTreeSet;

#[derive(Debug, PartialEq)]
pub enum Bad {
    Uncovered(usize),
    NoWitness(usize),
}

fn meet(a: &[usize], b: &[usize]) -> (usize, usize) {
    // (intersection size capped at 2, last common element)
    let (mut i, mut j, mut k, mut x) = (0, 0, 0, usize::MAX);
    while i < a.len() && j < b.len() && k < 2 {
        match a[i].cmp(&b[j]) {
            std::cmp::Ordering::Less => i += 1,
            std::cmp::Ordering::Greater => j += 1,
            std::cmp::Ordering::Equal => {
                k += 1;
                x = a[i];
                i += 1;
                j += 1;
            }
        }
    }
    (k, x)
}

/// Coverage of every row and a private witness row for every selected vertex.
pub fn check(rows: &[Vec<usize>], c: &[usize]) -> Result<Vec<usize>, Bad> {
    let mut w = vec![usize::MAX; c.len()];
    for (r, e) in rows.iter().enumerate() {
        match meet(e, c) {
            (0, _) => return Err(Bad::Uncovered(r)),
            (1, x) => {
                let i = c.binary_search(&x).unwrap();
                if w[i] == usize::MAX {
                    w[i] = r;
                }
            }
            _ => {}
        }
    }
    match w.iter().position(|&r| r == usize::MAX) {
        Some(i) => Err(Bad::NoWitness(c[i])),
        None => Ok(w),
    }
}

pub struct Cmp {
    pub missing: usize,
    pub extra: usize,
    pub dup_a: usize,
    pub dup_b: usize,
}

/// Exact family comparison of a (tested) against b (reference).
pub fn cmp(a: &[Vec<usize>], b: &[Vec<usize>]) -> Cmp {
    let sa: BTreeSet<&Vec<usize>> = a.iter().collect();
    let sb: BTreeSet<&Vec<usize>> = b.iter().collect();
    Cmp {
        missing: sb.difference(&sa).count(),
        extra: sa.difference(&sb).count(),
        dup_a: a.len() - sa.len(),
        dup_b: b.len() - sb.len(),
    }
}
