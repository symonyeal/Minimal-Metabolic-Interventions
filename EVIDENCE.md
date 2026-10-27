# Evidence

Every number below comes from a run on 2026-10-27 on the Windows 11 work laptop (Rust 1.99.1,
Python 3.14, SciPy 1.17.1, NumPy 2.4.4, HiGHS through `linprog`). Times are laptop smoke
measurements including process start, not benchmark-runner results. `scripts/prepare.py` fetches
and rebuilds every input below into `data/` and stops on a hash mismatch.

## Inputs

| Input | Source | SHA-256 |
|---|---|---|
| E. coli METATOOL modes, five substrates (lineage A) | `VeraLiconaResearchGroup/MHSGenerationAlgorithms` @ `f043cc3`, `benchmark_inputs/ecoli/*.json` → `data/a/` | acetate `bf5f50c5…`, succinate `8ad10e58…`, glycerol `a0b1afc2…`, glucose `a666e621…`, combined `141f94db…` |
| Growth-mode subsets of lineage A | `data/a/*-growth.json` (`scripts/prepare.py`) | acetate `68f10235…`, succinate `9245a42e…`, glycerol `8c1b504f…`, glucose `9f9b7d23…`, combined `1c8955f2…` |
| BiGG `e_coli_core` (Orth et al. 2010; lineage B) | `https://bigg.ucsd.edu/static/models/e_coli_core.json` → `data/b/` | `7bedec10576cfe935b19218dc881f3fb14f890a1871448fc19a9b4ee15b448d8` |
| Lineage B integer cone | `adapters/bigg_cone.py` | `025c8051…` |
| Lineage B elementary modes | `efm` | `a6439f96…` |
| SHD binary (Uno), MSVC `/std:clatest`, from `containers/shd/src` at the same commit | not bundled | `42cc79f1…` |

Lineage B is an independent reconstruction of the same organism. A different-organism family was
searched for and not found in tractable form: the smallest non-E. coli BiGG model
(`iAB_RBC_283`, human erythrocyte) has a 130-dimensional flux cone kernel against 24 for
`e_coli_core`, which puts its mode count far beyond enumeration.

## G1: baseline correctness

- `cargo test --release`: 5 test groups pass. They compare MMCS (index and unhit order), RS, twins
  on and off, census on, and k = 0–3 with a brute-force oracle on 406 hypergraphs (planted twins,
  degenerate and random), plus checker rejection of omitted, extra, duplicate and nonminimal
  solutions and the input contract.
- Python Berge oracle (`fka.transversal` in `symonyeal/Hypergraph-Dualization`; families in
  `results/oracle/`): exact family match for 5 configurations on small-example (5), Fano (7) and
  acetate-growth (617), checked with `mcs check FILE SOLS --ref results/oracle/NAME.oracle.json`.
- SHD (`bench/anchor.py`): exact family match with Rust MMCS on all five growth instances
  (617; 4,778; 9,346; 9,374; 73,760) and on four unfiltered instances (3,363; 14,136; 25,619;
  21,001). Unfiltered combined: only Rust with twins and unhit order finished within 15 s
  (275,914 solutions, 3.2 s); no comparison is claimed for it.
- Competitiveness anchor, count-only medians of 3 (growth) or 1 (unfiltered) runs: Rust MMCS
  without twins takes 1.81 s on combined-growth against SHD 2.50 s (mode 0) and 4.13 s (mode D),
  and 1.10 s against 3.31 s on unfiltered glucose. Small instances are dominated by process start.
- RS is not competitive (glycerol unfiltered: over 15 s against 0.51 s for MMCS); it lacks the
  published violator pruning and is a correctness comparator only.

## Lineage B: elementary modes and cut sets

- `efm`: 100,274 elementary modes of `e_coli_core` under its bound signs, in 70 s, peak 100,320
  rays; 46 two-cycles from split reversible reactions dropped. The count equals the 100,274 EFMs
  that Clement, Baalhuis, Teusink, Bruggeman, Planqué and de Groot report for `e_coli_core` with
  efmtool (bioRxiv 10.1101/2020.06.06.137554); their medium is not stated, ours is the model's
  default bound signs.
- `adapters/verify_modes.py`: all 100,274 supports are distinct and have nullity exactly 1 with a
  sign-feasible null vector nonzero on the support (smallest kept singular value 1.3e-3, largest
  dropped 6.0e-14); a 500-mode random sample also passes in exact rational arithmetic.
- Hypergraphs (`adapters/modes_to_hg.py`): growth (83,601 biomass modes, all reactions but biomass
  eligible); ethanol (48,770 `EX_etoh_e` modes, core policy: 69 eligible); all (100,274 modes).
- `adapters/lp_cuts.py`, independent of the modes: exhaustive LP over the cone finds 128 growth and
  67 ethanol minimal cut sets of size at most 2; `mcs check --ref` gives an exact family match for
  both.
- `mcs enumerate --check`: ethanol has 34,389 minimal cut sets (complete, all checked).
- Complete counts, `--count --ord unhit --twins`: growth 889,389 (3.4 s search), ethanol 34,389
  (0.2 s), all modes 3,545,843 (19.3 s).
- `bench/shd_family.py`: ethanol family equals SHD's exactly (34,389; Rust 2.5 s, SHD 9.1 s, both
  writing output). Growth: SHD did not finish within 300 s; Rust finished in 22.2 s (18.2 s of it
  input normalization of 83,601 edges), peak working set 349 MB. No growth comparison with SHD is
  claimed.

## G2: opportunity census

`bench/census.py`; full tables in `results/census-lineage-a.md` and `results/census-lineage-b.md`.
All 24 runs completed. Shares are of MMCS search nodes. "Exact" counts subtrees rooted at an exact
repeat of a completed restricted state. "WL bound" uses colour-refinement equality, which every
isomorphic pair satisfies, so it bounds isomorphic reuse from above. With twins and unhit order
(the comparator configuration):

| Instance | Nodes | Dead | Exact | WL bound | WL − exact |
|---|---:|---:|---:|---:|---:|
| A acetate-growth | 77 | 26.0% | 1.3% | 5.2% | 3.9% |
| A succinate-growth | 476 | 60.3% | 15.8% | 22.3% | 6.5% |
| A glycerol-growth | 1,162 | 62.8% | 18.2% | 25.3% | 7.1% |
| A glucose-growth | 3,901 | 76.6% | 21.0% | 28.0% | 6.9% |
| A combined-growth | 38,391 | 79.3% | 15.7% | 19.8% | 4.1% |
| B ethanol | 41,948 | 77.8% | 35.9% | 42.7% | 6.8% |

Without twins the WL − exact gap is 14–27%. Static twins already remove most of the symmetry an
isomorphism test would find. After twins, isomorphism adds at most 7.1 percentage points over
exact memoization in every case, below the 1/6 needed for the 1.2x target. By the PRD stop rule,
isomorphic continuation reuse (claim item 1) stops here. Most remaining opportunity is in exact
repeats and dead subtrees: on B ethanol 13,171 of 41,948 nodes (31%) lie in exact repeats of dead
states. The B growth census did not finish within 240 s and is not reported.

## Known limits

- Normalization's superset removal is quadratic in the worst case and dominates wide inputs
  (18 s of 22 s on lineage B growth).
- The census computes full restricted keys at every node: about 100 times the node cost of plain
  search. It is a measurement, not a usable cache.
- Peak memory is measured on Windows only. Checkpoint/resume (F12) is not implemented.
