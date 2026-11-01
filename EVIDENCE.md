# Evidence

Two kinds of measurement are recorded here. Items marked **v0.2.0** were run on 2026-10-31 with
the Rust `mcs` in this repository (Windows 10 laptop, Rust 1.98.1). Items marked **v0.1.0** were run
on 2026-10-27 on a Windows 11 work laptop with external tools (SHD, an LP oracle through SciPy
1.17.1/HiGHS, a NumPy mode check, a Python Berge implementation); their drivers are in the `v0.1.0`
tag and their outputs in `results/`. Times are laptop measurements, not benchmark-runner results.

## Inputs

`mcs prepare data` downloads the sources below, stops on a hash mismatch, and builds every derived
input. Derived hypergraphs are identified by `sha256_normalized` (sorted ids of each normalized
edge, in normalized order), which does not depend on file formatting.

| Input | Source | SHA-256 |
|---|---|---|
| E. coli METATOOL modes, five substrates (lineage A) | `VeraLiconaResearchGroup/MHSGenerationAlgorithms` @ `f043cc3`, `benchmark_inputs/ecoli/*.json` | acetate `bf5f50c5…`, succinate `8ad10e58…`, glycerol `a0b1afc2…`, glucose `a666e621…`, combined `141f94db…` |
| Growth subsets of lineage A | `data/a/*-growth.json` | normalized: acetate `a04c5155…`, succinate `7265d5b8…`, glycerol `5f50ff4f…`, glucose `62a6672f…`, combined `f8f4d7cc…` |
| BiGG `e_coli_core` (Orth et al. 2010; lineage B) | `https://bigg.ucsd.edu/static/models/e_coli_core.json` | `7bedec10576cfe935b19218dc881f3fb14f890a1871448fc19a9b4ee15b448d8` |
| Lineage B elementary modes | `data/b/modes.jsonl` (`mcs modes`) | `a6439f96…` |
| Lineage B hypergraphs | `data/b/ecc-{growth,ethanol,all}.json` | normalized: growth `f5ab6211…`, ethanol `867c98d8…`, all `bd8e7f53…` |

Lineage B is an independent reconstruction of the same organism. A different-organism family was
searched for and not found in tractable form: the smallest non-E. coli BiGG model
(`iAB_RBC_283`, human erythrocyte) has a 130-dimensional flux cone kernel against 24 for
`e_coli_core`, which puts its mode count far beyond enumeration.

## v0.2.0: reproduction of v0.1.0

The Python preparation scripts (`prepare.py`, `bigg_cone.py`, `modes_to_hg.py`) and the separate
`efm` binary were replaced by `mcs prepare` and `mcs modes`. Against the v0.1.0 pipeline run on the
same day:

- `modes.jsonl` is byte-identical (`a6439f96…`; 100,274 modes, 46 two-cycles, 133 columns, rank 63,
  70 free, peak 100,320 rays). Sign propagation prunes the same 8 reactions the LP marked blocked
  (`EX_fru_e`, `EX_fum_e`, `EX_gln__L_e`, `EX_mal__L_e`, `FRUpts2`, `FUMt2_2`, `GLNabc`, `MALt2_2`), and
  these are exactly the reactions in no mode ([`results/modes-ecoli-core.json`](results/modes-ecoli-core.json)).
- All eight derived hypergraphs parse to identical vertices, edges and eligible sets. The v0.1.0 and
  v0.2.0 binaries, each run on both input sets with three configurations (unhit order with twins,
  index order, k = 2), give identical status, counts, nodes, normalized hashes and twin reductions
  in every run that finished. Two runs hit the 120 s limit (index order on B growth and B all), where
  partial counts depend on wall time.
- `mcs census` reproduces all 24 census rows of v0.1.0 exactly (search times excepted); both tables
  are byte-identical.
- The EVIDENCE hashes of v0.1.0 for the growth subsets (`68f10235…` etc.) did not match the files the
  v0.1.0 `prepare.py` wrote; they are replaced by normalized hashes above.

## G1: baseline correctness

- `cargo test --release` (v0.2.0): MMCS in both orders, with and without twins, with the census on
  and off, and with k = 0–3, against a brute-force oracle on 406 hypergraphs (planted twins,
  degenerate and random); checker rejection of omitted, extra, duplicate and nonminimal solutions;
  the input contract and exit statuses; `mcs modes`' double description against brute-force
  elementary-mode supports on 403 cones; pruning, two-cycles and coefficient scaling on a toy network.
- Berge oracle (v0.1.0, `fka.transversal` in `symonyeal/Hypergraph-Dualization`; families in
  `results/oracle/`): exact family match on small-example (5), Fano (7) and acetate-growth (617).
  v0.2.0 repeats all three with `mcs check --ref`; the first two run in `cargo test`.
- SHD (v0.1.0, [`results/anchor-lineage-a-growth.json`](results/anchor-lineage-a-growth.json),
  [`results/anchor-lineage-a-all.json`](results/anchor-lineage-a-all.json)): exact family match
  with Rust MMCS on all five growth instances (617; 4,778; 9,346; 9,374; 73,760) and on four
  unfiltered instances (3,363; 14,136; 25,619; 21,001). Unfiltered combined: only Rust with twins
  and unhit order finished within 15 s (275,914 solutions, 3.2 s); no comparison is claimed for it.
- Competitiveness anchor (v0.1.0), count-only medians of 3 (growth) or 1 (unfiltered) runs: Rust
  MMCS without twins took 1.81 s on combined-growth against SHD 2.50 s (mode 0) and 4.13 s (mode D),
  and 1.10 s against 3.31 s on unfiltered glucose. Small instances are dominated by process start.

## Lineage B: elementary modes and cut sets

- `mcs modes` (v0.2.0): 100,274 elementary modes of `e_coli_core` under its bound signs, in 82 s.
  The count equals the 100,274 EFMs that Clement, Baalhuis, Teusink, Bruggeman, Planqué and de Groot
  report for `e_coli_core` with efmtool (bioRxiv 10.1101/2020.06.06.137554); their medium is not
  stated, ours is the model's default bound signs.
- Mode check (v0.1.0, [`results/verify-modes-ecoli-core.json`](results/verify-modes-ecoli-core.json)):
  all 100,274 supports are distinct and have nullity exactly 1 with a sign-feasible null vector
  nonzero on the support (smallest kept singular value 1.3e-3, largest dropped 6.0e-14); a 500-mode
  random sample also passes in exact rational arithmetic. The modes are unchanged in v0.2.0.
- Hypergraphs: growth (83,601 biomass modes, all reactions but biomass eligible); ethanol (48,770
  `EX_etoh_e` modes, core policy: 69 eligible); all (100,274 modes).
- LP oracle (v0.1.0), independent of the modes: exhaustive LP over the cone found 128 growth cuts
  (86 eligible, 3,742 LPs) and 67 ethanol cuts (69 eligible, 2,416 LPs) of size at most 2, stored in
  `results/oracle/ecc-{growth,ethanol}-k2.oracle.json`. v0.2.0: `mcs enumerate --k 2` then
  `mcs check --ref` gives an exact family match for both.
- `mcs enumerate --check` (v0.2.0): ethanol has 34,389 minimal cut sets (complete, all checked, no
  duplicates).
- Complete counts (v0.2.0, `--count --ord unhit --twins`, [`results/count-ecc-*.json`](results/)):
  growth 889,389 (4.9 s search), ethanol 34,389 (0.2 s), all modes 3,545,843 (19.4 s).
- SHD family (v0.1.0, [`results/shd-family-ecc-ethanol.json`](results/shd-family-ecc-ethanol.json)):
  ethanol family equals SHD's exactly (34,389; Rust 2.5 s, SHD 9.1 s, both writing output). Growth:
  SHD did not finish within 300 s; Rust finished in 22.2 s, 18.2 s of it input normalization of
  83,601 edges. No growth comparison with SHD is claimed.

## G2: opportunity census

`mcs census` (v0.2.0, identical to v0.1.0); full tables in
[`results/census-lineage-a.md`](results/census-lineage-a.md) and
[`results/census-lineage-b.md`](results/census-lineage-b.md). All 24 runs completed. Shares are of
MMCS search nodes. "Exact" counts subtrees rooted at an exact repeat of a completed restricted state.
"WL bound" uses colour-refinement equality, which every isomorphic pair satisfies, so it bounds
isomorphic reuse from above. With twins and unhit order (the comparator configuration):

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
exact memoization in every case, below the 1/6 needed for the 1.2x target. By the predeclared stop
rule ([REPORT.md](REPORT.md) §1), isomorphic continuation reuse stops here. Most remaining
opportunity is in exact repeats and dead subtrees: on B ethanol 13,171 of 41,948 nodes (31%) lie in
exact repeats of dead states. The B growth census did not finish within 240 s and is not reported.

## Known limits

- Normalization's superset removal is quadratic in the worst case and dominates wide inputs
  (10–18 s on lineage B growth, depending on load).
- The census computes full restricted keys at every node: about 100 times the node cost of plain
  search. It is a measurement, not a usable cache.
- Checkpoint and resume are not implemented.
