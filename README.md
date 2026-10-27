# Minimal metabolic interventions

`mcs` enumerates the inclusion-minimal hitting sets of an explicit hypergraph with MMCS or RS,
checks every solution independently, and records where MMCS spends its search. `efm` computes the
elementary modes of small metabolic networks exactly. Applied here to minimal reaction cut sets of
two E. coli network reconstructions.

## Background

A minimal cut set of a metabolic network is a set of reactions whose removal blocks a target
function, such as growth or ethanol export, and that contains no smaller blocking set. If the
elementary modes carrying the target are known, the minimal cut sets are exactly the minimal
hitting sets of the modes' reaction supports (Haus, Klamt and Stephen 2008). Enumerating them is
hypergraph dualization, for which MMCS (Murakami and Uno 2014) is the standard practical algorithm.

## Result

The census in [REPORT.md](REPORT.md) measures, without pruning, which MMCS search nodes lie in
subtrees that emit nothing (dead) or that repeat an earlier completed state. With static twins
merged and the unhit-degree candidate order:

| Instance | Nodes | Dead | Exact repeats | Isomorphic repeats (upper bound) |
|---|---:|---:|---:|---:|
| A acetate-growth | 77 | 26.0% | 1.3% | 5.2% |
| A succinate-growth | 476 | 60.3% | 15.8% | 22.3% |
| A glycerol-growth | 1,162 | 62.8% | 18.2% | 25.3% |
| A glucose-growth | 3,901 | 76.6% | 21.0% | 28.0% |
| A combined-growth | 38,391 | 79.3% | 15.7% | 19.8% |
| B ethanol | 41,948 | 77.8% | 35.9% | 42.7% |

Isomorphism adds at most 7.1 percentage points over exact repeats. The study's hypothesis, reuse of
isomorphic search states, needed 1/6 and was stopped at that gate. Most of the search lies in dead
subtrees. Lineage A is the METATOOL network of Haus, Klamt and Stephen on five substrates; lineage B
is BiGG `e_coli_core`, whose 100,274 elementary modes `efm` reproduces.

Every enumerated family used in these tables was checked against at least one independent
enumerator (brute force, a Berge implementation, SHD or an LP cut oracle); see
[EVIDENCE.md](EVIDENCE.md). On combined-growth, `mcs` without twins takes 1.81 s against 2.50 s for
SHD on the same laptop.

## Usage

```text
cargo build --release
cargo test --release
python -m pip install -r requirements.txt   # adapters only: numpy, scipy
python scripts/prepare.py                   # fetch and build every input into data/
```

```text
mcs inspect FILE
mcs enumerate FILE [--alg mmcs|rs] [--ord idx|unhit] [--twins] [--k K] [--count]
              [--out SOLS] [--summary JSON] [--census] [--time SECS] [--cap N] [--check]
mcs check FILE SOLS [--ref SOLS]
efm CONE.json MODES.jsonl SUMMARY.json
```

| Option | Effect |
|---|---|
| `--alg mmcs` | MMCS with the minimum-intersection branching edge and violator pruning (default) |
| `--alg rs` | RS reverse search, without violator pruning; a correctness comparator |
| `--ord unhit` | Candidates by increasing unhit degree (Hörmann and Schirneck 2026), reimplemented |
| `--twins` | Merge vertices with identical incidence; solutions are expanded over their classes |
| `--k K` | Only minimal hitting sets of at most K vertices (`bounded_complete`) |
| `--check` | Verify each expanded solution against the raw rows: coverage and a private witness row per vertex |
| `--census` | Record dead and repeated subtrees; MMCS only |

Input is `mcs-hg/1` JSON, or a publisher file `{"sets": [[1, 2], ...], "labels": {"1": "mue"}}`:

```json
{"schema": "mcs-hg/1",
 "vertices": [{"id": "PFK", "label": "Phosphofructokinase"}],
 "edges": [["PFK", "PGI"]],
 "eligible": ["PFK"],
 "meta": {}}
```

Unknown ids, conflicting labels, malformed numbers and other schema versions are rejected. Edges
are intersected with `eligible`. An empty edge means no solution exists; no edges means the empty
set is the only solution. Duplicate and superset edges are removed with their source rows recorded.

`--out` writes one JSON array of original ids per line. The summary reports `search_status`
(`complete`, `bounded_complete`, `resource_limited`, `invalid_input`), `verification_status`
(`not_run`, `passed`, `failed`), hashes of the input bytes and normalized form, settings,
representative and expanded counts, nodes, stage times, the census and peak memory (Windows peak
working set; not measured elsewhere). It is written last through a rename; a missing summary means
an unfinished run.

## Reproduction

| Table or claim | Command |
|---|---|
| Tests | `cargo test --release` |
| Inputs, lineage B modes | `python scripts/prepare.py` |
| Mode check | `python adapters/verify_modes.py data/b/e_coli_core.json data/b/modes.jsonl 500 OUT.json` |
| LP oracle, cuts of size ≤ 2 | `python adapters/lp_cuts.py data/b/e_coli_core.json data/b/ecc-growth.json 2 OUT.json` |
| SHD cross-check and timing | `python bench/anchor.py target/release/mcs SHD OUT 15 3 data/a/*.json` |
| SHD family, lineage B | `python bench/shd_family.py target/release/mcs SHD data/b/ecc-ethanol.json OUT 300` |
| Census tables | `python bench/census.py target/release/mcs OUT 200 data/a/*-growth.json` |

SHD (Uno) is not bundled; its source is in `MHSGenerationAlgorithms` under `containers/shd`, with
its own academic-use notice. It builds with MSVC `/std:clatest` or a C compiler accepting `typeof`.

## Project map

| Path | Holds |
|---|---|
| `src/hg.rs` | Input contract, normalization, static twins, expansion |
| `src/srch.rs` | MMCS, RS, census |
| `src/chk.rs` | Independent checker and family comparison |
| `src/main.rs` | `mcs` command line |
| `src/bin/efm.rs` | Elementary modes by nullspace double description in exact integers |
| `adapters/` | BiGG model to cone, mode check, modes to hypergraph, LP cut oracle |
| `bench/` | SHD comparison and timing, census tables |
| `scripts/prepare.py` | Pinned downloads with SHA-256 checks and input construction |
| `results/` | Frozen v0.1.0 measurements and oracle families |

## Limits

- Two reconstructions of one organism; no claim across organisms.
- Superset-edge removal is quadratic in the worst case: 18 of 22 s on the 83,601-edge B growth input.
- RS lacks violator pruning and is slow. Checkpoint and resume are not implemented.
- Times are laptop measurements including process start.

## Licence and data

MIT for the code in this repository. Input data are downloaded from their sources by
`scripts/prepare.py` and not redistributed: the METATOOL mode files from
`VeraLiconaResearchGroup/MHSGenerationAlgorithms` and the model from BiGG. Cite the original
authors when using them.
