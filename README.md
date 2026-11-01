# Minimal metabolic interventions

`mcs` computes minimal cut sets of metabolic networks. It enumerates the elementary modes of a
BiGG model in exact integer arithmetic, enumerates the minimal hitting sets of their supports with
MMCS, checks solutions independently and measures where the search spends its nodes. This
repository studies MMCS search on minimal reaction cut sets from two reconstructions of *E. coli*.

## Background

A minimal cut set of a metabolic network is a set of reactions whose removal blocks a target
function, such as growth or ethanol export, and that contains no smaller blocking set. If the
elementary modes carrying the target are known, the minimal cut sets are exactly the minimal
hitting sets of the modes' reaction supports (Haus, Klamt and Stephen 2008). Enumerating them is
hypergraph dualization, for which MMCS (Murakami and Uno 2014) is the standard practical algorithm.

## Result

The census in [REPORT.md](REPORT.md) measured search nodes without pruning. It records nodes
whose subtrees emit no solutions or repeat an earlier completed state. With static twins
(vertices with identical incidence) merged and candidates ordered by increasing unhit degree:

| Instance | Nodes | Dead | Exact repeats | Isomorphic repeats (upper bound) |
|---|---:|---:|---:|---:|
| A acetate-growth | 77 | 26.0% | 1.3% | 5.2% |
| A succinate-growth | 476 | 60.3% | 15.8% | 22.3% |
| A glycerol-growth | 1,162 | 62.8% | 18.2% | 25.3% |
| A glucose-growth | 3,901 | 76.6% | 21.0% | 28.0% |
| A combined-growth | 38,391 | 79.3% | 15.7% | 19.8% |
| B ethanol | 41,948 | 77.8% | 35.9% | 42.7% |

The possible gain from isomorphic-state reuse is at most 7.1 percentage points beyond exact
repeats, below the predeclared 1/6 threshold; the study stopped there. Dead subtrees contain most
search nodes. Lineage A uses the METATOOL network of Haus, Klamt and Stephen on five substrates.
Lineage B uses BiGG `e_coli_core`; `mcs modes` reproduces its 100,274 elementary modes.

Every enumerated family used in these tables was checked against at least one independent
enumerator (brute force, a Berge implementation, SHD or an LP cut oracle); see
[EVIDENCE.md](EVIDENCE.md).

## Usage

```text
cargo build --release
cargo test --release
target/release/mcs prepare data     # download (curl, pinned SHA-256) and build every input
```

```text
mcs enumerate FILE [--ord idx|unhit] [--twins] [--k K] [--count] [--census] [--check]
              [--time SECS] [--out SOLS] [--summary JSON]
mcs check FILE SOLS [--ref SOLS]
mcs census OUT SECS FILE...
mcs modes MODEL MODES
mcs prepare DIR [--skip-b]
```

| Option | Effect |
|---|---|
| `--ord unhit` | Candidates by increasing unhit degree (Hörmann and Schirneck 2026), reimplemented; default `idx` |
| `--twins` | Merge vertices with identical incidence; solutions are expanded over their classes |
| `--k K` | Only minimal hitting sets of at most K vertices (`bounded_complete`) |
| `--count` | Count solutions without writing them; exclusive with `--out` |
| `--census` | Record dead and repeated subtrees (REPORT.md §2) |
| `--check` | Verify each expanded solution against the raw rows: coverage and a private witness row per vertex |
| `--time SECS` | Stop the search after SECS (`resource_limited`) |

`enumerate` reads `mcs-hg/1` JSON, or a publisher file `{"sets": [[1, 2], ...], "labels": {"1": "mue"}}`:

```json
{"schema": "mcs-hg/1",
 "vertices": [{"id": "PFK", "label": "Phosphofructokinase"}],
 "edges": [["PFK", "PGI"]],
 "eligible": ["PFK"],
 "meta": {}}
```

Unknown ids, conflicting labels, malformed numbers, unknown options and other schema versions are
rejected. Edges are intersected with `eligible`. An empty edge means no solution exists; no edges
means the empty set is the only solution. Duplicate and superset edges are removed.

`--out` writes one JSON array of original ids per line. The summary reports `search_status`
(`complete`, `bounded_complete`, `resource_limited`, `invalid_input`), `verification_status`
(`not_run`, `passed`, `failed`), hashes of the input bytes and of the normalized form, settings,
representative and expanded counts, nodes, stage times and the census. It is written last through
a rename; a missing summary means an unfinished run. Exit status is 1 when verification fails, 2 on
a usage or I/O error and 3 on invalid input.

`check` verifies a solutions file (a JSON array of arrays or JSON lines) and, with `--ref`, compares
it with a reference family. `census` runs both candidate orders with and without twins on each file
and writes `census.json` and `census.md`.

`modes` reads a BiGG JSON model and writes one line per elementary mode,
`{"columns": [[reaction, ±1], ...], "support": [...]}`. Its cone is the model's stoichiometry under
the bound signs; reversible reactions are split, and the two-cycles this creates are dropped.
Reactions that sign propagation proves blocked are removed first; the summary lists them as
`pruned`, and as `blocked` the reactions that occur in no mode, which are exactly the blocked ones.

## Reproduction

| Table or claim | Command |
|---|---|
| Tests | `cargo test --release` |
| Inputs, lineage B modes | `mcs prepare data` |
| Counts | `mcs enumerate data/b/ecc-ethanol.json --count --ord unhit --twins` |
| Ethanol cut sets, all checked | `mcs enumerate data/b/ecc-ethanol.json --ord unhit --twins --check --count` |
| LP oracle, cuts of size ≤ 2 | `mcs enumerate data/b/ecc-growth.json --k 2 --out S` then `mcs check data/b/ecc-growth.json S --ref results/oracle/ecc-growth-k2.oracle.json` |
| Berge oracle, acetate-growth | `mcs enumerate data/a/ecoli-acetate-growth.json --out S` then `mcs check data/a/ecoli-acetate-growth.json S --ref results/oracle/acetate-growth.oracle.json` |
| Census tables | `mcs census OUT 240 data/a/ecoli-*-growth.json data/b/ecc-ethanol.json` |

The SHD comparison, the LP cut oracle and the floating-point mode check were external tools run
once for v0.1.0. Their outputs are in `results/`; their drivers are in the `v0.1.0` tag. SHD (Uno)
is in `MHSGenerationAlgorithms` under `containers/shd`, with its own academic-use notice.

## Project map

| Path | Holds |
|---|---|
| `src/bits.rs` | Fixed-length bitsets |
| `src/hg.rs` | Normalization, static twins, expansion |
| `src/srch.rs` | MMCS |
| `src/cen.rs` | Census keys and counters |
| `src/chk.rs` | Independent checker and family comparison |
| `src/dd.rs` | Extreme rays of an integer cone by double description |
| `src/net.rs` | Metabolic network to flux cone, modes and cut-set scenarios |
| `src/bin/mcs/` | Command line: input formats, downloads, census tables |
| `tests/` | Brute-force oracles for MMCS and for modes; the command line against `results/oracle/` |
| `results/` | Frozen measurements and oracle families |

## Limits

- Two reconstructions of one organism; no claim across organisms.
- Superset-edge removal is quadratic in the worst case: about 18 s on the 83,601-edge B growth input.
- Sign propagation can miss blocked reactions; any it misses stay in the cone, which leaves the
  modes unchanged and shows in the summary as `blocked` ≠ `pruned`.
- Checkpoint and resume are not implemented. Times are laptop measurements.

## Licence and data

MIT for the code in this repository. Input data are downloaded from their sources by
`mcs prepare` and not redistributed: the METATOOL mode files from
`VeraLiconaResearchGroup/MHSGenerationAlgorithms` and the model from BiGG. Cite the original
authors when using them.
