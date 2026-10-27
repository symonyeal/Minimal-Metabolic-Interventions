# Minimal metabolic interventions

A Rust tool that enumerates inclusion-minimal hitting sets of an explicit hypergraph, used here
for minimal reaction cut sets of metabolic networks given their elementary-mode supports. It is
the software half of the study specified in the PRD
(`Useful Misc/Claude Func Folder/hypergraph-project-evaluation/PRD.md`).

Status: G1 baseline and G2 census implemented and checked on a laptop. Not published; the
repository name, licence (proposed GPL-3.0-or-later) and release are owner decisions.

## Build

```text
cargo build --release
cargo test --release
```

## Commands

```text
mcs inspect FILE
mcs enumerate FILE [--alg mmcs|rs] [--ord idx|unhit] [--twins] [--k K] [--count]
              [--out SOLS] [--summary JSON] [--census] [--time SECS] [--cap N] [--check]
mcs check FILE SOLS [--ref SOLS]
efm CONE.json MODES.jsonl SUMMARY.json
```

- `--alg mmcs` is Murakami–Uno MMCS with the min-intersection branching edge and violator
  pruning; `--alg rs` is their RS (without violator pruning).
- `--ord unhit` processes branching candidates by increasing unhit degree (Hörmann–Schirneck 2026);
  this is a reimplementation, not the authors' code.
- `--twins` merges vertices with identical incidence and expands each solution over its classes.
- `--k` keeps the inclusion-minimal hitting sets of size at most K (status `bounded_complete`).
- `--check` verifies every expanded solution against the raw input rows: coverage and a private
  witness row for each selected vertex. It shares no code with the search.
- `--census` records, without pruning, how much of the MMCS tree lies in subtrees that are dead
  ends or repeat an earlier completed state (exactly, or up to colour-refinement equality).
- `efm` enumerates elementary modes of an integer cone written by `adapters/bigg_cone.py`.

## Input

`mcs-hg/1` JSON:

```json
{"schema": "mcs-hg/1",
 "vertices": [{"id": "PFK", "label": "Phosphofructokinase"}],
 "edges": [["PFK", "PGI"]],
 "eligible": ["PFK"],
 "meta": {}}
```

Publisher `{"sets": [[1, 2], ...], "labels": {"1": "mue"}}` files are read directly. Unknown ids,
conflicting labels, malformed numbers and other schema versions are rejected. Edges are intersected
with `eligible`. An empty edge means no solution exists; no edges means the empty set is the only
solution. Duplicate and superset edges are dropped with their provenance kept.

## Output

`--out` writes one JSON array of original ids per line. The summary carries `search_status`
(`complete`, `bounded_complete`, `resource_limited`, `invalid_input`), `verification_status`
(`not_run`, `passed`, `failed`), input hashes (bytes and normalized form), method settings,
representative and expanded counts, node count, stage times and the census. It is written last
through a rename, so a missing summary means an unfinished run. Peak memory is the process peak
working set on Windows and is not measured on other platforms yet.

## Layout

| Path | Holds |
|---|---|
| `src/bits.rs` | Dynamic bitset |
| `src/hg.rs` | Input contract, normalization, static twins, expansion |
| `src/srch.rs` | MMCS, RS, census |
| `src/chk.rs` | Independent solution checker and family comparison |
| `src/main.rs` | `mcs` CLI |
| `src/bin/efm.rs` | Elementary modes by nullspace double description (exact integers) |
| `adapters/` | BiGG model to cone, mode verification, modes to hypergraph, LP cut oracle |
| `bench/` | SHD cross-check and timing anchor, census table |
| `tests/fam.rs` | Brute-force family tests, degenerate inputs, checker and input-contract tests |

Input data are not bundled; see `EVIDENCE.md` for sources, hashes and the commands behind every
reported number.
