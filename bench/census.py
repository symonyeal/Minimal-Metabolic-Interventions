"""G2 census table: repeated and dead-end MMCS subtrees per instance and configuration.

Legend:
  ord    candidate order (idx | unhit)      tw   static twins on/off
  dead   share of nodes whose subtree emits nothing
  ex     share of nodes inside subtrees rooted at an exact repeat of a completed state
  wl     same with color-refinement equality: an upper bound on isomorphic repeats
  gap    wl - ex: the most isomorphism can add over exact memoization

Usage: python bench/census.py MCS_EXE OUT_DIR LIMIT_S INPUT.json...
Writes OUT_DIR/census.json and OUT_DIR/census.md. Census shares are deterministic; times are not.
"""
import json, subprocess, sys
from pathlib import Path

mcs, out, lim = str(Path(sys.argv[1]).resolve()), Path(sys.argv[2]), sys.argv[3]
out.mkdir(parents=True, exist_ok=True)
rows = []
for inp in sys.argv[4:]:
    for ord_ in ("idx", "unhit"):
        for tw in (False, True):
            cmd = [mcs, "enumerate", inp, "--count", "--census", "--ord", ord_, "--time", lim] + (["--twins"] if tw else [])
            d = json.loads(subprocess.run(cmd, capture_output=True, text=True, check=True).stdout)
            c = d["census"]
            rows.append({"instance": Path(inp).stem, "ord": ord_, "twins": tw, "status": d["search_status"],
                         "nodes": c["nodes"], "solutions": d["counts"]["expanded_solutions"],
                         "dead": c["dead_share"], "ex": c["exact"]["skippable_share"],
                         "ex_dead_nodes": c["exact"]["skippable_dead_nodes"], "wl": c["wl_upper_bound"]["skippable_share"],
                         "search_s": d["elapsed_s"]["search"]})
            print(json.dumps(rows[-1]), flush=True)
(out / "census.json").write_text(json.dumps(rows, indent=1) + "\n")
p = lambda x: f"{100 * x:.1f}%"
md = ["| instance | order | twins | status | nodes | dead | exact | WL bound | gap |", "|---|---|---|---|---:|---:|---:|---:|---:|"]
for r in rows:
    md.append(f"| {r['instance']} | {r['ord']} | {'yes' if r['twins'] else 'no'} | {r['status']} | {r['nodes']:,} | "
              f"{p(r['dead'])} | {p(r['ex'])} | {p(r['wl'])} | {p(r['wl'] - r['ex'])} |")
(out / "census.md").write_text("\n".join(md) + "\n")
