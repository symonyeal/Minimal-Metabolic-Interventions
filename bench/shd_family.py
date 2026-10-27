"""Full-family cross-check of mcs against SHD on an mcs-hg/1 instance.

Legend:
  el    eligible ids (all vertices when absent)      ix   id -> 1-based SHD integer
  fa    mcs family (twins + unhit order)              fb   SHD family mapped back to ids

Eligibility is applied before writing the SHD edge list, as mcs does. An empty edge makes the
family empty for both; this script refuses such inputs rather than feeding SHD an empty line.

Usage: python bench/shd_family.py MCS_EXE SHD_EXE HG.json OUT_DIR LIMIT_S
"""
import json, subprocess, sys, time
from pathlib import Path

mcs, shd = str(Path(sys.argv[1]).resolve()), str(Path(sys.argv[2]).resolve())
hg, out, lim = Path(sys.argv[3]), Path(sys.argv[4]), float(sys.argv[5])
out.mkdir(parents=True, exist_ok=True)
d = json.loads(hg.read_text())
el = set(d["eligible"]) if d.get("eligible") is not None else {v["id"] for v in d["vertices"]}
E = [sorted(set(e) & el) for e in d["edges"]]
assert all(E), "empty edge after eligibility"
ids = sorted({x for e in E for x in e})
ix = {x: i + 1 for i, x in enumerate(ids)}
dat = out / f"{hg.stem}.dat"
dat.write_text("".join(" ".join(str(ix[x]) for x in e) + "\n" for e in E))
p, q = out / f"{hg.stem}.mcs.jsonl", out / f"{hg.stem}.shd.out"
t = time.perf_counter()
subprocess.run([mcs, "enumerate", str(hg), "--ord", "unhit", "--twins", "--out", str(p), "--time", str(lim),
                "--summary", str(out / f"{hg.stem}.mcs.summary.json")], check=True, timeout=lim + 60)
ta = time.perf_counter() - t
s = json.loads((out / f"{hg.stem}.mcs.summary.json").read_text())
t = time.perf_counter()
x = subprocess.run([shd, "0", str(dat), str(q)], capture_output=True, text=True, timeout=lim)
tb = time.perf_counter() - t
assert x.returncode == 0, x.stderr
fa = {frozenset(json.loads(l)) for l in p.read_text().splitlines() if l.strip()}
fb = {frozenset(ids[int(v) - 1] for v in l.split()) for l in q.read_text().splitlines()}
r = {"instance": hg.stem, "mcs_status": s["search_status"], "mcs": len(fa), "shd": len(fb),
     "missing": len(fb - fa), "extra": len(fa - fb), "exact_match": fa == fb and s["search_status"] == "complete",
     "mcs_s_with_output": round(ta, 3), "shd_s_with_output": round(tb, 3)}
p.unlink(); q.unlink()
(out / f"{hg.stem}.shd-family.json").write_text(json.dumps(r, indent=1) + "\n")
print(json.dumps(r))
