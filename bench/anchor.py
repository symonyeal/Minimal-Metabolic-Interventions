"""SHD cross-check and competitiveness anchor (laptop smoke scale).

Legend:
  mcs    Rust binary              shd   original SHD binary (Uno), MSVC /std:clatest build
  inp    instance JSON (sets)      dat   SHD edge-list input written from inp
  lim    per-process time limit    rep   timed repetitions after one warmup
  fam    solution family as a set of frozensets of original ids

Usage: python bench/anchor.py MCS_EXE SHD_EXE OUT_DIR LIMIT_S REPS INPUT.json...
Writes OUT_DIR/anchor.json. Times include process start; they are laptop smoke measurements,
not benchmark-runner results.
"""
import json, subprocess, sys, time
from pathlib import Path

mcs, shd, out, lim, rep = sys.argv[1], sys.argv[2], Path(sys.argv[3]), float(sys.argv[4]), int(sys.argv[5])
mcs, shd = str(Path(mcs).resolve()), str(Path(shd).resolve())
out.mkdir(parents=True, exist_ok=True)


def timed(cmd):
    t = time.perf_counter()
    try:
        x = subprocess.run(cmd, capture_output=True, text=True, timeout=lim)
    except subprocess.TimeoutExpired:
        return None, None
    return time.perf_counter() - t, x


def fam_lines(p):
    return {frozenset(map(str, json.loads(l))) for l in Path(p).read_text().splitlines() if l.strip()}


res = []
for inp in sys.argv[6:]:
    inp = Path(inp)
    nm = inp.stem
    sets = json.loads(inp.read_text())["sets"]
    dat = out / f"{nm}.dat"
    dat.write_text("".join(" ".join(map(str, sorted(set(e)))) + "\n" for e in sets))
    r = {"instance": nm, "edges": len(sets), "runs": {}}
    cfgs = {
        "mcs_mmcs_idx": [mcs, "enumerate", str(inp), "--count", "--time", str(lim)],
        "mcs_mmcs_unhit": [mcs, "enumerate", str(inp), "--count", "--ord", "unhit", "--time", str(lim)],
        "mcs_mmcs_unhit_twins": [mcs, "enumerate", str(inp), "--count", "--ord", "unhit", "--twins", "--time", str(lim)],
        "mcs_rs": [mcs, "enumerate", str(inp), "--count", "--alg", "rs", "--time", str(lim)],
        "shd_0": [shd, "0", str(dat)],
        "shd_D": [shd, "D", str(dat)],
    }
    for k, cmd in cfgs.items():
        ts, cnt, st = [], None, "complete"
        for i in range(rep + 1):
            t, x = timed(cmd)
            if t is None:
                st = "timeout"
                break
            if k.startswith("mcs"):
                s = json.loads(x.stdout)
                if s["search_status"] != "complete":
                    st = s["search_status"]
                    break
                cnt = int(s["counts"]["expanded_solutions"])
            else:
                assert x.returncode == 0, x.stderr
                cnt = int(x.stdout.split()[0])
            if i:
                ts.append(t)
        r["runs"][k] = {"status": st, "solutions": cnt, "times_s": ts,
                        "median_s": sorted(ts)[len(ts) // 2] if ts else None}
    # Full-family comparison of Rust MMCS and SHD when both finished.
    a, b = r["runs"]["mcs_mmcs_idx"], r["runs"]["shd_0"]
    if a["status"] == b["status"] == "complete" and a["solutions"] <= 3_000_000:
        p, q = out / f"{nm}.mcs.jsonl", out / f"{nm}.shd.out"
        subprocess.run([mcs, "enumerate", str(inp), "--ord", "unhit", "--twins", "--out", str(p),
                        "--summary", str(out / f"{nm}.mcs.summary.json")], check=True, timeout=10 * lim)
        subprocess.run([shd, "0", str(dat), str(q)], check=True, capture_output=True, timeout=10 * lim)
        fa = fam_lines(p)
        fb = {frozenset(l.split()) for l in q.read_text().splitlines()}
        r["family_vs_shd"] = {"mcs": len(fa), "shd": len(fb), "missing": len(fb - fa), "extra": len(fa - fb),
                              "exact_match": fa == fb}
        p.unlink(); q.unlink()
    res.append(r)
    print(json.dumps(r), flush=True)
(out / "anchor.json").write_text(json.dumps({"limit_s": lim, "reps": rep, "results": res}, indent=1) + "\n")
