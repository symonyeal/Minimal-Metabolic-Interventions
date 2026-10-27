"""BiGG JSON model -> integer flux cone for elementary-mode enumeration.

Legend:
  N     stoichiometric matrix (metabolites x reactions)
  sg    bound-sign interval per reaction: (lo, hi) in {-1, 0} x {0, 1}
  blk   reactions zero in every cone vector (blocked), removed
  col   cone columns: (reaction id, direction); every column is irreversible (x >= 0)
  A     integer matrix: each column of N scaled by the lcm of its coefficient denominators

The cone is {v : N v = 0, sign constraints from the bounds}; finite bound magnitudes are ignored,
as in elementary-mode analysis. Reversible reactions are split into forward and reverse columns;
backward-only reactions become one reverse column. Splitting adds one two-cycle per reversible
reaction, which the enumerator drops.

Usage: python adapters/bigg_cone.py MODEL.json OUT.json
"""
import hashlib, json, math, sys
from pathlib import Path
from fractions import Fraction
import numpy as np
from scipy.optimize import linprog

src = sys.argv[1]
raw = open(src, "rb").read()
m = json.loads(raw)
mets = [x["id"] for x in m["metabolites"]]
mi = {k: i for i, k in enumerate(mets)}
rx = m["reactions"]
n = len(rx)
N = np.zeros((len(mets), n))
Q = [dict() for _ in range(n)]  # exact coefficients from their decimal text
for j, r in enumerate(rx):
    for k, c in r["metabolites"].items():
        N[mi[k], j] = c
        Q[j][mi[k]] = Fraction(repr(c))
sg = [(-1 if r["lower_bound"] < 0 else 0, 1 if r["upper_bound"] > 0 else 0) for r in rx]

blk = []
for j in range(n):
    if sg[j] == (0, 0):
        blk.append(j)
        continue
    live = False
    for s in (1, -1):
        if (s == 1 and sg[j][1] == 0) or (s == -1 and sg[j][0] == 0):
            continue
        c = np.zeros(n)
        c[j] = -s
        r = linprog(c, A_eq=N, b_eq=np.zeros(len(mets)), bounds=sg, method="highs")
        if r.status == 0 and -r.fun > 1e-9:
            live = True
            break
    if not live:
        blk.append(j)
bs = set(blk)

col = []
for j in range(n):
    if j in bs:
        continue
    if sg[j][1] == 1:
        col.append((j, 1))
    if sg[j][0] == -1:
        col.append((j, -1))

# Integer columns: scale each by the lcm of its denominators (supports are unchanged).
cols = []
for j, d in col:
    l = 1
    for q in Q[j].values():
        l = l * q.denominator // math.gcd(l, q.denominator)
    cols.append({i: int(q * l) * d for i, q in Q[j].items()})
rows = sorted({i for c in cols for i in c})
A = [[c.get(i, 0) for c in cols] for i in rows]
out = {
    "source": Path(src).name, "model_sha256": hashlib.sha256(raw).hexdigest(), "model_id": m.get("id"),
    "columns": [{"reaction": rx[j]["id"], "dir": d} for j, d in col],
    "metabolites": [mets[i] for i in rows],
    "A": A,
    "blocked": [rx[j]["id"] for j in blk],
    "note": "cone from bound signs only; reversible reactions split into dir 1 and dir -1 columns",
}
json.dump(out, open(sys.argv[2], "w"), indent=None)
print(json.dumps({"reactions": n, "blocked": len(blk), "columns": len(col),
                  "reversible_split": sum(1 for j, d in col if d == -1 and sg[j][1] == 1),
                  "rows": len(rows), "max_abs_coeff": max(abs(x) for r in A for x in r)}))
