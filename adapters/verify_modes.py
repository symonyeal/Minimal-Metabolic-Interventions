"""Independent check of elementary modes against the original BiGG model.

Legend:
  N      stoichiometric matrix of the model (float)     Q   exact coefficients (Fraction)
  S      a mode's reaction support                       v   null vector of N[:, S]
  sg     bound-sign interval per reaction
  marg   singular-value margin: smallest kept / largest dropped

A support S is an elementary mode iff N[:, S] has nullity exactly 1 and its null vector is
nonzero on all of S and sign-feasible. Every mode is tested in floating point with margins;
a random sample is re-tested in exact rational arithmetic.

Usage: python adapters/verify_modes.py MODEL.json MODES.jsonl SAMPLE OUT.json
"""
import json, random, sys
from fractions import Fraction
import numpy as np

m = json.load(open(sys.argv[1]))
mets = [x["id"] for x in m["metabolites"]]
mi = {k: i for i, k in enumerate(mets)}
rx = m["reactions"]
ri = {r["id"]: j for j, r in enumerate(rx)}
N = np.zeros((len(mets), len(rx)))
Q = [dict() for _ in rx]
for j, r in enumerate(rx):
    for k, c in r["metabolites"].items():
        N[mi[k], j] = c
        Q[j][mi[k]] = Fraction(repr(c))
sg = [(-1 if r["lower_bound"] < 0 else 0, 1 if r["upper_bound"] > 0 else 0) for r in rx]

modes = [json.loads(l) for l in open(sys.argv[2])]
sups = [tuple(sorted(md["support"])) for md in modes]
bad, worst = [], (np.inf, 0.0)
for t, S in enumerate(sups):
    J = [ri[x] for x in S]
    A = N[:, J]
    u, s, vt = np.linalg.svd(A)
    s = np.concatenate([s, np.zeros(max(0, len(J) - len(s)))])
    tol = 1e-9 * max(1.0, s[0])
    r = int((s > tol).sum())
    if len(J) - r != 1:
        bad.append((t, "nullity", len(J) - r))
        continue
    v = vt[-1]
    worst = (min(worst[0], s[r - 1]), max(worst[1], s[r] if r < len(s) else 0.0))
    if np.abs(v).min() <= 1e-9 * np.abs(v).max():
        bad.append((t, "zero entry in support", None))
        continue
    ok = any(all((sg[j][0] == -1 or x * d > 0) and (sg[j][1] == 1 or x * d < 0) for j, x in zip(J, v)) for d in (1, -1))
    if not ok:
        bad.append((t, "sign infeasible", None))


def rank(rows):
    rows = [r[:] for r in rows]
    k = 0
    for c in range(len(rows[0]) if rows else 0):
        p = next((i for i in range(k, len(rows)) if rows[i][c] != 0), None)
        if p is None:
            continue
        rows[k], rows[p] = rows[p], rows[k]
        for i in range(len(rows)):
            if i != k and rows[i][c] != 0:
                f = rows[i][c] / rows[k][c]
                rows[i] = [a - f * b for a, b in zip(rows[i], rows[k])]
        k += 1
    return k


random.seed(20261027)
smp = random.sample(range(len(sups)), min(int(sys.argv[3]), len(sups)))
xbad = []
for t in smp:
    J = [ri[x] for x in sups[t]]
    rows = [[Q[j].get(i, Fraction(0)) for j in J] for i in range(len(mets))]
    if len(J) - rank(rows) != 1:
        xbad.append(t)
out = {
    "modes": len(sups), "distinct_supports": len(set(sups)), "float_failures": len(bad),
    "first_failures": bad[:10], "min_kept_singular_value": worst[0], "max_dropped_singular_value": worst[1],
    "exact_sample": len(smp), "exact_failures": len(xbad), "sample_seed": 20261027,
}
json.dump(out, open(sys.argv[4], "w"), indent=1)
print(json.dumps(out))
