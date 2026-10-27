"""Exhaustive LP oracle: inclusion-minimal cut sets of size <= k blocking a target in the flux cone.

Legend:
  N    stoichiometric matrix         sg   bound-sign interval per reaction (cone semantics)
  T    target reaction               el   eligible reactions (from an mcs-hg/1 file)
  C    candidate cut                 blk  C blocks T: max v_T = 0 over {N v = 0, signs, |v| <= 1, v_C = 0}

Independent of elementary modes: blockage is decided by HiGHS on the model itself. C is a minimal
cut iff it blocks T and no C minus one element does; by monotonicity it suffices to test the
maximal proper subsets.

Usage: python adapters/lp_cuts.py MODEL.json HG.json K OUT.json
"""
import itertools, json, sys
import numpy as np
from scipy.optimize import linprog

m = json.load(open(sys.argv[1]))
hg = json.load(open(sys.argv[2]))
k = int(sys.argv[3])
mets = [x["id"] for x in m["metabolites"]]
mi = {x: i for i, x in enumerate(mets)}
rx = m["reactions"]
ri = {r["id"]: j for j, r in enumerate(rx)}
N = np.zeros((len(mets), len(rx)))
for j, r in enumerate(rx):
    for x, c in r["metabolites"].items():
        N[mi[x], j] = c
sg = [(-1.0 if r["lower_bound"] < 0 else 0.0, 1.0 if r["upper_bound"] > 0 else 0.0) for r in rx]
T = ri[hg["meta"]["target"]]
el = sorted(hg["eligible"])
memo = {}


def blk(C):
    C = frozenset(C)
    if C not in memo:
        b = list(sg)
        for x in C:
            b[ri[x]] = (0.0, 0.0)
        c = np.zeros(len(rx))
        c[T] = -1.0
        r = linprog(c, A_eq=N, b_eq=np.zeros(len(mets)), bounds=b, method="highs")
        assert r.status == 0, (sorted(C), r.message)
        memo[C] = -r.fun <= 1e-9
    return memo[C]


assert not blk(()), "target already blocked"
cuts = []
for s in range(1, k + 1):
    for C in itertools.combinations(el, s):
        if blk(C) and all(not blk(D) for D in itertools.combinations(C, s - 1)):
            cuts.append(sorted(C))
json.dump({"target": hg["meta"]["target"], "k": k, "eligible": len(el), "lps": len(memo), "cuts": cuts},
          open(sys.argv[4], "w"), indent=0)
print(json.dumps({"target": hg["meta"]["target"], "k": k, "eligible": len(el), "lps": len(memo), "cuts": len(cuts)}))
