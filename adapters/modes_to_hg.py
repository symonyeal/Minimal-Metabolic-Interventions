"""Elementary-mode supports -> mcs-hg/1 hypergraph for one scenario.

Legend:
  T      target reaction (or "all" for every mode)       el  eligible interventions
  pol    eligibility policy: "all-but-target" (S2 analogue) or "core" (S3 policy: no exchange,
         demand, sink, biomass, ATP maintenance, or reaction whose bounds exclude zero)

Edges are the supports of modes that use T, intersected with the eligible set by the tool. An
empty edge (a target mode with no eligible reaction) is kept: it means no solution exists.

Usage: python adapters/modes_to_hg.py MODEL.json MODES.jsonl TARGET|all POLICY OUT.json
"""
import hashlib, json, sys

mp, mo, T, pol, out = sys.argv[1:6]
m = json.load(open(mp))
rx = {r["id"]: r for r in m["reactions"]}
modes = [json.loads(l)["support"] for l in open(mo)]
E = [s for s in modes if T == "all" or T in s]
V = sorted({x for s in E for x in s})


def core_ok(r):
    i = r["id"]
    if i.startswith(("EX_", "DM_", "SK_")) or i.startswith("BIOMASS") or i == "ATPM":
        return False
    return r["lower_bound"] <= 0 <= r["upper_bound"]


if pol == "all-but-target":
    el = [x for x in V if x != T]
elif pol == "core":
    el = [x for x in V if core_ok(rx[x])]
else:
    sys.exit(f"unknown policy {pol}")
doc = {
    "schema": "mcs-hg/1",
    "vertices": [{"id": x, "label": rx[x].get("name") or None} for x in V],
    "edges": E,
    "eligible": el,
    "meta": {"model": mp, "model_sha256": hashlib.sha256(open(mp, "rb").read()).hexdigest(),
             "modes": mo, "modes_sha256": hashlib.sha256(open(mo, "rb").read()).hexdigest(),
             "target": T, "policy": pol, "target_modes": len(E), "all_modes": len(modes)},
}
json.dump(doc, open(out, "w"))
print(json.dumps({"target": T, "policy": pol, "edges": len(E), "vertices": len(V), "eligible": len(el)}))
