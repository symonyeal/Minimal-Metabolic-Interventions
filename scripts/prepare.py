"""Fetch and build every input used in EVIDENCE.md and REPORT.md.

Legend:
  A     lineage A: METATOOL E. coli elementary-mode supports from MHSGenerationAlgorithms @ REV
  B     lineage B: BiGG e_coli_core and its elementary modes computed by efm
  sha   expected SHA-256 of each downloaded file; a mismatch stops the run
  grow  A growth instance: modes containing mue (id 1), with mue removed

Writes data/ (not tracked: third-party redistribution terms are unresolved).
Usage: python scripts/prepare.py [--skip-b]
"""
import hashlib, json, subprocess, sys
from pathlib import Path
from urllib.request import urlopen

root = Path(__file__).resolve().parent.parent
dat = root / "data"
REPO, REV = "VeraLiconaResearchGroup/MHSGenerationAlgorithms", "f043cc3fae8088191307d59b987ce5fe9de4a326"
A = {
    "acetate": "bf5f50c582c8631af024dfcfe574f6892721bb72887c039e275f0ebf0b9b5717",
    "succinate": "8ad10e583218aa19b91530c6eb26d4f417b1fa6e7f9b6f3f2b639511ecea44fc",
    "glycerol": "a0b1afc2eec0852df9fbd5b98a57b2bbef25d78bbc00ae50a31d080ddadac413",
    "glucose": "a666e621bd597c6a8ce97931c0ebb49a9464063459267d0741be436d2284ef89",
    "combined": "141f94db708f88ff12781f5aac213ffe9aa7718c3a84dd03179428f7194a27db",
}
B_URL = "https://bigg.ucsd.edu/static/models/e_coli_core.json"
B_SHA = "7bedec10576cfe935b19218dc881f3fb14f890a1871448fc19a9b4ee15b448d8"


def get(u, sha, out):
    if out.exists() and hashlib.sha256(out.read_bytes()).hexdigest() == sha:
        return
    with urlopen(u, timeout=60) as x:
        b = x.read()
    h = hashlib.sha256(b).hexdigest()
    if h != sha:
        sys.exit(f"{u}: sha256 {h}, expected {sha}")
    out.write_bytes(b)


def run(*a):
    print(" ".join(map(str, a)), flush=True)
    subprocess.run([str(x) for x in a], check=True)


(dat / "a").mkdir(parents=True, exist_ok=True)
for n, sha in A.items():
    # The JSON files are stored with Git LFS; media.githubusercontent.com serves their content.
    p = dat / "a" / f"ecoli-{n}.json"
    get(f"https://media.githubusercontent.com/media/{REPO}/{REV}/benchmark_inputs/ecoli/ecoli-{n}.json", sha, p)
    s = json.loads(p.read_text())["sets"]
    g = [[x for x in e if x != 1] for e in s if 1 in e]
    (dat / "a" / f"ecoli-{n}-growth.json").write_text(json.dumps({
        "sets": g, "scenario": "growth blocking: modes containing mue (id 1), mue removed", "source": f"{REPO}@{REV}"}))
    print(n, "modes", len(s), "growth modes", len(g))

if "--skip-b" in sys.argv:
    sys.exit(0)
b = dat / "b"
b.mkdir(exist_ok=True)
mdl = b / "e_coli_core.json"
get(B_URL, B_SHA, mdl)
exe = root / "target" / "release" / ("efm.exe" if sys.platform == "win32" else "efm")
py = sys.executable
run(py, root / "adapters/bigg_cone.py", mdl, b / "cone.json")
run(exe, b / "cone.json", b / "modes.jsonl", b / "efm-summary.json")
for t, pol, nm in [("BIOMASS_Ecoli_core_w_GAM", "all-but-target", "growth"), ("EX_etoh_e", "core", "ethanol"), ("all", "all-but-target", "all")]:
    run(py, root / "adapters/modes_to_hg.py", mdl, b / "modes.jsonl", t, pol, b / f"ecc-{nm}.json")
