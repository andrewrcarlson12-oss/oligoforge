"""RDML 1.3 export gate. Offline. Run from repo root:  python tests/test_rdml.py"""
import base64, io, sys, zipfile
import xml.etree.ElementTree as ET
sys.path.insert(0, ".")
from oligoforge import rdml

NS = "{http://www.rdml.org}"
_fails = []
def check(name, cond, detail=""):
    print(("  PASS " if cond else "  FAIL ") + name + (f"  [{detail}]" if detail and not cond else ""))
    if not cond:
        _fails.append(name)

# a representative panel: two reference genes (RPL13, YWHAZ), a target of interest with validation,
# and a SYBR (no-probe) assay
panel = [
    dict(name="FSJ IFNG", gene="IFNG", organism="Aphelocoma coerulescens",
         forward="AGTCATTCTGATGTCGCTGATG", reverse="ACCTGTCAGTGTTTTCAAGCA",
         probe="TCATTTCTCTCTGTCCAGCCTGATAGCTTCTCT", amplicon=136, dye="FAM",
         chemistry="IDT PrimeTime (ZEN/FAM)", validation=dict(efficiency_pct=99.2, r2=0.998)),
    dict(name="FSJ RPL13", gene="RPL13", organism="Aphelocoma coerulescens", target_type="ref",
         forward="TCGCTGGCATCAACAAGAAG", reverse="TCGGGAAGAGGATGAGCTTG",
         probe="AACAAGTCCACCGAGTCCCTGCA", amplicon=138, dye="FAM", chemistry="ZEN/FAM"),
    dict(name="FSJ YWHAZ", gene="YWHAZ", target_type="ref", forward="CCGTTACTTGGCTGAGGTTG",
         reverse="GATGGGATGTGTTGGTTGCA", probe="CCACTATCCCTTTCTTGTCATCTCCAGCAG", dye="FAM"),
    dict(name="Plas cytb", gene="Plasmodium cytb", forward="TACCTGGACTWGTTTCATGG",
         reverse="AAAGGATTTGTGCTACCTTG", chemistry="low-Tm SYBR"),   # no probe -> SYBR dye
]

out = rdml.build(panel, meta={"title": "OligoForge FSJ panel"})
check("build returns xml + b64 + filename", bool(out.get("xml") and out.get("rdml_b64")), list(out))
check("n_assays == 4", out["n_assays"] == 4, out["n_assays"])
check("filename ends .rdml", out["filename"].endswith(".rdml"), out["filename"])

# 1) well-formed XML
root = None
try:
    root = ET.fromstring(out["xml"])
    check("XML is well-formed", True)
except ET.ParseError as e:
    check("XML is well-formed", False, str(e))

if root is not None:
    check("root is rdml v1.3", root.tag == NS + "rdml" and root.get("version") == "1.3",
          (root.tag, root.get("version")))
    targets = root.findall(NS + "target")
    dyes = root.findall(NS + "dye")
    check("one <target> per assay", len(targets) == 4, len(targets))
    check("dyes declared before targets (order)",
          list(root).index(dyes[0]) < list(root).index(targets[0]) if dyes and targets else False)
    check("FAM + SYBR dyes present", {d.get("id") for d in dyes} >= {"FAM", "SYBR"},
          [d.get("id") for d in dyes])

    by = {t.get("id"): t for t in targets}
    ifng = by.get("FSJ_IFNG"); rpl13 = by.get("FSJ_RPL13"); ywhaz = by.get("FSJ_YWHAZ"); plas = by.get("Plas_cytb")
    check("all target ids resolved", all([ifng is not None, rpl13 is not None, ywhaz is not None, plas is not None]),
          list(by))

    def ttype(t): return (t.find(NS + "type").text if t is not None and t.find(NS + "type") is not None else None)
    check("reference gene RPL13 typed 'ref'", ttype(rpl13) == "ref", ttype(rpl13))
    check("reference gene YWHAZ typed 'ref'", ttype(ywhaz) == "ref", ttype(ywhaz))
    check("IFNG typed 'toi'", ttype(ifng) == "toi", ttype(ifng))

    # 2) sequences present and correct
    if ifng is not None:
        seqs = ifng.find(NS + "sequences")
        fwd = seqs.find(NS + "forwardPrimer/" + NS + "sequence") if seqs is not None else None
        prb = seqs.find(NS + "probe1/" + NS + "sequence") if seqs is not None else None
        check("IFNG forward primer sequence preserved",
              fwd is not None and fwd.text == "AGTCATTCTGATGTCGCTGATG", fwd.text if fwd is not None else None)
        check("IFNG probe1 present", prb is not None and prb.text.startswith("TCATTT"), prb.text if prb is not None else None)

    # 3) efficiency conversion: 99.2% -> E = 1.992
    eff = ifng.find(NS + "amplificationEfficiency") if ifng is not None else None
    check("efficiency converted to RDML E (99.2%% -> 1.992)",
          eff is not None and abs(float(eff.text) - 1.992) < 1e-6, eff.text if eff is not None else None)

    # 4) SYBR assay: dyeId SYBR, no probe1
    if plas is not None:
        did = plas.find(NS + "dyeId")
        seqs = plas.find(NS + "sequences")
        prb = seqs.find(NS + "probe1") if seqs is not None else None
        check("SYBR assay dyeId == SYBR", did is not None and did.get("id") == "SYBR", did.get("id") if did is not None else None)
        check("SYBR assay has no probe", prb is None)

# 5) the base64 .rdml decodes to a zip whose XML round-trips
try:
    raw = base64.b64decode(out["rdml_b64"])
    z = zipfile.ZipFile(io.BytesIO(raw))
    inner = z.read(z.namelist()[0]).decode("utf-8")
    ET.fromstring(inner)
    check("rdml_b64 is a valid zip containing well-formed XML", True)
    check("zip XML matches returned XML", inner.strip() == out["xml"].strip())
except Exception as e:
    check("rdml_b64 is a valid zip containing well-formed XML", False, repr(e))

print("")
if _fails:
    print("RDML GATE FAILED:", ", ".join(_fails)); sys.exit(1)
print("ALL RDML ASSERTS PASS")
