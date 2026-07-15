import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from oligoforge import expression as EXP

fails = []
def chk(name, cond, detail=""):
    print(("  PASS " if cond else "  FAIL ") + name + (f"  [{detail}]" if detail and not cond else ""))
    if not cond:
        fails.append(name)
def approx(a, b, t=1e-6):
    return a is not None and abs(a - b) <= t

# 1) Livak basic + multi-reference + control == 1
samples = [
    {"sample": "c1", "group": "control",  "cq": {"IFNG": 25.0, "RPL13": 20.0, "YWHAZ": 20.0}},
    {"sample": "c2", "group": "control",  "cq": {"IFNG": 25.0, "RPL13": 20.0, "YWHAZ": 20.0}},
    {"sample": "i1", "group": "infected", "cq": {"IFNG": 23.0, "RPL13": 20.0, "YWHAZ": 20.0}},
    {"sample": "i2", "group": "infected", "cq": {"IFNG": 23.0, "RPL13": 20.0, "YWHAZ": 20.0}},
]
r = EXP.analyze(samples, ["RPL13", "YWHAZ"], "control")
row = next(x for x in r["results"] if x["target"] == "IFNG" and x["group"] == "infected")
chk("IFNG infected ddCq = -2.0", approx(row["ddcq"], -2.0), row["ddcq"])
chk("IFNG infected fold (Livak) = 4.0", approx(row["fold_livak"], 4.0, 1e-3), row["fold_livak"])
chk("IFNG infected log2 fold = 2.0", approx(row["log2_fold"], 2.0, 1e-3), row["log2_fold"])
ctrl = next(x for x in r["results"] if x["target"] == "IFNG" and x["group"] == "control")
chk("control group fold = 1.0", approx(ctrl["fold_livak"], 1.0), ctrl["fold_livak"])
chk("targets exclude reference genes", r["targets"] == ["IFNG"], str(r["targets"]))

# Pfaffl reduces to Livak when all efficiencies are 100%
r2 = EXP.analyze(samples, ["RPL13", "YWHAZ"], "control", {"IFNG": 100, "RPL13": 100, "YWHAZ": 100})
row2 = next(x for x in r2["results"] if x["target"] == "IFNG" and x["group"] == "infected")
chk("Pfaffl == Livak at 100% efficiency",
    approx(row2.get("fold_pfaffl"), row2["fold_livak"], 1e-3),
    f'{row2.get("fold_pfaffl")} vs {row2["fold_livak"]}')

# 2) Pfaffl diverges from Livak with a low-efficiency target (E_target = 1.8 -> ratio 1.8^2 = 3.24)
s2 = [
    {"sample": "c1", "group": "control",  "cq": {"IFNG": 25.0, "RPL13": 20.0}},
    {"sample": "i1", "group": "infected", "cq": {"IFNG": 23.0, "RPL13": 20.0}},
]
r3 = EXP.analyze(s2, ["RPL13"], "control", {"IFNG": 80, "RPL13": 100})
row3 = next(x for x in r3["results"] if x["target"] == "IFNG" and x["group"] == "infected")
chk("Livak fold still 4.0 (efficiency-independent)", approx(row3["fold_livak"], 4.0, 1e-3), row3["fold_livak"])
chk("Pfaffl fold = 3.24 at E_target=1.8", approx(row3.get("fold_pfaffl"), 3.24, 1e-2), row3.get("fold_pfaffl"))

# 3) parse_table (CSV)
txt = "sample,group,IFNG,RPL13\nc1,control,25,20\ni1,infected,23,20\n"
sp = EXP.parse_table(txt)
chk("parse_table parses CSV",
    len(sp) == 2 and sp[0]["cq"]["IFNG"] == 25.0 and sp[1]["group"] == "infected", str(sp))
# tab-delimited + an undetermined cell
txt2 = "sample\tgroup\tIFNG\tRPL13\nc1\tcontrol\tUndetermined\t20\n"
sp2 = EXP.parse_table(txt2)
chk("parse_table drops undetermined", "IFNG" not in sp2[0]["cq"] and sp2[0]["cq"]["RPL13"] == 20.0, str(sp2))

# 4) validation guards
def raises(fn):
    try:
        fn(); return False
    except ValueError:
        return True
    except Exception:
        return False
chk("no reference genes -> error", raises(lambda: EXP.analyze(samples, [], "control")))
chk("bad control group -> error", raises(lambda: EXP.analyze(samples, ["RPL13"], "nope")))
chk("reference not in data -> error", raises(lambda: EXP.analyze(samples, ["FOO"], "control")))
chk("bad header -> error", raises(lambda: EXP.parse_table("a,b\n1,2")))

print("")
if fails:
    print("EXPRESSION GATE FAILED:", ", ".join(fails)); sys.exit(1)
print("ALL EXPRESSION ASSERTS PASS")
