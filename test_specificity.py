import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from oligoforge import specificity as SP

fails = []
def chk(name, cond, detail=""):
    print(("  PASS " if cond else "  FAIL ") + name + (f"  [{detail}]" if detail and not cond else ""))
    if not cond:
        fails.append(name)

# epcr: two convergent on-subject hits in range -> one product (composition sanity, no network)
hits = [
    {"primer": "F", "subject": "NM_t", "lo": 100, "hi": 119, "strand": "+", "q3": True},
    {"primer": "R", "subject": "NM_t", "lo": 230, "hi": 249, "strand": "-", "q3": True},
]
prods = SP.epcr(hits, min_product=40, max_product=3000)
chk("epcr predicts one product", len(prods) == 1 and prods[0]["subject"] == "NM_t", str(prods))
chk("epcr product size correct", prods and prods[0]["size"] == 150, str(prods))

# verdict: off-size product present (480 vs modal 150), no probe binding there
prod3 = [
    {"subject": "NM_target", "size": 150, "left": "F", "right": "R", "span": [100, 249]},
    {"subject": "NM_target", "size": 150, "left": "F", "right": "R", "span": [500, 649]},
    {"subject": "XR_off",     "size": 480, "left": "F", "right": "R", "span": [10, 489]},
]
v = SP.assay_verdict(prod3, [{"subject": "NM_target", "lo": 150, "hi": 170}])
chk("modal size = 150", v["modal_size"] == 150, v["modal_size"])
chk("one off-size product flagged", v["n_off_size"] == 1, v["n_off_size"])
chk("no probe cross-binding", v["n_probe_binding"] == 0, v["n_probe_binding"])
p1 = next(p for p in v["products"] if p["span"] == [100, 249])
chk("probe binds the on-size product", p1["probe_binds"] is True)
poff = next(p for p in v["products"] if p["subject"] == "XR_off")
chk("off-size product marked not on_size", poff["on_size"] is False)
chk("verdict mentions off-size", "off-size" in v["verdict"], v["verdict"])

# verdict: probe binds inside an OFF-size product -> cross-reactivity (the serious case)
v2 = SP.assay_verdict(prod3, [{"subject": "XR_off", "lo": 50, "hi": 70}])
chk("probe cross-binding detected", v2["n_probe_binding"] == 1, v2["n_probe_binding"])
chk("verdict flags cross-reactivity", "cross-react" in v2["verdict"], v2["verdict"])

# clean assay: single on-size product with probe inside
v3 = SP.assay_verdict([{"subject": "NM_t", "size": 150, "span": [100, 249]}],
                      [{"subject": "NM_t", "lo": 150, "hi": 170}])
chk("clean verdict: no off-size", v3["n_off_size"] == 0 and v3["n_probe_binding"] == 0
    and "no off-size" in v3["verdict"], v3["verdict"])

# probe hit on a DIFFERENT subject than the product does not count as binding
v4 = SP.assay_verdict([{"subject": "NM_t", "size": 150, "span": [100, 249]}],
                      [{"subject": "OTHER", "lo": 150, "hi": 170}])
chk("probe on other subject does not bind", v4["products"][0]["probe_binds"] is False)

# empty product set
v5 = SP.assay_verdict([], [])
chk("empty -> no products verdict", v5["n_products"] == 0 and "no predicted products" in v5["verdict"])

print("")
if fails:
    print("SPECIFICITY GATE FAILED:", ", ".join(fails)); sys.exit(1)
print("ALL SPECIFICITY ASSERTS PASS")
