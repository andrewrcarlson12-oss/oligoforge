"""Factory reset (server half): /api/factory_reset deletes every saved panel and project, touches
only *.json inside the two managed dirs, leaves anything else alone, and is idempotent."""
import sys; sys.path.insert(0, ".")
import os
import app as APP
from fastapi.testclient import TestClient

c = TestClient(APP.app)

# seed saved panels + a project, plus a non-json file that must survive
open(os.path.join(APP.PANELS_DIR, "t1.json"), "w").write("[]")
open(os.path.join(APP.PANELS_DIR, "t2.json"), "w").write("[]")
open(os.path.join(APP.PROJECTS_DIR, "p1.json"), "w").write("{}")
keep = os.path.join(APP.PANELS_DIR, "keep.txt"); open(keep, "w").write("x")

r = c.post("/api/factory_reset", json={})
assert r.status_code == 200, r.status_code
d = r.json()
assert d["ok"] is True
assert d["panels"] >= 2 and d["projects"] >= 1, d
assert not [f for f in os.listdir(APP.PANELS_DIR) if f.endswith(".json")], "panels not cleared"
assert not [f for f in os.listdir(APP.PROJECTS_DIR) if f.endswith(".json")], "projects not cleared"
assert os.path.exists(keep), "a non-json file must be left untouched"
os.remove(keep)

# idempotent: running again on empty dirs returns 0 and no error
r2 = c.post("/api/factory_reset", json={})
assert r2.status_code == 200 and r2.json()["ok"] is True and r2.json()["panels"] == 0

print("ALL RESET ASSERTS PASS")
