"""Release gate for deterministic, non-mutating ranking-truth artifacts."""
import hashlib, json, os, sys, tempfile
from pathlib import Path
sys.path.insert(0, os.path.dirname(__file__))
from run_ranking_benchmark import main, OUTPUT_NAMES


ROOT = Path(__file__).resolve().parent
CANONICAL = ROOT / "benchmark"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


before = {name: digest(CANONICAL / name) for name in OUTPUT_NAMES}
with tempfile.TemporaryDirectory() as tmp:
    first = Path(tmp) / "first"
    second = Path(tmp) / "second"
    if main(first) != 0 or main(second) != 0:
        raise SystemExit("ranking benchmark acceptance gate failed")
    for name in OUTPUT_NAMES:
        left = first / name
        right = second / name
        if not left.is_file() or left.read_bytes() != right.read_bytes():
            raise SystemExit("ranking benchmark output is not byte deterministic: %s" % name)
        if (CANONICAL / name).read_bytes() != left.read_bytes():
            raise SystemExit("committed ranking benchmark artifact is stale: %s" % name)
    result = json.loads((first / "ranking_truth_results.json").read_text(encoding="utf-8"))
    if "runtime_seconds" in json.dumps(result, sort_keys=True):
        raise SystemExit("frozen benchmark results contain volatile runtime fields")

after = {name: digest(CANONICAL / name) for name in OUTPUT_NAMES}
if after != before:
    raise SystemExit("ranking benchmark regression mutated committed evidence")

print("PASS ranking benchmark outputs are byte deterministic")
print("PASS committed ranking benchmark evidence matches a fresh run")
print("PASS ranking benchmark regression leaves committed evidence unchanged")
