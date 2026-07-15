#!/usr/bin/env bash
set -u
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"
PYTHON_BIN="${PYTHON:-python}"
PY_TIMEOUT="${OLIGOFORGE_TEST_TIMEOUT:-300}"
NODE_TIMEOUT="${OLIGOFORGE_NODE_TEST_TIMEOUT:-120}"
export OLIGOFORGE_EMAIL="${OLIGOFORGE_EMAIL:-ci@example.com}"
export PYTHONPATH="$ROOT${PYTHONPATH:+:$PYTHONPATH}"
mode="${1:-}"
fail=0; py_pass=0; py_fail=0; node_pass=0; node_fail=0
run_with_timeout(){
  local seconds="$1"; shift
  if command -v timeout >/dev/null 2>&1; then timeout "${seconds}s" "$@"; else "$@"; fi
}
if [[ "$mode" == "" || "$mode" == "--python" ]]; then
  echo "== Python suite =="
  shopt -s nullglob
  tests=(tests/test_*.py)
  for t in "${tests[@]}"; do
    name="$(basename "$t")"
    log="$(mktemp)"
    if PYTHONUNBUFFERED=1 run_with_timeout "$PY_TIMEOUT" "$PYTHON_BIN" -u "$t" >"$log" 2>&1; then
      echo "  PASS $name"; py_pass=$((py_pass+1))
    else
      rc=$?; echo "  FAIL $name (exit $rc)"; tail -n 12 "$log" | sed 's/^/      /'; py_fail=$((py_fail+1)); fail=1
    fi
    rm -f "$log"
  done
  echo "PYTHON: $py_pass passed / $py_fail failed / $((py_pass+py_fail)) total"
fi
if [[ "$mode" == "" || "$mode" == "--node" ]]; then
  echo "== Node UI harnesses =="
  if command -v node >/dev/null 2>&1; then
    harnesses=(tests/ui_*.js)
    for h in "${harnesses[@]}"; do
      name="$(basename "$h")"; log="$(mktemp)"
      if run_with_timeout "$NODE_TIMEOUT" node "$h" >"$log" 2>&1; then
        echo "  PASS $name"; node_pass=$((node_pass+1))
      else
        rc=$?; echo "  FAIL $name (exit $rc)"; tail -n 12 "$log" | sed 's/^/      /'; node_fail=$((node_fail+1)); fail=1
      fi
      rm -f "$log"
    done
  else
    echo "  SKIP node unavailable"
  fi
  echo "NODE:   $node_pass passed / $node_fail failed / $((node_pass+node_fail)) total"
fi
if [[ "$fail" -ne 0 ]]; then echo "CI FAILED"; exit 1; fi
echo "CI PASS: all requested suites green"
