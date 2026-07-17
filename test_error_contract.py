#!/usr/bin/env python3
"""Focused structured API-error, request-correlation, and diagnostics gates."""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from fastapi.testclient import TestClient

import app as APP
from oligoforge import api_errors as AE
from oligoforge import design_contract as DC


passed = failed = 0


def check(name, condition, detail=None):
    global passed, failed
    if condition:
        passed += 1
        print("PASS", name)
    else:
        failed += 1
        print("FAIL", name, detail if detail is not None else "")


def check_problem(name, response, *, status, code, category, request_id=None):
    body = response.json()
    problem = body.get("problem") or {}
    ok = (
        response.status_code == status and
        problem.get("schema_version") == AE.PROBLEM_SCHEMA and
        problem.get("code") == code and
        problem.get("category") == category and
        isinstance(problem.get("retryable"), bool) and
        isinstance(problem.get("recovery"), list) and
        bool(problem.get("request_id")) and
        body.get("request_id") == problem.get("request_id") == response.headers.get("x-request-id")
    )
    if request_id is not None:
        ok = ok and problem.get("request_id") == request_id
    check(name, ok, {"status": response.status_code, "headers": dict(response.headers), "body": body})
    return body


client = TestClient(APP.app)
TARGET = (
    "GGCCCGGATTCAGACTGATAGTGTAGTTATGATGCTCCGTGAGCTATACCCCGACCTCTGCTTTGAGATTGTGGCC"
    "ATGTCAACAACTGGGGACAAGATCTTGGATACAGCGCTTTCCAAGATTGGAGAGAAGAGTCTCTTCACCAAAGAGTTGGAAAATGCACTTGAAAGAA"
)

rid = "contract-test_request.42"
unknown_profile = client.post(
    "/api/design",
    json={"template": TARGET, "profile": "does_not_exist", "objective": "balanced"},
    headers={"X-Request-ID": rid},
)
body = check_problem(
    "unknown chemistry profile is a typed input problem",
    unknown_profile,
    status=422,
    code="input.unknown_profile",
    category="input",
    request_id=rid,
)
check("unknown-profile problem identifies the failed stage and recovery",
      body.get("problem", {}).get("stage") == "profile_resolution" and
      any("/api/profiles" in item for item in body.get("problem", {}).get("recovery", [])), body)

# A valid explicit profile previously crossed the validation boundary as a dict
# and was then (incorrectly) reused as a registry key, producing an internal 500.
original_design_from_sequences = APP._AD_design.design_from_sequences
try:
    APP._AD_design.design_from_sequences = lambda *_args, **_kwargs: {"error": "synthetic no solution"}
    explicit_profile = client.post(
        "/api/design",
        json={"template": TARGET, "profile": "idt_taqman", "objective": "balanced"},
        headers={"X-Request-ID": "explicit-profile-contract-test"},
    )
finally:
    APP._AD_design.design_from_sequences = original_design_from_sequences
check_problem(
    "valid explicit chemistry reaches design instead of raising an internal profile-key failure",
    explicit_profile,
    status=422,
    code="design.no_solution",
    category="scientific_no_solution",
    request_id="explicit-profile-contract-test",
)

secret = "SECRET_INPUT_VALUE_" + "X" * 100
invalid = client.post(
    "/api/qc",
    json={"role": {"invalid": secret}},
    headers={"X-Request-ID": "validation-contract-test"},
)
validation_body = check_problem(
    "request validation is a typed problem",
    invalid,
    status=422,
    code="input.validation_failed",
    category="input",
    request_id="validation-contract-test",
)
check("validation diagnostics list fields without echoing rejected values",
      bool(validation_body.get("problem", {}).get("field_errors")) and
      validation_body.get("details") == validation_body["problem"]["field_errors"] and
      secret not in str(validation_body), validation_body)

missing = client.get("/api/definitely-not-a-route", headers={"X-Request-ID": "route-contract-test"})
check_problem(
    "unknown API route is a typed request problem",
    missing,
    status=404,
    code="request.not_found",
    category="request",
    request_id="route-contract-test",
)

# Unsafe caller-controlled IDs are replaced with opaque IDs, never reflected.
unsafe = "../private path/APIKEY=SECRET"
unsafe_response = client.get("/api/definitely-not-a-route", headers={"X-Request-ID": unsafe})
unsafe_body = unsafe_response.json()
generated_rid = unsafe_body.get("problem", {}).get("request_id", "")
check("unsafe request ID is replaced and correlated consistently",
      generated_rid.startswith("ofreq_") and unsafe not in str(unsafe_body) and
      unsafe_response.headers.get("x-request-id") == generated_rid, unsafe_body)

diagnostics = client.get("/api/system/diagnostics", headers={"X-Request-ID": "diagnostics-contract-test"})
diag = diagnostics.json()
check("system diagnostics expose a stable non-sensitive readiness schema",
      diagnostics.status_code == 200 and
      diagnostics.headers.get("x-request-id") == "diagnostics-contract-test" and
      diag.get("schema_version") == "oligoforge-system-diagnostics/v1" and
      isinstance(diag.get("ok"), bool) and
      diag.get("core_design", {}).get("policy_id") == DC.POLICY_ID and
      diag.get("core_design", {}).get("contract_schema") == DC.CONTRACT_SCHEMA and
      isinstance(diag.get("capabilities"), dict) and
      isinstance(diag.get("limits"), dict), diag)
check("diagnostics make their privacy boundary explicit",
      "No sequences" in diag.get("privacy", "") and
      "credentials" in diag.get("privacy", "") and
      TARGET not in str(diag), diag.get("privacy"))

contract = DC.build_contract({"candidates": []}, workflow="api_contract_test", targets=[TARGET])
verified = client.post("/api/design-contracts/verify", json={"contract": contract},
                       headers={"X-Request-ID": "verify-contract-test"})
check("contract verification endpoint accepts an untampered contract",
      verified.status_code == 200 and verified.json().get("valid") is True and
      verified.headers.get("x-request-id") == "verify-contract-test", verified.json())
compared = client.post("/api/design-contracts/compare", json={"left": contract, "right": contract})
check("contract comparison endpoint reports equivalent reproduction",
      compared.status_code == 200 and compared.json().get("state") == "equivalent_reproduction",
      compared.json())

# Internal exceptions are logged server-side but emit only a typed, sanitized
# support envelope to the caller.
original_report_build = APP.RPT.build
try:
    APP.RPT.build = lambda *_args, **_kwargs: (_ for _ in ()).throw(
        RuntimeError("/tmp/private/APIKEY=SECRET_INTERNAL_DETAIL")
    )
    internal = client.post(
        "/api/report",
        json={"panel": []},
        headers={"X-Request-ID": "internal-contract-test"},
    )
finally:
    APP.RPT.build = original_report_build

internal_body = check_problem(
    "internal endpoint failure is typed and uses an error HTTP status",
    internal,
    status=500,
    code="internal.report_failed",
    category="internal",
    request_id="internal-contract-test",
)
check("internal problem omits exception detail and gives support recovery",
      "SECRET_INTERNAL_DETAIL" not in str(internal_body) and
      internal_body.get("problem", {}).get("stage") == "report" and
      any("diagnostic request ID" in item for item in internal_body["problem"].get("recovery", [])),
      internal_body)

print("ERROR_CONTRACT: %d passed / %d failed / %d total" % (passed, failed, passed + failed))
raise SystemExit(1 if failed else 0)
