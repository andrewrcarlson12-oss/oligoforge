"""Stable, privacy-preserving API problem responses.

The legacy top-level ``error`` string is retained for older clients.  New code
should treat ``problem`` as authoritative: it distinguishes correctable input,
scientific no-solution outcomes, transient dependencies, capacity, and defects.
Request identifiers deliberately contain no host, user, or sequence material.
"""
from __future__ import annotations

from contextvars import ContextVar
import re
import secrets
from typing import Any, Iterable, Mapping, Optional

from fastapi.responses import JSONResponse


PROBLEM_SCHEMA = "oligoforge-problem/v1"
_REQUEST_ID: ContextVar[Optional[str]] = ContextVar("oligoforge_request_id", default=None)
_SAFE_REQUEST_ID = re.compile(r"^[A-Za-z0-9._-]{1,64}$")


def request_id(candidate: Optional[str] = None) -> str:
    """Return a safe caller correlation ID or create an opaque local one."""
    value = str(candidate or "").strip()
    if value and _SAFE_REQUEST_ID.fullmatch(value):
        return value
    return "ofreq_" + secrets.token_hex(12)


def bind_request_id(value: str):
    return _REQUEST_ID.set(value)


def reset_request_id(token) -> None:
    _REQUEST_ID.reset(token)


def current_request_id() -> Optional[str]:
    return _REQUEST_ID.get()


def problem_payload(message: str, *, code: str, category: str,
                    retryable: bool = False, recovery: Optional[Iterable[str]] = None,
                    stage: Optional[str] = None, retry_after_seconds: Optional[int] = None,
                    field_errors: Optional[Iterable[Mapping[str, Any]]] = None,
                    request_id_value: Optional[str] = None,
                    extra: Optional[Mapping[str, Any]] = None) -> dict[str, Any]:
    rid = request_id_value or current_request_id()
    problem = {
        "schema_version": PROBLEM_SCHEMA,
        "code": str(code),
        "category": str(category),
        "message": str(message),
        "stage": str(stage) if stage else None,
        "retryable": bool(retryable),
        "retry_after_seconds": (int(retry_after_seconds)
                                if retry_after_seconds is not None else None),
        "field_errors": [dict(x) for x in (field_errors or [])],
        "recovery": [str(x) for x in (recovery or []) if str(x).strip()],
        "request_id": rid,
    }
    payload: dict[str, Any] = {"error": str(message), "problem": problem}
    if rid:
        payload["request_id"] = rid
    if extra:
        payload.update(dict(extra))
    return payload


def problem_response(message: str, *, code: str, category: str,
                     status_code: int, retryable: bool = False,
                     recovery: Optional[Iterable[str]] = None,
                     stage: Optional[str] = None,
                     retry_after_seconds: Optional[int] = None,
                     field_errors: Optional[Iterable[Mapping[str, Any]]] = None,
                     extra: Optional[Mapping[str, Any]] = None) -> JSONResponse:
    headers = {}
    if retry_after_seconds is not None:
        headers["Retry-After"] = str(max(0, int(retry_after_seconds)))
    return JSONResponse(
        problem_payload(
            message, code=code, category=category, retryable=retryable,
            recovery=recovery, stage=stage, retry_after_seconds=retry_after_seconds,
            field_errors=field_errors, extra=extra,
        ),
        status_code=int(status_code), headers=headers,
    )

