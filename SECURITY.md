# Security

## Supported deployment modes

OligoForge is safest as a local, single-user application. Public multi-user deployment requires `OLIGOFORGE_HOSTED=1`.

Hosted mode disables by default:

- server-side project and panel persistence;
- process-wide reaction-condition changes;
- local BLAST database paths supplied by clients.

The API applies request-size limits, security headers, sanitized validation responses, and generic public error messages while retaining full diagnostics in server logs.

## Operator responsibilities

- Set a real `OLIGOFORGE_EMAIL` for NCBI Entrez and keep `OLIGOFORGE_NCBI_KEY` in environment secrets.
- Put the service behind TLS and a trusted reverse proxy.
- Add authentication before enabling server persistence or shared mutable settings.
- Do not mount sensitive directories where a local BLAST process or application worker can read them.
- Restrict network egress and resource limits appropriate to computationally expensive endpoints.
- Rotate credentials if they appear in logs, uploaded projects, or support bundles.

## Reporting a vulnerability

Do not open a public issue containing credentials, private sequence data, filesystem paths, or exploit details. Contact the repository owner privately and include the affected version, endpoint, minimal reproduction, and impact.
