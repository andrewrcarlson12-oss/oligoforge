# OligoForge 1.35.0 deployment

OligoForge is safest as a local, single-user service. The shipped container and Render blueprint can run a public web service, but the application has no built-in authentication, authorization, tenant isolation, durable job queue, or regulatory electronic-record controls. A public or sensitive-data deployment must add those controls outside the application.

## Shipped runtime

The repository `Dockerfile` is the deployment authority. It:

- builds from `python:3.12-slim`;
- installs the exact direct versions in `requirements.txt`;
- copies the repository into `/app`;
- runs as non-root user/UID 10001;
- binds Uvicorn to `0.0.0.0:${PORT:-8111}`;
- limits Uvicorn concurrency to 24 and keep-alive to 30 seconds; and
- disables proxy-header trust with `--no-proxy-headers`.

The image does not add TLS, authentication, a reverse proxy, a persistent disk, a shared queue, or an external database.

## Local container deployment

Build from the repository root:

```bash
docker build --pull -t oligoforge:1.35.0 .
```

Run locally and bind only to loopback on the host:

```bash
docker run --rm \
  --name oligoforge \
  --publish 127.0.0.1:8111:8111 \
  --env OLIGOFORGE_EMAIL=operator@example.org \
  oligoforge:1.35.0
```

Use a secret manager or a permission-restricted environment file for `OLIGOFORGE_NCBI_KEY`; avoid shell history, image layers, source files, logs, and request bodies. `OLIGOFORGE_EMAIL` is an NCBI client identity and should identify a monitored contact.

Check the service:

```bash
curl --fail --silent --show-error http://127.0.0.1:8111/healthz
curl --fail --silent --show-error http://127.0.0.1:8111/api/autodesign/limits
curl --fail --silent --show-error http://127.0.0.1:8111/openapi.json >/dev/null
```

`/healthz` reports readiness details including application version/commit, boot time, Primer3 availability, writable data directory, and whether NCBI credentials are configured. Do not publish its configuration details to untrusted users without review.

## Local persistent data

Set `OLIGOFORGE_DATA_PATH` to an explicit mounted directory when NCBI cache, projects, or panels must survive container replacement:

```bash
docker run --rm \
  --name oligoforge \
  --publish 127.0.0.1:8111:8111 \
  --env OLIGOFORGE_DATA_PATH=/var/lib/oligoforge \
  --mount type=bind,src=/srv/oligoforge-data,dst=/var/lib/oligoforge \
  oligoforge:1.35.0
```

The mounted host directory must be writable by container UID 10001. Back it up and test restore according to the organization's retention policy. Do not mount a broad home, credential, or system directory. Project and panel JSON storage is simple instance-local file storage, not a concurrent database.

Automatic-design jobs never persist to this directory. They remain process memory only and are lost on restart, redeploy, crash, or process replacement.

## Render blueprint

`render.yaml` defines one Docker web service named `oligoforge` with:

- the free plan and automatic deploys;
- `/healthz` as the health-check path;
- `OLIGOFORGE_HOSTED=1`;
- server project/panel storage disabled;
- shared reaction-condition mutation disabled; and
- dashboard-supplied `OLIGOFORGE_EMAIL` and `OLIGOFORGE_NCBI_KEY` secrets.

The blueprint declares no persistent disk. Do not rely on filesystem state or NCBI cache durability across deploys. Review automatic deployment policy for a controlled release; `autoDeploy: true` deploys repository changes rather than waiting for an OligoForge-specific approval gate.

The Docker command accepts the platform-provided `PORT`. Keep `/healthz` reachable by Render. If additional authentication is placed in front of all routes, exempt the platform health check only as narrowly as necessary and do not expose other APIs anonymously.

## Hosted-mode boundary

`OLIGOFORGE_HOSTED=1` is required for an internet-facing instance. By default it disables:

- server-side project and panel persistence;
- process-wide reaction-condition changes; and
- client-supplied local BLAST database paths.

Do not override `OLIGOFORGE_ALLOW_SERVER_STORAGE` or `OLIGOFORGE_ALLOW_SHARED_CONDITIONS` on a shared service until authentication, authorization, tenant isolation, concurrency behavior, audit, and recovery have been designed and tested.

Remote NCBI operations create network egress. Remote BLAST sends sequences to NCBI; the automatic-design workflow sends only the winning primer pair and only when the caller explicitly requests remote BLAST. Assurance CLI commands and Assurance calculations are offline and do not retrieve data.

## Automatic-design job deployment constraints

The job manager has one scientific worker and a bounded in-memory waiting queue. It is not a distributed job system.

- Job IDs are random capability-style identifiers; there is no job-list endpoint.
- Default queue capacity is 8, terminal retention is 1,800 seconds, primary timeout is 240 seconds, and optional BLAST timeout is 360 seconds.
- Submit and BLAST-retry routes accept `Idempotency-Key`.
- Cancellation is observed at stage boundaries. A native Primer3 or network/BLAST call already in progress may continue to drain before another scientific stage begins.
- A process restart causes subsequent polling to return “not found, expired, or lost on restart.”
- Job snapshots omit credentials, local database paths, fetched corpora, and sequence-like internal values.

Run one application process per service instance for this backend. Multiple processes or horizontally scaled instances create independent queues and capability namespaces; a polling request routed to a different process cannot find the job. Do not scale this implementation horizontally without sticky routing for the entire job lifetime or, preferably, replacing the backend with an authenticated shared durable queue and store.

## Environment variables

| Variable | Default | Deployment guidance |
|---|---:|---|
| `PORT` | `8111` in the Docker command | Platform listen port. |
| `OLIGOFORGE_HOSTED` | `0` | Set to `1` for any internet-facing deployment. |
| `OLIGOFORGE_ALLOW_SERVER_STORAGE` | on locally; off hosted | Leave off for public/shared service. |
| `OLIGOFORGE_ALLOW_SHARED_CONDITIONS` | on locally; off hosted | Leave off for public/shared service. |
| `OLIGOFORGE_DATA_PATH` | application/temp default | Mount an explicit least-privilege path when local state must persist. |
| `OLIGOFORGE_EMAIL` | unset | Required monitored identity for responsible NCBI use. |
| `OLIGOFORGE_NCBI_KEY` | unset | Optional secret; never bake into the image. |
| `OLIGOFORGE_NCBI_TIMEOUT` | 30 s | Clamped to 3–300 seconds. |
| `OLIGOFORGE_NCBI_CACHE` | `1` | Set `0` to disable response caching. |
| `OLIGOFORGE_NCBI_CACHE_TTL` | 604,800 s | Cache retention in seconds. |
| `OLIGOFORGE_NCBI_RETRIES` | `3` | Bounded to 1–5. |
| `OLIGOFORGE_JOB_QUEUE` | `8` | In-memory waiting capacity; minimum 1. |
| `OLIGOFORGE_JOB_TTL_SECONDS` | `1800` | Terminal in-memory retention; minimum 60 seconds in `app.py`. |
| `OLIGOFORGE_DESIGN_TIMEOUT_SECONDS` | `240` | Primary-stage deadline; minimum 30 seconds. |
| `OLIGOFORGE_BLAST_TIMEOUT_SECONDS` | `360` | Optional specificity-stage deadline; minimum 30 seconds. |
| `OLIGOFORGE_MAX_REQUEST_BYTES` | 5 MiB | Declared body-size ceiling. Enforce an equal or smaller proxy limit. |
| `OLIGOFORGE_LOG_LEVEL` | `INFO` | Avoid verbose logs containing sensitive operational context. |

## Reverse proxy and security controls

The shipped Uvicorn command uses `--no-proxy-headers`, so it deliberately does not trust forwarded client/protocol headers. Keep that boundary unless the exact proxy chain and trusted IPs are configured and tested. TLS can terminate at a trusted platform or reverse proxy; the application itself listens with plain HTTP.

Before exposing the service, supply at least:

1. TLS with certificate lifecycle management;
2. authentication and least-privilege authorization for every API and UI route;
3. request, concurrency, CPU, memory, and execution-rate limits appropriate to expensive scientific endpoints;
4. network egress policy for NCBI and denial of unintended internal/private destinations;
5. secret management, rotation, and log redaction;
6. access logging, monitoring, alerting, and incident response without logging sequence payloads or credentials;
7. dependency/container scanning and a software SBOM for the exact built image;
8. backup, restore, retention, and deletion for any mounted state;
9. tenant isolation if more than one trust domain is served; and
10. deployment-specific validation and change control before regulated reliance.

Security headers and a 5 MiB application body limit are defense-in-depth only. They are not a web application firewall, denial-of-service protection, identity system, or compliance control set.

## Release, upgrade, and rollback

For a controlled deployment:

1. Build from an approved source revision and record the source commit, image digest, base-image digest, dependency inventory, test results, and configuration.
2. Scan the exact image and review third-party licenses/notices; top-level pins do not enumerate all transitive/native components.
3. Deploy to a staging environment with the same hosted/storage/job settings and external services.
4. Verify `/healthz`, `/openapi.json`, representative local calculations, intended NCBI behavior, hosted restrictions, job submission/polling/cancellation, restart-loss messaging, logs, and resource limits.
5. Promote the immutable image digest, not a floating tag.
6. Retain the previous approved image/configuration and a compatible backup of mounted data.
7. During rollback, expect all in-memory jobs to disappear. Communicate that callers must resubmit rather than treating loss as scientific failure.

There is no automated data migration or distributed queue migration in the shipped configuration. Review schema compatibility and make a recoverable copy before changing any persisted project/panel data.

## Not supplied by this deployment

The container or Render blueprint does not implement Aegis, an Assurance Repair subsystem, FutureProof design, continuous sequence surveillance, scheduled retrieval, a durable notification service, a regulatory evidence repository, electronic signatures, or a validated audit trail.

## Implementation and official sources

Repository sources reviewed: `Dockerfile`, `render.yaml`, `app.py`, `oligoforge/jobs.py`, `oligoforge/ncbi.py`, `SECURITY.md`, and `API.md`.

Official external sources, all accessed 2026-07-15:

- Docker, building images: <https://docs.docker.com/get-started/docker-concepts/building-images/>
- Render, Blueprint YAML reference: <https://render.com/docs/blueprint-spec>
- Render, Docker deployment: <https://render.com/docs/docker>
- Uvicorn settings: <https://www.uvicorn.org/settings/>

